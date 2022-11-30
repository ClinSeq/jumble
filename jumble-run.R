#!/usr/bin/env Rscript

jumble_version <- '0.2.0'

# Markus Mayrhofer 2022

# Dependencies ------------------------------------------------------------
{
    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(bamsignals))
    suppressPackageStartupMessages(library(Rsamtools))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(MASS))
    suppressPackageStartupMessages(library(VariantAnnotation))
    suppressPackageStartupMessages(library(PSCBS))
    suppressPackageStartupMessages(library(ggplot2)); theme_set(theme_bw())
    suppressPackageStartupMessages(library(patchwork))
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    suppressPackageStartupMessages(library(BSgenome))
    suppressPackageStartupMessages(library(Repitools))
}


# Options ------------------------------------------------------------
option_list <- list(
    make_option(c("-r", "--reference_file"), action = "store", type = "character",default = NULL,
                help = "reference file e.g. <design>.bed.references.RDS"),
    make_option(c("-b", "--input_bam"), action = "store", type = "character",default = NULL,
                help = "input bam or bam.counts.RDS file"),
    make_option(c("-v", "--snp_vcf"), action = "store", type = "character",default = NULL,
                help = "optional vcf file with SNPs, including any matched normal"),
    make_option(c("-o", "--output_dir"), action = "store", type = "character",default = '.',
                help = "directory for output")
)

opt <- parse_args(OptionParser(option_list = option_list))



# Reference file ------------------------------------------------------------

reference <- readRDS(opt$reference_file)
counts_template <- reference[c("target_bed_file","chromlength","ranges")]

wgs <- reference$target_bed_file=='wgs'

#save.image('ws.Rdata')


# Fragment counts ------------------------------------------------------------

# function for bam process
countsFromBam <- function(counts,bampath) {
    counts$date_count <- date()
    counts$input_bam_file <- bampath
    ranges <- counts$ranges
    counts[['count']] <- bamCount(bampath, ranges, paired.end="midpoint",
                                  mapq=20, filteredFlag=1024, verbose=F) # todo: 3840?
    counts[['count_short']] <- bamCount(bampath, ranges, paired.end="midpoint",
                                        mapq=20, filteredFlag=1024, tlenFilter=c(0,150), verbose=F)
    return(counts)
}

# parse this sample
input <- opt$input_bam
if (str_detect(input,'.RDS$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- readRDS(input)
} else if (str_detect(input,'.[bB][aA][mM]$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- countsFromBam(counts_template,input)
} else stop("No input file?")



# Tables of bins ------------------------------------------------------------


targets <- reference$targets

name <- str_remove(opt$input_bam,'.*/')
name <- str_remove(name,'\\.counts.RDS')

clinbarcode <- str_remove(name, "[_-]nodups.bam")

targets$sample <- clinbarcode


# Mark tiled
targets[,left:=c(Inf,abs(diff(mid)))]
targets[,right:=c(abs(diff(mid)),Inf)]
targets[,farthest:=left][right>left,farthest:=right]
targets[,is_tiled:=F][farthest<250,is_tiled:=T]
targets[,left:=NULL][,right:=NULL][,farthest:=NULL]
tiled_genes <- as.data.table(sort(table(targets[is_tiled==T]$gene),decreasing = T))

# Add counts
targets[,count:=counts$count]
alltargetcount <- targets[is_target==T & bin %in% reference$keep]$count

# same, short
targets[,count_short:=counts$count_short]



# For SNPs if present
targets[,snps:=0]
targets[,allele_count_correction:=0]


target_ranges <- makeGRangesFromDataFrame(targets[is_target==T])
ranges <- makeGRangesFromDataFrame(targets)



# Modify backbone ----------------------------------------------------------
set.seed(25) # <------------------ To be reproducible.
max_backbone_in_gene <- 100 #  <--- applies to some

if (!wgs) {
    genes <- c('ATM','BRCA1','BRCA2','PTEN','RB1')
    targets[,is_backbone:=chromosome %in% 1:22 & !gene %in% genes]
    for (g in genes) {
        ix=targets[gene==g & chromosome %in% as.character(1:22),.I]
        n <- length(ix)
        if (n>max_backbone_in_gene) ix=ix[order(rnorm(n))][1:max_backbone_in_gene]
        targets[ix,is_backbone:=T]
    }
}



# SNP allele ratio ------------------------------------------------------------
#save.image('ws.Rdata')

peakx <- function(data) {
    if(length(data)==1) return(data)
    maxx <- NA
    try({
        d <- density(data[!is.na(data)])
        maxy <- max(d$y)
        maxx <- d$x[d$y==maxy][1]
    }, silent=T)
    return(maxx)
}

snp_allele_ratio <- FALSE
input <- opt$snp_vcf
if (!is.null(input)) {
    snp_allele_ratio <- TRUE
    
    if (is.null(reference$snp_rlm_model)) 
        warning('SNP data in sample but not in reference object. Raw SNP allele ratio will be plotted but not used.')
    
    if (!str_detect(input,'.[vV][cC][fF]$') & !str_detect(input,'.[vV][cC][fF].[gG][zZ]$')) 
        stop("SNP vcf file appears incorrect")
    
    vcf=readVcf(input)
    
    # If there is more than one sample in the VCF, which one to use?
    ix <- 1 # base assumption: first is this sample
    names <- colnames(vcf)
    if (length(names)>1) {
        # but if one name in the VCF fits the sample name better, use it:
        check <- which(str_detect(name,names))
        if (length(check)==1) ix <- check
        
        # Keep only this sample
        vcf <- vcf[,ix]
        
    }
    
    
    
    # remove SNPs that do not have 2 alleles
    alleles <- as.data.table(table(as.data.table(alt(vcf))$group))$N
    vcf <- vcf[alleles==1 & sapply(geno(vcf)$AD,length)==2]
    
    name <- colnames(vcf)
    
    g <- geno(vcf)
    
    # Prepare table of SNPs for this sample
    snp_table <- as.data.table(rowRanges(vcf))[,1:3][,chromosome:=seqnames][,c(4,2,3)]
    snp_table <- cbind(
        data.table(
            sample=name,
            id=names(vcf)),
        snp_table
    )
    
    # Get alleles
    snp_table$ref_allele <- as.character(ref(vcf))
    alt_allele <- as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
    snp_table$n_alt_alleles <- sapply(alt_allele, length)
    snp_table$alt_allele <- sapply(alt_allele, "[[", 1)
    
    # Get SNP type
    # from <- rep('C/G',length(vcf)); from[snp_table$ref_allele %in% c('A','T')] <- 'A/T'
    # to <- rep('C/G',length(vcf)); to[snp_table$alt_allele %in% c('A','T')] <- 'A/T'
    snp_table[!str_detect(ref_allele,'^[ACGT]$'),ref_allele:='other']
    snp_table[!str_detect(alt_allele,'^[ACGT]$'),alt_allele:='other']
    snp_table[,type:=paste0(ref_allele,'>',alt_allele)]
    snp_table[str_length(ref_allele)!=1 | str_length(alt_allele)!=1]$type <- 'other'
    
    # Get read counts
    snp_table$AD <- sapply(g$AD, "[[", 2)
    snp_table$RD <- sapply(g$AD, "[[", 1)
    snp_table$DP=snp_table$AD+snp_table$RD
    snp_table[,logDP:=log2(DP)]
    
    # Compute raw allele ratio
    raw_allele_ratio <- unname(round(snp_table$AD/snp_table$DP,4))
    raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
    snp_table$allele_ratio <- raw_allele_ratio
    
    # SNP id with details
    snp_table[,snp:=paste(id,ref_allele,alt_allele)]
    
    # Overlap with bins, drop SNPs not on targets
    overlap <- findOverlaps(makeGRangesFromDataFrame(snp_table),ranges)
    snp_table[queryHits(overlap),bin:=subjectHits(overlap)]
    snp_table <- snp_table[bin %in% targets[is_target==T]$bin]
    
    
    
    # SNP correction ----------------------------------------------------------
    
    if (!is.null(reference$snp_rlm_model)) {
        
        # Compute GC content for SNPs
        gc_ranges <- makeGRangesFromDataFrame(unique(snp_table[,.(chromosome,start,end)]))
        seqlevelsStyle(gc_ranges) <- "UCSC"
        gc_ranges <- gc_ranges[width(gc_ranges)==1]
        start(gc_ranges) <- start(gc_ranges)-1
        end(gc_ranges) <- end(gc_ranges)+1
        gc_ranges$gc3 <- gcContentCalc(gc_ranges, organism=Hsapiens)
        start(gc_ranges) <- start(gc_ranges)-1
        end(gc_ranges) <- end(gc_ranges)+1
        gc_ranges$gc5 <- gcContentCalc(gc_ranges , organism=Hsapiens)
        start(gc_ranges) <- start(gc_ranges)-48
        end(gc_ranges) <- end(gc_ranges)+48
        gc_ranges$gc101 <- gcContentCalc(gc_ranges , organism=Hsapiens)
        start(gc_ranges) <- start(gc_ranges)+50
        end(gc_ranges) <- end(gc_ranges)-50
        seqlevelsStyle(gc_ranges) <- "NCBI"
        gc_table <- as.data.table(gc_ranges)[,-c('strand','width')][,chromosome:=seqnames][,-'seqnames']
        snp_table <- merge(snp_table,gc_table,all=T,by=c('chromosome','start','end'))
        
        
        # Compute features for ref allele bias
        snp_table[,hom_alt:=F][chromosome %in% 1:22 & allele_ratio > .95,hom_alt:=T]
        
        t <- c('C>T','G>A')
        snp_table[,peak_ld:=peakx(log2(DP[type %in% t & hom_alt==T & chromosome %in% 1:22])),by=c('sample')]
        
        snp_table[,homC_T:=peakx(log2(DP[hom_alt==T & type=='C>T']))-peak_ld,by=c('sample')]
        snp_table[,homT_C:=peakx(log2(DP[hom_alt==T & type=='T>C']))-peak_ld,by=c('sample')]
        snp_table[,homA_G:=peakx(log2(DP[hom_alt==T & type=='A>G']))-peak_ld,by=c('sample')]
        snp_table[,homG_A:=peakx(log2(DP[hom_alt==T & type=='G>A']))-peak_ld,by=c('sample')]
        
        # Compute general model based bias estimates
        snp_rlm_model <- reference$snp_rlm_model
        snp_table[,bias_estimate:=.04] # default 
        for (ref in unique(snp_table$ref_allele)) 
            for (alt in unique(snp_table$alt_allele)) 
                if(ref!=alt & ref != 'other') {
                    cat(ref,'>',alt,'\n')
                    thismodel <- snp_rlm_model[[paste0(ref,'_',alt)]]
                    
                    ix <- snp_table$ref_allele==ref & snp_table$alt_allele==alt
                    cat(sum(ix),'total \n')
                    predicted_bias <- predict(thismodel,snp_table[ix])
                    cat(length(predicted_bias),'predicted \n')
                    #hist(predicted_bias)
                    snp_table[ix,bias_estimate:=predicted_bias]
                }
        
        snp_table[,AD_corrected:=AD*(2^bias_estimate)]
        snp_table[,DP_corrected:=RD+AD_corrected]
        snp_table[,allele_ratio_corrected:=AD_corrected/DP_corrected]
        
        # Compute residual bias estimates where available
        snp_coeff_table <- reference$snp_coeff_table
        ix <- match(snp_table$snp,snp_coeff_table$snp)
        snp_table[,residual_bias_estimate:=0]
        snp_table[,residual_bias_estimate:=
                      snp_coeff_table$constant[ix]+
                      snp_coeff_table$logDP[ix]*logDP+
                      snp_coeff_table$peak_ld[ix]*peak_ld+
                      snp_coeff_table$homC_T[ix]*homC_T+
                      snp_coeff_table$homT_C[ix]*homT_C+
                      snp_coeff_table$homA_G[ix]*homA_G+
                      snp_coeff_table$homG_A[ix]*homG_A]
        snp_table[is.na(residual_bias_estimate),residual_bias_estimate:=0]
        
        # Correct for residual bias
        snp_table[,AD_corrected2:=AD_corrected*(2^residual_bias_estimate)]
        snp_table[,DP_corrected2:=RD+AD_corrected2]
        snp_table[,allele_ratio_corrected2:=AD_corrected2/DP_corrected2]
        
        # Compute the correction factor and cap it at 2/3 to 3/2. Note that only the model based, first correction factor is used for that.
        snp_table[,correct_factor:=2^(bias_estimate+residual_bias_estimate)]
        snp_table[correct_factor>1.5,correct_factor:=1.5]
        snp_table[correct_factor<2/3,correct_factor:=2/3]
        
        # Compute the correction number
        snp_table[,correct_number:=(correct_factor-1)*AD]
        snp_table[,correct_number_check:=AD_corrected2-AD] # should be same
        
        # Add correction number to target table
        snps_by_bin <- snp_table[,.N,by=bin]
        targets[snps_by_bin$bin,snps:=snps_by_bin$N]
        correction_by_bin <- snp_table[,max(correct_number),by=bin]
        targets[correction_by_bin$bin,allele_count_correction:=correction_by_bin$V1]
        
    }
    
}




# LogR and genotype correction ------------------------------------------------------------




min1 <- function(data) {
    data[data<1] <- 1
    data[is.na(data)] <- 1
    return(data)
}

# Basic logR, targets (with correction for SNPs)
targets[,rawLR:=log2(min1(count+allele_count_correction))]
targets[,rawLR_short:=log2(min1(count_short+allele_count_correction*(count_short/count)))] # allelic correction scaled for short-fragments

# median correct to backbone
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by='is_target']
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by='is_target']

# correct by reference median
targets[reference$keep,rawLR:=rawLR-reference$median]
targets[reference$keep,rawLR_short:=rawLR_short-reference$median_short]



# Remove outliers 1 ------------------------------------------------------------

alltargets <- copy(targets) # to get the full set back later

# keep only bins "ok" in this reference set
targets <- targets[bin %in% reference$keep]


# PCA outliers based on score SD
set.seed(25)
pca <- as.data.table(prcomp(reference$targets_ref[,-1],center = F,scale. = F)$x)
targets[,keep:=T]

for (pc in colnames(pca)) {
    fact <- ifelse(pc %in% c('PC1','PC2'),4,4)
    sd <- sd(pca[[pc]])
    targets[pca[[pc]] < -sd*fact, keep:=F]
    targets[pca[[pc]] > sd*fact, keep:=F]
}
targets <- targets[keep==T]
targets[,keep:=NULL]




# PCA v1 ------------------------------------------------------------

ix <- targets[is_target==T]$bin # the ontarget

set.seed(25)
tpca <- as.data.table(prcomp(reference$targets_ref[bin %in% ix,-1],center = F,scale. = F)$x)
tpca_short <- as.data.table(prcomp(reference$targets_ref_short[bin %in% ix,-1],center = F,scale. = F)$x)

if (!wgs) if (any(targets$is_target==F)) {
    ix <- targets[is_target==F]$bin # the offtarget
    bgpca <- as.data.table(prcomp(reference$targets_ref[bin %in% ix,-1],center = F,scale. = F)$x)
    bgpca_short <- as.data.table(prcomp(reference$targets_ref_short[bin %in% ix,-1],center = F,scale. = F)$x)
}



# Reference data correction v1 ------------------------------------------------------------

# correct using reference PCA
jcorrect <- function(temp,train_ix=NULL) {
    
    if (is.logical(train_ix[1])) train_ix <- which(train_ix)
    if (is.null(train_ix)) train_ix <- which(rep(TRUE,nrow(temp)))
    
    set.seed(25)
    if (length(train_ix)>10e3) {
        new_ix <- train_ix[order(rnorm(n=length(train_ix)))][1:10e3]
        train_ix <- new_ix
    }
    
    pcs <- ncol(temp)-2
    bins <- nrow(temp)
    
    if (pcs>=3 & bins<50e3) {
        # loess_temp <- loess(lr ~ PC1+PC2+PC3, data = temp,
        #                     subset = train_ix,
        #                     family="symmetric", control = loess.control(surface = "direct"))
        # 
        # temp[,lr:=lr-predict(loess_temp,temp)]
        loess_temp <- rlm(lr ~ PC1+PC2+PC3, data=temp,
                          subset = train_ix)
        temp[,lr:=lr-predict(loess_temp,temp)]
    }
    
    
    for (i in 1:(pcs)) {
        temp$thispc=temp[[paste0('PC',i)]]
        loess_temp <- rlm(lr ~ thispc, data=temp,
                          subset = train_ix)
        temp[,lr:=lr-predict(loess_temp,temp)]
    }
    
    # apply GC correct as linear function, where there is gc, if many bins
    if (bins > 50e3 & any(!is.na(temp$gc))) {
        ix <- which(!is.na(temp$gc))
        loess_temp <- rlm(lr ~ gc, data=temp[ix],
                          subset = intersect(train_ix,ix))
        temp[ix,lr:=lr-predict(loess_temp,temp[ix])]
        # or loess if fewer
    } else if (any(!is.na(temp$gc))) {
        ix <- which(!is.na(temp$gc))
        loess_temp=loess(lr ~ gc, data = temp[ix],
                         subset = intersect(train_ix,ix),
                         family="symmetric", control = loess.control(surface = "direct"))
        temp[ix,lr:=lr-predict(loess_temp,temp[ix])]
    }
    
    return(temp$lr)
}


targets[,log2:=rawLR]
targets[,log2_short:=rawLR_short]

# standard
ix <- targets$is_target
temp <- cbind(data.table(
    lr=targets[ix]$rawLR),
    gc=targets[ix]$gc,
    tpca)
targets[ix,log2:=jcorrect(temp,targets[ix]$is_backbone)]

# short
temp <- cbind(data.table(
    lr=targets[ix]$rawLR_short),
    gc=targets[ix]$gc,
    tpca_short)
targets[ix,log2_short:=jcorrect(temp,targets[ix]$is_backbone)]

if (!wgs) if (any(targets$is_target==F)) {
    ix <- targets$is_target
    # standard bg
    temp <- cbind(data.table(
        lr=targets[!ix]$rawLR),
        gc=targets[!ix]$gc,
        bgpca)
    targets[!ix,log2:=jcorrect(temp,targets[!ix]$is_backbone)]
    
    # bg short
    temp <- cbind(data.table(
        lr=targets[!ix]$rawLR_short),
        gc=targets[!ix]$gc,
        bgpca_short)
    targets[!ix,log2_short:=jcorrect(temp,targets[!ix]$is_backbone)]
    
}







# Alternative correction ------------------------------------------------------------

# These are performed only for "targets".

if (!wgs) {
    
    # Simple reference + GC
    ix <- which(!is.na(targets$gc))
    loess_temp=loess(rawLR ~ gc, data = targets[ix],
                     family="symmetric", control = loess.control(surface = "direct"))
    targets[ix,log2_simple:=rawLR-predict(loess_temp,targets[ix])]
    
    
    # Pca projection + GC
    ix <- targets[is_target==T]$bin # the ontarget
    set.seed(25)
    pca <- prcomp(t(reference$targets_ref[bin %in% ix,-1]),center = F,scale. = F)
    npcs <- ceiling(ncol(pca$x)/2)
    query_x <- t(targets[bin %in% ix,.(rawLR)]) %*% pca$rotation[,1:npcs]
    projection <- t(query_x %*% t(pca$rotation[,1:npcs]))
    targets[bin %in% ix,log2_pca:=rawLR - projection]
    loess_temp=loess(log2_pca ~ gc, data = targets[bin %in% ix],
                     family="symmetric", control = loess.control(surface = "direct"))
    targets[bin %in% ix,log2_pca:=log2_pca-predict(loess_temp,targets[bin %in% ix])]
    
    # Jumble correct incl GC but no subselect
    
    ix <- targets$is_target
    temp <- cbind(data.table(
        lr=targets[ix]$rawLR),
        gc=targets[ix]$gc,
        tpca)
    targets[ix,log2_nosub:=jcorrect(temp)]
}


# Remove outliers 2 ------------------------------------------------------------

# deviation <- function(vector) {
#     m <- runmed(vector,k = 7)
#     d <- abs(vector-m)
#     return(d)
# }
# 
# targets[,dev:=deviation(log2),by=chromosome]
#targets <- targets[dev<1]

targets[log2 < -4, log2:=-4]
targets[log2 > 7, log2:=7]
#targets[log2x < -4, log2x:=-4]
#targets[log2x > 7, log2x:=7]



# Adjust X to background ------------------------------------------------------------

if (!wgs) try( 
    {
        # X-chromosome correction factor to targeted bins, based on background bins
        temp <- targets[chromosome=='X' & !is.na(log2)]
        temp[,bg:=as.numeric(NA)][is_target==F,bg:=2^log2]
        temp[,bg_median:=runmed(bg,51,na.action = 'na.omit')]
        temp[,tg:=as.numeric(NA)][is_target==T & is_tiled==F & gene=='',tg:=2^log2]
        temp[,tg_median:=runmed(tg,51,na.action = 'na.omit')]
        temp[,dif:=bg_median-tg_median]
        x_correct <- median(temp$dif,na.rm = T)
        targets[chromosome=='X' & is_target==T, log2:=log2+x_correct]
    }, silent=T
)

# Segmentation ------------------------------------------------------------


getsegs <- function(targets, logratio) {
    
    alpha <- .01
    if (wgs) { 
        alpha <- 1e-5
        #if (median(targets$end-targets$start) >= 5000) alpha <- 1e-5
    }
    
    segments <- segmentByCBS(y=logratio,avg='median',
                             chromosome=targets$chromosome,
                             alpha = alpha,undo=1)
    segments <- as.data.table(segments)[!is.na(chromosome),-1][!is.na(mean)]
    segments[,start_pos:=targets$start[ceiling(start)]]
    segments[,end_pos:=targets$end[floor(end)]]
    segments[,genes:='']
    
    targets[,segment:=as.numeric(NA)]
    
    for (i in 1:nrow(segments)) {
        ix <- ceiling(segments[i]$start):floor(segments[i]$end)
        targets[ix,segment:=i]
        segments[i,mean:=median(logratio[ix],na.rm = T)] # <----------------- here, value per segment is set
        genes <- paste0(unique(targets[ix][gene!='']$gene),collapse = ',')
        genes <- unique(strsplit(genes,',')[[1]])
        genes <- genes[!genes %in% c('','Background')]
        if (length(genes)>0) segments[i]$genes <- paste(genes,collapse = ', ')
        
        # adjust segment start and end to bin number, rather than bin order
        ix <- ceiling(segments[i]$start)
        segments[i,start:=targets[ix]$bin]
        ix <- floor(segments[i]$end)
        segments[i,end:=targets[ix]$bin]
    }
    
    
    
    
    return(segments)
}



# Do the segmentation
targets[,chromosome:=str_replace(chromosome,'Y','24')][,chromosome:=str_replace(chromosome,'X','23')][,chromosome:=as.numeric(chromosome)]
segments <- getsegs((targets), targets$log2)
targets[,chromosome:=as.character(chromosome)][chromosome=='23',chromosome:='X'][chromosome=='24',chromosome:='Y']
segments[,chromosome:=as.character(chromosome)][chromosome=='23',chromosome:='X'][chromosome=='24',chromosome:='Y']


# Gene+segment table ------------------------------------------------------------

segments_temp <- segments[,.(segment=paste(1:.N),type='segment',chromosome,start=start_pos,end=end_pos,
                             length=end_pos-start_pos,
                             bins=nbrOfLoci,genes,mean)]

genes <- targets[is_target==T,.(segment='',type='gene',chromosome,start,end,length=NA,bins=0,
                                genes=gene,log2,
                                mean=0)]


for (i in 1:nrow(segments)) {
    ix <- ceiling(segments[i]$start):floor(segments[i]$end)
    genes[ix]$segment <- i
}
genes <- genes[genes!='']
genes[,segment:=paste(unique(segment),collapse = ','),by=genes]
genes[,start:=min(start),by=genes]
genes[,end:=max(end),by=genes]
genes[,length:=end-start]
genes[,bins:=.N,by=genes]
suppressWarnings(genes[,mean:=round(median(log2,na.rm=T),3),by=genes])

genes[,log2:=NULL]

suppressWarnings(
    segments_genes <- rbind(segments_temp,unique(genes))[order(as.numeric(chromosome),start)]
)




# Table output ------------------------------------------------------------

targets <- merge(alltargets,targets,by=colnames(alltargets),all=T)[order(bin)]

# The combined segments and genes table (skipped for now)
#fwrite(x = segments_genes,file = paste0(opt$output_dir,'/',clinbarcode,'.segments.csv'))

# Jumble targets and background
saveRDS(targets,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble.RDS'))
saveRDS(snp_table,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble_snps.RDS'))

# for compatibility with CNVkit.
# cnr:  chromosome      start   end     gene    depth   log2    weight ()
cnr <- targets[!is.na(log2),.(chromosome=as.character(chromosome),start,end,gene,
                              depth=round(count/width*200,3),log2,weight=1,
                              gc,count)][gene=='',gene:='-']

#fwrite(x = cnr,file = paste0(opt$output_dir,'/',clinbarcode,'.cnr'),sep = '\t')

# cns:  chromosome      start   end     gene    log2    depth   probes  weight
cns <- segments_genes[type=='segment',.(chromosome,start,end,gene=genes,log2=mean,depth=mean,probes=bins)]
#fwrite(x = cns,file = paste0(opt$output_dir,'/',clinbarcode,'.cns'),sep = '\t')


# DNAcopy segment file:
# ID    chrom   loc.start       loc.end num.mark        seg.mean        C
seg <- segments[,.(ID=name,chrom=chromosome,loc.start=start_pos,loc.end=end_pos,num.mark=nbrOfLoci,seg.mean=mean,C=NA)]
seg[,chrom:=str_replace(chrom,'Y','24')][,chrom:=str_replace(chrom,'X','23')][,chrom:=as.numeric(chrom)]
#fwrite(x = seg,file = paste0(opt$output_dir,'/',clinbarcode,'_dnacopy.seg'),sep = '\t')


# Count file output ------------------------------------------------------------
# (not overwrite, not if input was counts.RDS)
if (!file.exists(paste0(opt$output_dir,'/',clinbarcode,'.*counts.RDS')))
    if (!str_detect(opt$input_bam,'counts.RDS$'))
        saveRDS(counts,paste0(opt$output_dir,'/',clinbarcode,'.counts.RDS'))


# Save workspace? ------------------------------------------------------------
#save.image(paste0(opt$output_dir,'/',clinbarcode,'.jumble_workspace.Rdata'))



# QC metrics ------------------------------------------------------------
#save.image(paste0(opt$output_dir,'/',clinbarcode,'.jumble_workspace.Rdata'))


mapd <- function(data) {
    return(median(abs(diff(data)),na.rm = T))
} 

noise <- function(data) {
    m <- mapd(data)
    f <- 2^m-1
    return(round(100*f,1))
}

if (!wgs) {
    
    # tiled gene bias
    tiled_genes <- as.data.table(sort(table(targets[is_tiled==T]$gene),decreasing = T))
    include_genes <- tiled_genes[N>50]$V1
    targets[,chr_median:=as.numeric(NA)][,chr_median:=peakx(log2[is_tiled==F & is_target==T]),by=chromosome]
    targets[,tiled_median:=as.numeric(NA)][,tiled_median:=peakx(log2[gene %in% include_genes & is_tiled==T]),by=gene]
    targets[,tiled_bias:=tiled_median-chr_median]
    tiled <- unique(targets[!is.na(tiled_bias),.(sample,gene,chr_median,tiled_median,tiled_bias)])
    tiled_bias <- round(100*2^median(abs(tiled$tiled_bias))-100,1)
    
    # waviness 11 to 51 (non-tiled) target bins
    temp <- targets[is_tiled==F & is_target==T][!is.na(log2)]
    temp[,peak11:=frollapply(log2, 11, peakx)]
    temp[,peak51:=frollapply(log2, 51, peakx)]
    waviness <- round(100*2^median(abs(temp$peak11-temp$peak51),na.rm=T)-100,1)
    
    
    stats <- paste0('Fragments per target: ',
                    paste(round(quantile(alltargetcount,c(.025,.975))),collapse = '-'),
                    ', Noise: ',
                    noise(targets[is_target==T]$log2),'%/', noise(targets[is_target==F]$log2),'%',
                    ', Bias: ', tiled_bias,'%',
                    ', Waviness: ', waviness,'%'
    )
}



if (wgs) stats <- paste0('Fragments per target: ',
                         paste(round(quantile(alltargetcount,c(.025,.975))),collapse = '-'),
                         ', Noise: ',
                         noise(targets$log2),'%'
)




# Plot ------------------------------------------------------------



if (F) {
    
    
    
    p <- NULL
    targets[,smooth_log2:=runmed(log2,k=7),by=chromosome]
    ylims <- c(.4,max(2,max(2^targets$smooth_log2)))
    
    size <- 1; if (wgs) size <- 2
    
    if (wgs) {
        label_genes <- c('BRCA2','PTEN','RB1')
        targets[,label:=as.character(NA)]
        for (g in label_genes) {
            targets[str_detect(gene,paste0('^',g,'$')),`label`:=g]
            targets[str_detect(gene,paste0(',',g,'$')),`label`:=g]
            targets[str_detect(gene,paste0('^',g,',')),`label`:=g]
            targets[str_detect(gene,paste0(',',g,',')),`label`:=g]
        }
        
    } else {
        #label_genes <- c('AR','ATM','BRCA2','PTEN','RB1','ERG','CDK12','TMPRSS2')
        label_genes <- c('BRCA2','PTEN')
        targets[,label:=as.character(NA)]
        targets[is_tiled==T,label:='other_tiled']
        for (g in label_genes) {
            targets[str_detect(gene,paste0('^',g,'$')),`label`:=g]
            targets[str_detect(gene,paste0(',',g,'$')),`label`:=g]
            targets[str_detect(gene,paste0('^',g,',')),`label`:=g]
            targets[str_detect(gene,paste0(',',g,',')),`label`:=g]
        }
    }
    
    #targets <- targets[is_target==T]
    
    ## Grid ------------------------------------------------------------
    
    if (snp_allele_ratio) { 
        
        # defaults to using raw allele ratio:
        snp_table[,allele_ratio_use:=allele_ratio] 
        # but if there is a corrected allele ratio, use it:
        if (!is.null(snp_table$allele_ratio_corrected2)) snp_table[,allele_ratio_use:=allele_ratio_corrected2]
        
        snp_table <- snp_table[type!='other'][allele_ratio_use < .99][allele_ratio_use > .01]
        snp_table <- snp_table[DP > median(DP)/3][DP < median(DP)*3]
        
        targets[,allele_ratio:=as.double(NA)][match(snp_table$bin,bin),allele_ratio:=snp_table$allele_ratio_use]
        targets[,maf:=abs(allele_ratio-.5)+.5]
        targets[!is.na(maf),maf:=runmed(maf,7)]
        # snp (grid) smooth-to-allele-ratio plot
        p$grid <- ggplot(targets) + xlim(c(.2,1.8)) + ylim(c(.5,1)) + xlab('Corrected depth (smooth)') + ylab('Major allele ratio (smooth)') +
            geom_point(data=targets[,.(smooth_log2,maf)],aes(x=2^smooth_log2,y=maf),col='lightgrey',alpha=.2) +
            geom_point(aes(x=2^smooth_log2,y=maf),fill='#60606090',col='#20202090',shape=21) +
            geom_point(data=targets[label!=''],aes(x=2^smooth_log2,y=maf,fill=label),shape=21,col='#00000050',size=size) +
            facet_wrap(facets = vars(factor(chromosome,levels=unique(chromosome),ordered=T)),ncol = 8) +
            theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8))
        # snp (all) smooth-to-allele-ratio plot
        temp <- targets[!is.na(label),median(log2),by=label]
        p$nogrid <- ggplot(targets) + xlim(c(0.2,1.8)) + ylim(c(.5,1)) + xlab('Corrected depth (smooth)') + ylab('Major allele ratio (smooth)') +
            geom_point(data=targets[,.(smooth_log2,maf)],aes(x=2^smooth_log2,y=maf),fill='#60606050',col='#20202050',shape=21) +
            geom_point(data=targets[label!=''],aes(x=2^smooth_log2,y=maf,fill=label),shape=21,col='#00000050',size=size) +
            geom_point(data=temp,mapping=aes(x=2^V1,y=1,fill=label),size=2,shape=25,show.legend=F)
    }
    
    # By pos ------------------------------------------------------------
    
    # chroms object by genomic pos
    chroms <- data.table(chromosome=names(reference$chromlength),length=reference$chromlength)
    chroms[,start:=as.double(0)] 
    chroms[,stop:=as.double(length)] 
    chroms[,mid:=as.double(round(length/2))] 
    for (i in 2:nrow(chroms)) {
        chroms[i,start:=chroms$stop[i-1]]
        chroms[i,stop:=chroms$stop[i-1]+length]
        chroms[i,mid:=chroms$stop[i-1]+round(length/2)]
    }    
    # fix positions by genome
    targets[,gpos:=mid]
    for (chr in unique(targets$chromosome)[-1]) targets[chromosome==chr,gpos:=gpos+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
    segments[,gstart:=as.double(start_pos)][,gstop:=as.double(end_pos)]
    for (chr in unique(targets$chromosome)[-1]) {
        segments[chromosome==chr,gstart:=gstart+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
        segments[chromosome==chr,gstop:=gstop+sum(chroms[1:(which(chromosome==chr)-1)]$length)]
    }
    # logR by pos + segments (2nd left)
    p$pos_log2 <- ggplot(targets) + xlab('Genomic position') + ylab('Corrected depth') +
        geom_point(data=targets[is.na(label)],mapping = aes(x=gpos,y=2^log2),fill='#60606070',col='#20202070',shape=21,size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=2^log2,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=gstart,xend=gstop,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    
    # depth by pos 
    limits <- c(1e2,1e4)
    p$pos_rawdepth <- ggplot(targets) + xlab('Genomic position') + ylab('Read count') +
        geom_point(data=targets[is_target==T],mapping = aes(x=gpos,y=count),fill='#60606050',col='#20202050',size=1,shape=21) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=count,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=limits) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line())
    
    if (snp_allele_ratio) {
        # allele ratio by pos
        p$pos_alleleratio <- ggplot(targets) + xlab('Genomic position') + ylab('Allele ratio') +
            geom_point(data=targets[is.na(label)],mapping = aes(x=gpos,y=allele_ratio),fill='#60606080',col='#20202080',shape=21,size=1) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=allele_ratio,fill=label),shape=21,col='#00000050',size=1) +
            scale_fill_hue() + ylim(0:1) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line()) 
    }
    
    # By order ------------------------------------------------------------
    
    # chroms object by order
    chroms=data.table(chromosome=unique(targets$chromosome),
                      start=0,
                      end=0,
                      mid=0)
    for (chr in chroms$chromosome) {
        chroms$start[chr==chroms$chromosome]=targets[chromosome==chr,min(bin)]
        chroms$end[chr==chroms$chromosome]=targets[chromosome==chr,max(bin)]
        chroms$mid[chr==chroms$chromosome]=targets[chromosome==chr,mean(bin)]
    }
    # logR by order
    p$order_log2 <- ggplot(targets) + xlab('Order of genomic position') + ylab('Corrected depth') +
        geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2),fill='#60606070',col='#20202070',size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc
    p$gc_log2 <- ggplot(targets) + xlab('Target GC content') + ylab('Corrected depth') + xlim(c(.2,.78)) +
        geom_point(data=targets,mapping = aes(x=gc,y=2^log2,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # fill='#60606040'
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    
    ##### short frags
    p$order_log2x <- ggplot(targets) + xlab('Order of genomic position') + ylab('Corrected depth (<150)') +
        geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2_short),fill='#60606070',col='#20202070',size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2_short,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc, short
    p$gc_log2x <- ggplot(targets) + xlab('Target GC content') + ylab('Corrected depth (<150)') + xlim(c(.2,.78)) +
        geom_point(data=targets,mapping = aes(x=gc,y=2^log2_short,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # fill='#60606040'
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    
    ##### testing alt. correction:
    p$order_log2_simple <- ggplot(targets) + xlab('Order of genomic position') + ylab('Simple corr.') +
        geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2_simple),fill='#60606070',col='#20202070',size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2_simple,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc
    p$gc_log2_simple <- ggplot(targets) + xlab('Target GC content') + ylab('Simple corr.') + xlim(c(.2,.78)) +
        geom_point(data=targets,mapping = aes(x=gc,y=2^log2_simple,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # fill='#60606040'
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    ##
    p$order_log2_pca <- ggplot(targets) + xlab('Order of genomic position') + ylab('PCA corr.') +
        geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2_pca),fill='#60606070',col='#20202070',size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2_pca,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc
    p$gc_log2_pca <- ggplot(targets) + xlab('Target GC content') + ylab('PCA corr.') + xlim(c(.2,.78)) +
        geom_point(data=targets,mapping = aes(x=gc,y=2^log2_pca,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # fill='#60606040'
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    ##
    p$order_log2_nosub <- ggplot(targets) + xlab('Order of genomic position') + ylab('Jumble, nosub') +
        geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2_nosub),fill='#60606070',col='#20202070',size=1) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2_nosub,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=ylims) +
        geom_segment(data=segments,col='green',size=1,
                     mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # logR by gc
    p$gc_log2_nosub <- ggplot(targets) + xlab('Target GC content') + ylab('Jumble, nosub') + xlim(c(.2,.78)) +
        geom_point(data=targets,mapping = aes(x=gc,y=2^log2_nosub,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # fill='#60606040'
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        scale_fill_hue() + scale_y_log10(limits=ylims) 
    ###
    
    m <- targets[is_target==T,median(count)]
    if (snp_allele_ratio) {
        # allele ratio by order 
        p$order_alleleratio <- ggplot(targets) + xlab('Order of genomic position') + ylab('Allele ratio') +
            geom_point(data=targets[is.na(label)],mapping = aes(x=bin,y=allele_ratio),fill='#60606050',col='#20202050',shape=21,size=1) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=allele_ratio,fill=label),shape=21,col='#00000050',size=size) +
            scale_fill_hue() + ylim(0:1) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line()) 
        # allele ratio by depth
        p$depth_alleleratio <- ggplot(targets) + xlab('Depth') + ylab('Allele ratio') +
            #geom_point(data=targets[is.na(label)],mapping = aes(x=count,y=allele_ratio),fill='#60606070',col='#20202070',shape=21,size=1) +
            geom_point(data=targets,mapping = aes(x=count,y=allele_ratio,fill=label),
                       shape=21,col='#00000050',size=1,alpha=.7,show.legend = F) +
            facet_wrap(facets = vars(label),ncol = 2) +
            theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
            scale_fill_hue() + ylim(0:1) + scale_x_log10(limits=c(m/3,m*3))
    }
    
    # Depth ------------------------------------------------------------
    
    # depth by order 
    p$order_rawdepth <- ggplot(targets) + xlab('Order of genomic position') + ylab('Depth') +
        geom_point(data=targets[is_target==T],mapping = aes(x=bin,y=count),fill='#60606050',col='#20202050',size=1,shape=21) +
        geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=count,fill=label),shape=21,col='#00000050',size=size) +
        scale_fill_hue() + scale_y_log10(limits=limits) +
        scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                           expand = c(.01,.01),labels = chroms$chromosome) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_line(),
              panel.grid.minor.x = element_line(color = 'black'),
              axis.line = element_line(),
              axis.ticks = element_line()) 
    # depth by GC 
    p$gc_rawdepth <- ggplot(targets) + xlab('Target GC content') + ylab('Depth') + xlim(c(.2,.78)) +
        geom_point(data=targets[is_target==T],mapping = aes(x=gc,y=count,fill=label),
                   col='#20202040',shape=21,size=1,alpha=.7,show.legend = F) + # 
        facet_wrap(facets = vars(label),ncol = 2) +
        theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
        # geom_smooth(data=targets[!is.na(label)],
        #             mapping = aes(x=gc,y=count,col=label),size=.5,se=F,show.legend = F,method = 'loess') +
        scale_fill_hue() + scale_y_log10(limits=limits)
    
    
    for (i in 1:length(p)) p[[i]] <- p[[i]] + guides(fill=guide_legend(override.aes=list(shape=21,size=3)))
    
    
    pa <- plot_annotation(
        title = paste(clinbarcode,'         ',stats),
        caption = paste('Jumble',jumble_version,'on',format(Sys.time(), "%a %b %e %Y, %H:%M"))
    )
    
    # Layout ------------------------------------------------------------
    
    if (snp_allele_ratio & !wgs) {
        
        layout <-  "ABBBB
                    CDDDD
                    EFFFF
                    GHHHH
                    IJJJJ
                    KLLLL"
        fig <-
            p$gc_rawdepth+p$order_rawdepth+
            p$gc_log2+p$order_log2+
            p$gc_log2_nosub+p$order_log2_nosub+
            p$gc_log2_simple+p$order_log2_simple+
            p$gc_log2_pca+p$order_log2_pca+
            p$depth_alleleratio+p$order_alleleratio+
            plot_layout(design = layout,guides = 'collect')
    }
    
    # if (snp_allele_ratio & !wgs) {
    #     
    #     layout <-  "ABBBB
    #             CDDDD
    #             EFFFF
    #             GGGGG
    #             HHHHH
    #             IJJJJ
    #             IJJJJ"
    #     fig <- 
    #         p$gc_rawdepth+p$order_rawdepth+
    #         p$gc_log2+p$order_log2+
    #         p$depth_alleleratio+p$order_alleleratio+
    #         p$pos_log2+
    #         p$pos_alleleratio+
    #         p$nogrid+p$grid+
    #         plot_layout(design = layout,guides = 'collect')
    # }
    if (!snp_allele_ratio & !wgs) {
        
        layout <-  "ABBBB
                CDDDD
                EFFFF
                GGGGG
                "
        fig <-
            p$gc_rawdepth+p$order_rawdepth+
            p$gc_log2+p$order_log2+
            p$gc_log2x+p$order_log2x+
            p$pos_log2+
            plot_layout(design = layout,guides = 'collect')
    }
    
    if (wgs) {
        
        layout <-  "ABBBB
                CDDDD
                "
        fig <- 
            p$gc_rawdepth+p$pos_rawdepth+
            p$gc_log2+p$pos_log2+
            plot_layout(design = layout,guides = 'collect')
    }
    
    png(file = paste0(opt$output_dir,'/',clinbarcode,'.png'),width = 1800,height=1400,res=100)
    print(fig+pa)
    
    
    
    
    
    # Close image ------------------------------------------------------------
    
    
    dev.off()
    
}


