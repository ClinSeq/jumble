#!/usr/bin/env Rscript

jumble_version <- '0.3'

# Markus Mayrhofer 2022-2023

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

this_is_wgs <- reference$target_bed_file=='wgs'



# save.image('ws.Rdata')
# stop()


# Target template ----------

{
    targets <- reference$target_template
    
    ## Mark tiled ---------
    targets[,left:=c(Inf,abs(diff(mid)))]
    targets[,right:=c(abs(diff(mid)),Inf)]
    targets[,farthest:=left][right>left,farthest:=right]
    targets[,is_tiled:=F][farthest<250,is_tiled:=T]
    targets[,left:=NULL][,right:=NULL][,farthest:=NULL]
    tiled_genes <- as.data.table(sort(table(targets[is_tiled==T]$gene),decreasing = T))
    targets[is_tiled==T,type:='tiled']
    
}

## Gene annotation -----------

{
    allgenes <- reference$allgenes
    allexons <- reference$allexons
    cancergenes_clinseq <- reference$cancergenes_clinseq
    
    
    cancergenes <- allgenes[`Gene stable ID` %in% cancergenes_clinseq$ensembl_gene_id_version]
    
    missing <- cancergenes_clinseq[!hugo_symbol %in% cancergenes$`Gene name` & 
                                       !ensembl_gene_id_version %in% cancergenes$`Gene stable ID`]
    cancergenes <- cancergenes[,.(ensembl_id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                                  start=`Gene start (bp)`,end=`Gene end (bp)`,
                                  symbol=`Gene name`)]
    cancergenes <- unique(cancergenes)
    cancergenes[,type:='A']
    cancergenes[ensembl_id %in% cancergenes_clinseq[ANNOT=='ONCO']$ensembl_gene_id_version,type:='O']
    cancergenes[ensembl_id %in% cancergenes_clinseq[ANNOT=='TSG']$ensembl_gene_id_version,type:='T']
    
    cancerexons <- allexons[`Gene stable ID` %in% cancergenes$ensembl_id,.(ensembl_id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                                                                           exon=`Exon rank in transcript`,
                                                                           start=`Exon region start (bp)`,end=`Exon region end (bp)`,
                                                                           symbol=`Gene name`)]
    
    
    allgenes <- allgenes[,.(ensembl_id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                            start=`Gene start (bp)`,end=`Gene end (bp)`,
                            symbol=`Gene name`)]
    allgenes <- unique(allgenes)
    
    
    # Add AR enhancer to gene/exon table
    ar_enh <- data.table(ensembl_id='enh_AR',chromosome='X',start=66100404,end=66160987,symbol='enh_AR',type='O')
    cancergenes <- rbind(cancergenes,ar_enh)
    
    ar_enh <- data.table(ensembl_id='enh_AR',chromosome='X',exon=1,start=66100404,end=66160987,symbol='enh_AR')
    cancerexons <- rbind(cancerexons,ar_enh)
    
    
    
    # Make ranges objects
    generanges <- makeGRangesFromDataFrame(cancergenes)
    exonranges <- makeGRangesFromDataFrame(cancerexons)
    binranges <- makeGRangesFromDataFrame(targets)
    
    
    # Gene overlap table
    gene_overlap <- as.data.table(findOverlaps(binranges,generanges))
    # Add symbol to overlap table
    gene_overlap[,symbol:=cancergenes[subjectHits]$symbol]
    
    for (i in unique(gene_overlap$queryHits)) {
        targets[i,gene:=paste0(gene_overlap[queryHits==i]$symbol,collapse = ',')]
    }
    
    
    # Exons overlap table
    gene_overlap <- as.data.table(findOverlaps(binranges,exonranges))
    # Add symbol to overlap table
    gene_overlap[,symbol:=cancerexons[subjectHits]$symbol]
    
    for (i in unique(gene_overlap$queryHits)) {
        targets[i,type:='exonic']
    }
    
    targets[is_target==F,type:='background']
    
}

## Backbone definition ------------------------------------------------------------

set.seed(25) # <------------------ To be reproducible.
max_backbone_in_gene <- 10 #  <--- applies to some

if (!this_is_wgs) {
    genes <- c('ATM','BRCA2','PTEN','RB1')
    targets[,is_backbone:=chromosome %in% 1:22 & !gene %in% genes]
    for (g in genes) {
        genebins=targets[gene==g & chromosome %in% as.character(1:22)]$bin
        n <- length(genebins)
        if (n>max_backbone_in_gene) genebins=genebins[order(rnorm(n))][1:max_backbone_in_gene]
        targets[bin %in% genebins,is_backbone:=T]
    }
}


## Template complete ------

target_template <- targets
#table(targets$type,targets$is_target)


# Manage reference set ----------


## SNP allele bias model ----------

if (!is.null(reference$snp_table)) {
    
    snp_table <- reference$snp_table
    
    train_table <- snp_table[allele_ratio>.25]
    
    # Check numbers per type
    n_per_type <- as.data.table(table(train_table$type))[order(N,decreasing = T)]
    
    
    # Compute features
    train_table[,hom_alt:=F][chromosome %in% 1:22 & allele_ratio > .95,hom_alt:=T]
    snp_table[,hom_alt:=F][chromosome %in% 1:22 & allele_ratio > .95,hom_alt:=T]
    train_table[,het:=F][allele_ratio > .35 & allele_ratio < .6,het:=T]
    
    # Depth peaks by sample of log2 (DP of all het C>T and G>A, nonX)
    t <- c('C>T','G>A')
    train_table[,peak_ld:=peakx(log2(DP[type %in% t & hom_alt==T & chromosome %in% 1:22])),by=c('sample')]
    snp_table[,peak_ld:=peakx(log2(DP[type %in% t & hom_alt==T & chromosome %in% 1:22])),by=c('sample')]
    
    # SNP type specific homozygous peak depth by sample, as log2 (this)-log2(sample_median)
    c('C>T','G>A','A>G','T>C')
    train_table[,homC_T:=peakx(log2(DP[hom_alt==T & type=='C>T']))-peak_ld,by=c('sample')]
    train_table[,homT_C:=peakx(log2(DP[hom_alt==T & type=='T>C']))-peak_ld,by=c('sample')]
    train_table[,homA_G:=peakx(log2(DP[hom_alt==T & type=='A>G']))-peak_ld,by=c('sample')]
    train_table[,homG_A:=peakx(log2(DP[hom_alt==T & type=='G>A']))-peak_ld,by=c('sample')]
    
    snp_table[,homC_T:=peakx(log2(DP[hom_alt==T & type=='C>T']))-peak_ld,by=c('sample')]
    snp_table[,homT_C:=peakx(log2(DP[hom_alt==T & type=='T>C']))-peak_ld,by=c('sample')]
    snp_table[,homA_G:=peakx(log2(DP[hom_alt==T & type=='A>G']))-peak_ld,by=c('sample')]
    snp_table[,homG_A:=peakx(log2(DP[hom_alt==T & type=='G>A']))-peak_ld,by=c('sample')]
    
    # Look at how they correlate
    temp <- unique(train_table[,.(peak_ld,homC_T,homT_C,homA_G,homG_A)])
    #pairs(temp)
    
    # Heterozygous SNPs table
    het_table <- train_table[het==T]
    #plot(density(het_table$allele_ratio))
    
    
    # Individual reference bias (het, near-normal) as log2
    het_table[,bias:=round(log2(RD/AD),4)]
    #plot(density(het_table$bias))
    
    # Allele bias model, general
    
    
    # build test model of bias from selected features, using RLM
    biasmodel <- rlm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=het_table)
    #summary(biasmodel)
    
    # same using lm
    lm_model <- lm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=het_table)
    #summary(lm_model)
    
    # But we really need separate models per variant type
    snp_rlm_model <- NULL
    for (ref in unique(het_table$ref_allele))
        for (alt in unique(het_table$alt_allele))
            if(ref!=alt & ref!='other') {
                data <- het_table[ref_allele==ref & alt_allele==alt]
                if (alt=='other') data <- het_table[ref_allele==ref]
                #if (ref=='other') data <- het_table[alt_allele==alt]
                biasmodel <- rlm(bias ~ logDP+gc3+gc5+gc101+peak_ld+homC_T+homT_C+homA_G+homG_A, data=data)
                cat(ref,'>',alt,':\n')
                print(summary(biasmodel))
                snp_rlm_model[[paste0(ref,'_',alt)]] <- biasmodel
            }
    
    # compute bias estimates
    snp_table[,bias_estimate:=peakx(het_table$bias)] # default
    het_table[,bias_estimate:=peakx(het_table$bias)] # default
    for (ref in unique(snp_table$ref_allele))
        for (alt in unique(snp_table$alt_allele))
            if(ref!=alt & ref != 'other') {
                cat(ref,'>',alt,'\n')
                thismodel <- snp_rlm_model[[paste0(ref,'_',alt)]]
                
                # in het snp table
                ix <- snp_table$ref_allele==ref & snp_table$alt_allele==alt
                cat(sum(ix),'total \n')
                predicted_bias <- predict(thismodel,snp_table[ix])
                cat(length(predicted_bias),'predicted \n')
                hist(predicted_bias)
                snp_table[ix,bias_estimate:=predicted_bias]
                
                # in snp table
                ix <- het_table$ref_allele==ref & het_table$alt_allele==alt
                cat(sum(ix),'total \n')
                predicted_bias <- predict(thismodel,het_table[ix])
                cat(length(predicted_bias),'predicted \n')
                hist(predicted_bias)
                het_table[ix,bias_estimate:=predicted_bias]
            }
    
    snp_table[,AD_corrected:=AD*(2^bias_estimate)]
    snp_table[,DP_corrected:=RD+AD_corrected]
    snp_table[,allele_ratio_corrected:=AD_corrected/DP_corrected]
    
    het_table[,AD_corrected:=AD*(2^bias_estimate)]
    het_table[,DP_corrected:=RD+AD_corrected]
    het_table[,allele_ratio_corrected:=AD_corrected/DP_corrected]
    
    # Allele bias model, common SNPs
    
    # How many heterozygous in reference of each SNP?
    het_table[,het_in_ref:=.N,by=c('id','ref_allele','alt_allele')]
    
    hist(het_table$het_in_ref[1:100000],100)
    
    het_table[,residual_bias:=log2((1-allele_ratio_corrected)/allele_ratio_corrected)]
    
    
    het_table[,constant:=0]
    het_table[het_in_ref>=5,constant:=mean(residual_bias),by=snp]
    #het_table[het_in_ref>=20,constant:=peakx(residual_bias),by=snp]
    
    
    snp_coeff_table <- unique(het_table[het_in_ref>=5,.(snp=paste(id,ref_allele,alt_allele),het_in_ref,constant,logDP=0,peak_ld=0,homC_T=0,homT_C=0,homA_G=0,homG_A=0)])
    
    for (i in which(snp_coeff_table$het_in_ref>20)) {
        data <- het_table[snp==snp_coeff_table$snp[i]]
        model <- rlm(residual_bias ~ logDP+peak_ld+homC_T+homT_C+homA_G+homG_A, data=data)
        coeff <- model$coefficients
        snp_coeff_table[i,constant:=coeff['(Intercept)']]
        snp_coeff_table[i,logDP:=coeff['logDP']]
        snp_coeff_table[i,peak_ld:=coeff['peak_ld']]
        snp_coeff_table[i,homC_T:=coeff['homC_T']]
        snp_coeff_table[i,homT_C:=coeff['homT_C']]
        snp_coeff_table[i,homA_G:=coeff['homA_G']]
        snp_coeff_table[i,homG_A:=coeff['homG_A']]
    }
    
    # Now predict the residual bias for all SNPs in training data (for sanity)
    ix <- match(het_table$snp,snp_coeff_table$snp)
    het_table[,residual_bias_estimate:=0]
    het_table[,residual_bias_estimate:=
                  snp_coeff_table$constant[ix]+
                  snp_coeff_table$logDP[ix]*logDP+
                  snp_coeff_table$peak_ld[ix]*peak_ld+
                  snp_coeff_table$homC_T[ix]*homC_T+
                  snp_coeff_table$homT_C[ix]*homT_C+
                  snp_coeff_table$homA_G[ix]*homA_G+
                  snp_coeff_table$homG_A[ix]*homG_A]
    het_table[is.na(residual_bias_estimate),residual_bias_estimate:=0]
    
    # And in full data for use
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
    het_table[,AD_corrected2:=AD_corrected*(2^residual_bias_estimate)]
    het_table[,DP_corrected2:=RD+AD_corrected2]
    het_table[,allele_ratio_corrected2:=AD_corrected2/DP_corrected2]
    
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
    
}



## Iterate reference samples ------------------------------------------------------------

allcounts <- reference$allcounts
targetlist <- NULL
for (i in 1:length(allcounts)) {
    
    counts <- allcounts[[i]]
    targets <- copy(target_template)
    name <- basename(counts$input_bam_file)
    targets$sample <- name
    
    
    # Add counts
    targets[,count:=counts$count]
    targets[,count_short:=counts$count_short]
    
    targets[,snps:=0]
    targets[,allele_count_correction:=0]
    
    # Add SNP info if present
    if (exists('snp_table')) {
        
        # Match vcf names in SNP table to this BAM file name
        this_table <- snp_table[str_detect(name,sample)]
        # Proceed if one matching name
        if (nrow(this_table)>0) if (length(unique(this_table$sample))==1) {
            snps_by_bin <- this_table[,.N,by=bin]
            targets[snps_by_bin$bin,snps:=snps_by_bin$N]
            correction_by_bin <- this_table[,max(correct_number),by=bin]
            targets[correction_by_bin$bin,allele_count_correction:=correction_by_bin$V1]
        }
        
    }
    
    targetlist[[i]] <- targets
}

targets <- rbindlist(targetlist)
rm(targetlist)





## Remove worst bins ------------------------------------------------------------

# Low coverage threshold
threshold <- median(targets$count) * 0.05
keep_targets <- targets[,median(count),by=bin][V1 > threshold]
targets <- targets[bin %in% keep_targets$bin]

# High coverage threshold
threshold <- median(targets$count) / 0.05
keep_targets <- targets[,median(count),by=bin][V1 < threshold]
targets <- targets[bin %in% keep_targets$bin]

# Low mappability removed
targets <- targets[map>0.6]


min1 <- function(data) {
    data[data<1] <- 1
    data[is.na(data)] <- 1
    return(data)
}

# Quantify variability for "blacklisting" high-variability regions.
if (this_is_wgs) {
    # sdev along each sample
    targets[,rollsd:= frollapply(log2(min1(count)), 10, sd, align = 'center',na.rm=T), by = sample]
    # median over samples, for each bin
    targets[,rollsd_median:= median(rollsd), by = bin]
    # standard deviation of that
    sdev <- sd(targets$rollsd_median,na.rm = T)
    # drop bins with more noise
    targets <- targets[!rollsd_median > sdev*4]
}


## SNP and Median correct ------------------------------------------------------------

# Basic logR, targets (with correction for SNPs if available)
targets[,rawLR:=log2(min1(count+allele_count_correction))]
targets[,rawLR_short:=log2(min1(count_short+allele_count_correction*(count_short/count)))] # allelic correction scaled for short-fragments

# Median correct, separated on sample and target/background status.
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by=c('sample','is_target')]
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by=c('sample','is_target')]


## X-Y chromosome correct ------------------------------------------------------------

# Double X values in samples where their median implies male
targets[,xmedian:=median(rawLR[chromosome=='X']),by=sample]
targets[,male:=2^xmedian < .75] # assign gender

targets[,nonPA:=chromosome %in% c('X') & end>2.70e6 & start<154.93e6] # hard coded PA

if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
    targets[chromosome=='X' & male==T & nonPA,rawLR:=rawLR+1]
    targets[chromosome=='X' & male==T & nonPA,rawLR_short:=rawLR_short+1]
}

# Double Y values where their median implies male
targets[,ymedian:=median(rawLR[chromosome=='Y']),by=sample] # median by sample
targets[,male:=2^ymedian > .25] # assign gender based on Y

if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
    targets[male==T & chromosome=='Y' & end<28.79e6,rawLR:=rawLR+1] # hard coded PA
    targets[male==T & chromosome=='Y' & end<28.79e6,rawLR_short:=rawLR_short+1]
}

# Y values set to NA where X median implied female
targets[chromosome=='Y' & male==F,rawLR:=NA]
targets[chromosome=='Y' & male==F,rawLR_short:=NA]


## Bin median correct ------------------------------------------------------------
targets[,refmedian:=median(rawLR,na.rm=T),by=bin][,rawLR:=rawLR-refmedian]
targets[,refmedian_short:=median(rawLR_short,na.rm=T),by=bin][,rawLR_short:=rawLR_short-refmedian_short]

bins_for_mediancorrect <- unique(targets[,.(bin,refmedian,refmedian_short)])

## Impute missing ----------------------------------------------------------

# If any missing, replace with random value near 0 (which is already the median)
targets[is.na(rawLR),rawLR:=rnorm(n = .N,mean = 0,sd = .1)]
targets[is.na(rawLR_short),rawLR_short:=rnorm(n = .N,mean = 0,sd = .1)]


## Matrix form ------------------------------------------------------------
mat <- dcast(data = targets[,.(bin,sample,rawLR)],formula = bin ~ sample, value.var = 'rawLR')
mat_short <- dcast(data = targets[,.(bin,sample,rawLR_short)],formula = bin ~ sample, value.var = 'rawLR_short')



## PCA 1: outliers ------------------------------------------------------------

set.seed(25)
pca <- as.data.table(prcomp(mat[,-1],center = F,scale. = F)$x)
pca$keep <- T

for (pc in colnames(pca)) {
    fact <- ifelse(pc %in% c('PC1','PC2'),4,4)
    sd <- sd(pca[[pc]])
    pca[pca[[pc]] < -sd*fact, keep:=F]
    pca[pca[[pc]] > sd*fact, keep:=F]
}

# keep certain bins in matrices and for query sample
mat <- mat[pca$keep==T]
mat_short <- mat_short[pca$keep==T]

bins_outliers_removed <- mat$bin

## PCA 2: latent features ------------------------------------------------------------

targetbins <- target_template[is_target==T]$bin # the ontarget

set.seed(25)
tpca <- as.data.table(prcomp(mat[bin %in% targetbins,-1],center = F,scale. = F)$x)
tpca_short <- as.data.table(prcomp(mat_short[bin %in% targetbins,-1],center = F,scale. = F)$x)

if (!this_is_wgs) if (any(targets$is_target==F)) {
    backgroundbins <- target_template[is_target==F]$bin # the offtarget
    bgpca <- as.data.table(prcomp(mat[bin %in% backgroundbins,-1],center = F,scale. = F)$x)
    bgpca_short <- as.data.table(prcomp(mat_short[bin %in% backgroundbins,-1],center = F,scale. = F)$x)
}





# Manage query sample ------------------------------------------------------------

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

## Parse this sample ---------------
input <- opt$input_bam
if (str_detect(input,'.RDS$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- readRDS(input)
} else if (str_detect(input,'.[bB][aA][mM]$')) {
    if (!file.exists(input)) stop(paste("Cannot find",input))
    counts <- countsFromBam(counts_template,input)
} else stop("No input file?")



## Table of bins ------------------------------------------------------------


# "targets" will now be the query sample.

targets <- target_template

name <- str_remove(opt$input_bam,'.*/')
name <- str_remove(name,'\\.counts.RDS')

clinbarcode <- str_remove(name, "[_-]nodups.bam")

targets$sample <- clinbarcode

# Add counts
targets[,count:=counts$count]

# same, short
targets[,count_short:=counts$count_short]

# for stats
alltargetcount <- targets[is_target==T]$count

# for retaining all bins in output
alltargets <- copy(targets)

# For SNPs if present
targets[,snps:=0]
targets[,allele_count_correction:=0]


target_ranges <- makeGRangesFromDataFrame(targets[is_target==T])
ranges <- makeGRangesFromDataFrame(targets)




## SNP allele ratio ------------------------------------------------------------
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
    
    
    
    ### SNP correction ----------------------------------------------------------
    
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




## LogR and genotype correction ------------------------------------------------------------


targets <- targets[bin %in% bins_for_mediancorrect$bin]


min1 <- function(data) {
    data[data<1] <- 1
    data[is.na(data)] <- 1
    return(data)
}

# Basic logR, targets (with correction for SNPs)
targets[,rawLR:=log2(min1(count+allele_count_correction))]
targets[,rawLR_short:=log2(min1(count_short+allele_count_correction*(count_short/count)))] # allelic correction scaled for short-fragments

# median correct to backbone, separately for targets/background
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by='is_target']
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by='is_target']

# correct by reference median
targets[,rawLR:=rawLR-bins_for_mediancorrect$refmedian]
targets[,rawLR_short:=rawLR_short-bins_for_mediancorrect$refmedian_short]




## Reference data correction ------------------------------------------------------------

targets <- targets[bin %in% bins_outliers_removed]

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

if (!this_is_wgs) if (any(targets$is_target==F)) {
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







## Alternative correction ------------------------------------------------------------

# These are performed only for "targets".

if (!this_is_wgs) {
    
    # Simple reference (done), adding GC:
    ix <- which(!is.na(targets$gc))
    loess_temp=loess(rawLR ~ gc, data = targets[ix],
                     family="symmetric", control = loess.control(surface = "direct"))
    targets[ix,log2_simple:=rawLR-predict(loess_temp,targets[ix])]
    
    
    # # Pca projection + GC
    # ix <- targets[is_target==T]$bin # the ontarget
    # set.seed(25)
    # pca <- prcomp(t(tmat[bin %in% ix,-1]),center = F,scale. = F)
    # npcs <- ceiling(ncol(pca$x)/2)
    # query_x <- t(targets[bin %in% ix,.(rawLR)]) %*% pca$rotation[,1:npcs]
    # projection <- t(query_x %*% t(pca$rotation[,1:npcs]))
    # targets[bin %in% ix,log2_pca:=rawLR - projection]
    # loess_temp=loess(log2_pca ~ gc, data = targets[bin %in% ix],
    #                  family="symmetric", control = loess.control(surface = "direct"))
    # targets[bin %in% ix,log2_pca:=log2_pca-predict(loess_temp,targets[bin %in% ix])]
    
    # Jumble correct incl GC but no subselect
    
    # ix <- targets$is_target
    # temp <- cbind(data.table(
    #     lr=targets[ix]$rawLR),
    #     gc=targets[ix]$gc,
    #     tpca)
    # targets[ix,log2_nosub:=jcorrect(temp)]
}


## Set min/max ------------------------------------------------------------

targets[log2 < -4, log2:=-4]
targets[log2 > 7, log2:=7]
targets[log2_short < -4, log2_short:=-4]
targets[log2_short > 7, log2_short:=7]



## Adjust X to background ------------------------------------------------------------

if (!this_is_wgs) try( 
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

## Segmentation ------------------------------------------------------------


getsegs <- function(targets, logratio) {
    
    alpha <- .01
    if (this_is_wgs) { 
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
    segments[,relevance:='']
    
    targets[,segment:=as.numeric(NA)]
    
    for (i in 1:nrow(segments)) {
        ix <- ceiling(segments[i]$start):floor(segments[i]$end)
        targets[ix,segment:=i]
        segments[i,mean:=median(logratio[ix],na.rm = T)] # <----------------- here, value per segment is set
        # genes <- paste0(unique(targets[ix][gene!='']$gene),collapse = ',')
        # genes <- unique(strsplit(genes,',')[[1]])
        # genes <- genes[!genes %in% c('','Background')]
        # if (length(genes)>0) segments[i]$genes <- paste(genes,collapse = ', ')
        
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

# Add gene related information to segments

segranges <- makeGRangesFromDataFrame(segments[,.(chromosome,start=start_pos,end=end_pos)])
generanges <- makeGRangesFromDataFrame(cancergenes)
exonranges <- makeGRangesFromDataFrame(cancerexons)

gene_overlap <- as.data.table(findOverlaps(segranges,generanges))
exon_overlap <- as.data.table(findOverlaps(segranges,exonranges))

for (i in 1:nrow(segments)) {
    gene_ix <- gene_overlap[queryHits==i]$subjectHits
    
    if (length(gene_ix)>0) {
        genetable <- cancergenes[gene_ix,.(ensembl_id,symbol,type,exonic=F,label='')]
        label <- paste(genetable$symbol,collapse=',')
        segments[i,genes:=label]
        
        exon_ix <- exon_overlap[queryHits==i]$subjectHits
        
        if (length(exon_ix)>0) {
            exontable <- cancerexons[exon_ix]
            genetable[ensembl_id %in% exontable$ensembl_id,exonic:=T]
        }
        genetable[exonic==T,label:=paste0(symbol,'|',type)]
        
        label <- paste(genetable[exonic==T]$label,collapse=',')
        segments[i,relevance:=label]
    }
}





# Output ------------------------------------------------------------



targets <- merge(alltargets,targets,by=colnames(alltargets),all=T)[order(bin)]

# Jumble targets and background
saveRDS(targets,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble.RDS'))
if (snp_allele_ratio) saveRDS(snp_table,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble_snps.RDS'))


## CNR/CNS ------
# cnr:  chromosome      start   end     gene    depth   log2    weight ()
cnr <- targets[!is.na(log2),.(chromosome=as.character(chromosome),start,end,gene,
                              depth=round(count/width*200,3),log2,weight=1,
                              gc,count,type)][gene=='',gene:='-']

fwrite(x = cnr,file = paste0(opt$output_dir,'/',clinbarcode,'.cnr'),sep = '\t')

# cns:  chromosome      start   end     gene    log2    depth   probes  weight
cns <- segments[,.(chromosome,start=start_pos,end=end_pos,
                   gene=genes,
                   log2=round(mean,2),
                   depth=round(2^mean,2),
                   probes=nbrOfLoci,
                   relevance)]
fwrite(x = cns,file = paste0(opt$output_dir,'/',clinbarcode,'.cns'),sep = '\t')


## DNAcopy segment file ----
# ID    chrom   loc.start       loc.end num.mark        seg.mean        C
seg <- segments[,.(ID=name,chrom=chromosome,loc.start=start_pos,loc.end=end_pos,num.mark=nbrOfLoci,seg.mean=mean,C=NA)]
seg[,chrom:=str_replace(chrom,'Y','24')][,chrom:=str_replace(chrom,'X','23')][,chrom:=as.numeric(chrom)]
fwrite(x = seg,file = paste0(opt$output_dir,'/',clinbarcode,'_dnacopy.seg'),sep = '\t')


## Count file ------
# (not overwrite, not if input was counts.RDS)
if (!file.exists(paste0(opt$output_dir,'/',clinbarcode,'.*counts.RDS')))
    if (!str_detect(opt$input_bam,'counts.RDS$'))
        saveRDS(counts,paste0(opt$output_dir,'/',clinbarcode,'.counts.RDS'))


#save.image(paste0(opt$output_dir,'/',clinbarcode,'.jumble_workspace.Rdata'))



# QC ------------

mapd <- function(data) {
    return(median(abs(diff(data)),na.rm = T))
} 

noise <- function(data) {
    m <- mapd(data)
    f <- 2^m-1
    return(round(100*f,1))
}

if (!this_is_wgs) {
    
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



if (this_is_wgs) stats <- paste0('Fragments per target: ',
                                 paste(round(quantile(alltargetcount,c(.025,.975))),collapse = '-'),
                                 ', Noise: ',
                                 noise(targets$log2),'%'
)




# Plot ------------------------------------------------------------



if (T) {
    
    
    
    p <- NULL
    targets[,smooth_log2:=runmed(log2,k=7),by=chromosome]
    ylims <- c(.4,max(2,max(2^targets$smooth_log2)))
    
    size <- 1; if (this_is_wgs) size <- 2
    
    if (this_is_wgs) {
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
    
    # if (snp_allele_ratio & !this_is_wgs) {
    #     
    #     layout <-  "ABBBB
    #                 CDDDD
    #                 EFFFF
    #                 GHHHH
    #                 IJJJJ
    #                 KLLLL"
    #     fig <-
    #         p$gc_rawdepth+p$order_rawdepth+
    #         p$gc_log2+p$order_log2+
    #         p$gc_log2_nosub+p$order_log2_nosub+
    #         p$gc_log2_simple+p$order_log2_simple+
    #         p$gc_log2_pca+p$order_log2_pca+
    #         p$depth_alleleratio+p$order_alleleratio+
    #         plot_layout(design = layout,guides = 'collect')
    # }
    
    if (snp_allele_ratio & !this_is_wgs) {

        layout <-  "ABBBB
                CDDDD
                EFFFF
                GGGGG
                HHHHH
                IJJJJ
                IJJJJ"
        fig <-
            p$gc_rawdepth+p$order_rawdepth+
            p$gc_log2+p$order_log2+
            p$depth_alleleratio+p$order_alleleratio+
            p$pos_log2+
            p$pos_alleleratio+
            p$nogrid+p$grid+
            plot_layout(design = layout,guides = 'collect')
    }
    if (!snp_allele_ratio & !this_is_wgs) {
        
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
    
    if (this_is_wgs) {
        
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


