#!/usr/bin/env Rscript

# Markus Mayrhofer 2022-2024.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

jumble_version <- 0.6

# Dependencies ------------------------------------------------------------

mylibrary <- function(package_names, cran_mirror = "https://cloud.r-project.org") {
    # # Set default CRAN mirror
    # options(repos = c(CRAN = cran_mirror))
    #
    # # Install BiocManager if not already installed
    # if (!requireNamespace("BiocManager", quietly = TRUE))
    #     install.packages("BiocManager")
    #
    for (package_name in package_names) {
        suppressPackageStartupMessages({
            # if (!require(package_name, character.only = TRUE)) {
            #     # Attempt to install from Bioconductor first
            #     tryCatch({
            #         BiocManager::install(package_name, ask = FALSE)
            #     }, error = function(e) {
            #         # If Bioconductor installation fails, try CRAN
            #         install.packages(package_name, type = "binary")
            #     })
            #     # Check if the package is loaded successfully after installation
            #     if (!require(package_name, character.only = TRUE)) {
            #         stop("Package ", package_name, " could not be loaded.")
            #     }
            # }
            # Load the package
            library(package_name, character.only = TRUE)
        })
    }
}



mylibrary(c(
    'optparse',
    'data.table',
    'stringr',
    'bamsignals',
    'Rsamtools',
    'GenomicRanges',
    'MASS',
    'fastICA',
    'VariantAnnotation',
    'PSCBS',
    'ggplot2',
    'patchwork',
    'BSgenome.Hsapiens.UCSC.hg19',
    'BSgenome',
    'Repitools'
)); theme_set(theme_bw())




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


#save.image('ws.Rdata')

# Reference file ------------------------------------------------------------

reference <- readRDS(opt$reference_file)


counts_template <- reference[c("target_bed_file","chromlength","ranges")]

this_is_wgs <- reference$target_bed_file=='wgs'


use_medium_fraglength <- FALSE
use_medium_fraglength <-
    !this_is_wgs &
    !is.null(reference$allcounts[[1]]$count_medium) &
    str_detect(opt$input_bam,'clip')



# Target template ----------

{
    targets <- reference$target_template
    targets[,chromosome:=str_remove(chromosome,'^chr')]

    ## Mark tiled
    targets[,left:=c(Inf,abs(diff(mid)))]
    targets[,right:=c(abs(diff(mid)),Inf)]
    targets[,farthest:=left][right>left,farthest:=right]
    targets[,is_tiled:=F][farthest<250,is_tiled:=T]
    targets[,left:=NULL][,right:=NULL][,farthest:=NULL]
    tiled_genes <- as.data.table(sort(table(targets[is_tiled==T]$gene),decreasing = T))
    if (sum(targets$is_tiled==T) > 500) targets[is_tiled==T,type:='tiled']

}

#save.image('ws.Rdata')

## Gene annotation -----------
{
    allgenes <- reference$allgenes
    allexons <- reference$allexons
    cancergenes_clinseq <- reference$cancergenes_clinseq
    
    # Extract the required cancer genes from all genes
    cancergenes <- allgenes[`Gene stable ID` %in% cancergenes_clinseq$ensembl_gene_id_version]
    
    # Check
    # missing <- cancergenes_clinseq[!hugo_symbol %in% cancergenes$`Gene name` &
    #                                    !hugo_symbol %in% cancergenes$`Gene Synonym` &
    #                                    !ensembl_gene_id_version %in% cancergenes$`Gene stable ID`]
    missing <- cancergenes_clinseq[!hugo_symbol %in% cancergenes$`Gene name`]
    


    # Remove duplicate rows
    cancergenes <- cancergenes[,.(id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                                  start=`Gene start (bp)`,end=`Gene end (bp)`,
                                  symbol=`Gene name`)]
    cancergenes <- unique(cancergenes)

    # Add ambigous, oncogene or tsg
    cancergenes[,type:='A']
    cancergenes[id %in% cancergenes_clinseq[ANNOT=='ONCO']$ensembl_gene_id_version,type:='O']
    cancergenes[symbol %in% cancergenes_clinseq[ANNOT=='ONCO']$hugo_symbol,type:='O']
    cancergenes[id %in% cancergenes_clinseq[ANNOT=='TSG']$ensembl_gene_id_version,type:='T']
    cancergenes[symbol %in% cancergenes_clinseq[ANNOT=='TSG']$hugo_symbol,type:='T']

    # Extract related exons
    cancerexons <- unique(allexons[`Gene stable ID` %in% cancergenes$id,
                            .(id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                              #exon=`Exon rank in transcript`,
                              start=`Exon region start (bp)`,end=`Exon region end (bp)`)])
    allexons <- unique(allexons[,
                            .(id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                              #exon=`Exon rank in transcript`,
                              start=`Exon region start (bp)`,end=`Exon region end (bp)`)])

    # All genes
    allgenes <- allgenes[,.(id=`Gene stable ID`,chromosome=`Chromosome/scaffold name`,
                            start=`Gene start (bp)`,end=`Gene end (bp)`,
                            symbol=`Gene name`)]
    allgenes <- unique(allgenes)


    # Make ranges objects
    cancergeneranges <- makeGRangesFromDataFrame(cancergenes)
    generanges <- makeGRangesFromDataFrame(allgenes)
    exonranges <- makeGRangesFromDataFrame(allexons)
    binranges <- makeGRangesFromDataFrame(targets)


    # Gene overlap table
    gene_overlap <- as.data.table(findOverlaps(binranges,generanges))
    # Add symbol to overlap table
    gene_overlap[,symbol:=allgenes[subjectHits]$symbol]

    for (i in unique(gene_overlap$queryHits)) {
        targets[i,gene:=paste0(gene_overlap[queryHits==i]$symbol,collapse = ',')]
    }

    # Exons overlap table
    exon_overlap <- as.data.table(findOverlaps(binranges,exonranges))
    targets[unique(exon_overlap$queryHits),type:='exonic']
    targets[is_target==F,type:='background']

}

## Backbone definition ------------------------------------------------------------

set.seed(25) # <------------------ To be reproducible.
max_in_b <- 25

targets[,is_backbone:=chromosome %in% c(1:22,'X')]

if (!this_is_wgs) {
    genes <- unique(targets[gene!='' & chromosome %in% 1:22]$gene)
    genes <- genes %in% c('BRCA2','PTEN') # To select just some genes for this restriction.
    for (g in genes) {
        genebins <- targets[gene==g]$bin
        n <- length(genebins)
        if (n>max_in_b) {
            targets[bin %in% genebins,is_backbone:=F]
            genebins=sample(genebins,max_in_b,replace = F)
            targets[bin %in% genebins,is_backbone:=T]
        }
    }
}

## Template complete ------

target_template <- targets
#table(targets$type,targets$is_target)


# Manage reference set ----------


## Iterate reference samples ------------------------------------------------------------

allcounts <- reference$allcounts
targetlist <- NULL
for (i in 1:length(allcounts)) {

    counts <- allcounts[[i]]
    targets <- copy(target_template)
    name <- paste0('ref_',i)
    targets$sample <- name


    # Add counts
    targets[,count:=counts$count]
    if (use_medium_fraglength) targets[,count:=counts$count_medium]
    targets[,count_short:=counts$count_short]

    targetlist[[i]] <- targets
}

targets <- rbindlist(targetlist)
rm(targetlist)





## Remove worst bins ------------------------------------------------------------

# Low coverage threshold
threshold <- median(targets[is_target==T]$count) * 0.01
keep_targets <- targets[,median(count),by=bin][V1 > threshold]
targets <- targets[bin %in% keep_targets$bin]

# High coverage threshold
threshold <- median(targets[is_target==T]$count) / 0.05
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


## Median correct ------------------------------------------------------------

# Basic logR, targets (with correction for SNPs if available)
targets[,rawLR:=log2(min1(count))]
targets[,rawLR_short:=log2(min1(count_short))]

# Median correct, separated on sample and target/background status.
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by=c('sample','is_target')]
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by=c('sample','is_target')]


## X-Y chromosome correct ------------------------------------------------------------

if (T) {
    # Doctor the X values in samples where their median implies male
    targets[,xmedian:=median(rawLR[chromosome=='X' & is_tiled==F]),by=.(sample,is_target)]
    targets[,xmedian_short:=median(rawLR_short[chromosome=='X' & is_tiled==F]),by=.(sample,is_target)]
    targets[,male:=2^xmedian < .75] # assign gender

    targets[,nonPA:=chromosome %in% c('X') & end>2.70e6 & start<154.93e6] # hard coded PA

    # adjust to the median
    if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
        targets[chromosome=='X' & male==T & nonPA,rawLR:=rawLR+1]
        targets[chromosome=='X' & male==T & nonPA,rawLR_short:=rawLR_short+1]
    }

    # Doctor Y values where their median implies male
    targets[,ymedian:=median(rawLR[chromosome=='Y']),by=sample] # median by sample
    targets[,male:=2^ymedian > .25] # assign gender based on Y

    if (length(unique(targets[male==T]$sample))>0) { # if at least 1 male
        targets[male==T & chromosome=='Y' & end<28.79e6,rawLR:=rawLR+1] # hard coded PA
        targets[male==T & chromosome=='Y' & end<28.79e6,rawLR_short:=rawLR_short+1]
    }

    # Y values set to NA where X median implied female
    targets[chromosome=='Y' & male==F,rawLR:=NA]
    targets[chromosome=='Y' & male==F,rawLR_short:=NA]

}

## Bin median correct ------------------------------------------------------------
targets[,refmedian:=median(rawLR,na.rm=T),by=bin][,rawLR:=rawLR-refmedian]
targets[,refmedian_short:=median(rawLR_short,na.rm=T),by=bin][,rawLR_short:=rawLR_short-refmedian_short]

bins_for_mediancorrect <- unique(targets[,.(bin,refmedian,refmedian_short)])

## Impute missing ----------------------------------------------------------

# If any missing, replace with random value near 0 (which is already the median)
targets[is.na(rawLR),rawLR:=rnorm(n = .N,mean = 0,sd = .01)]
targets[is.na(rawLR_short),rawLR_short:=rnorm(n = .N,mean = 0,sd = .01)]



reference_targets <- targets





# Manage query sample ------------------------------------------------------------

# function for bam process
countsFromBam <- function(counts,bampath) {
    counts$date_count <- date()
    counts$input_bam_file <- bampath
    ranges <- counts$ranges
    flag <- 1024
    if (!is.null(reference$flag)) flag <- reference$flag
    counts[['count']] <- bamCount(bampath, ranges,
                                  paired.end="midpoint",
                                  mapq=20,
                                  filteredFlag=flag,
                                  verbose=F)
    counts[['count_short']] <- bamCount(bampath, ranges,
                                        paired.end="midpoint",
                                        mapq=20,
                                        filteredFlag=flag,
                                        tlenFilter=c(0,150),
                                        verbose=F)
    counts[['count_medium']] <- bamCount(bampath, ranges,
                                         paired.end="midpoint",
                                         mapq=20,
                                         filteredFlag=flag,
                                         tlenFilter=c(0,300),
                                         verbose=F)
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

clinbarcode <- name #str_remove(name, "[_-]nodups.bam")
#clinbarcode <- str_remove(name, "[_-]clipoverlap.bam")

targets$sample <- clinbarcode

# Add counts
targets[,count:=counts$count]
if (use_medium_fraglength) targets[,count:=counts$count_medium]

# same, short
targets[,count_all:=counts$count]
if (!is.null(counts$count_medium)) targets[,count_medium:=counts$count_medium]
targets[,count_short:=counts$count_short]

# for stats
alltargetcount <- targets[is_target==T]$count


# For SNPs if present
targets[,snps:=0]


target_ranges <- makeGRangesFromDataFrame(targets[is_target==T])
ranges <- makeGRangesFromDataFrame(targets)




## SNP allele ratio ------------------------------------------------------------
#save.image('ws.Rdata')

peakx <- function(data,weights=NULL,na.rm = T) {
    if(length(data)==1) return(data)
    maxx <- NA
    try({
        if (is.null(weights)) d <- density(data[!is.na(data)])
        if (!is.null(weights)) d <- density(data[!is.na(data)],weights=weights[!is.na(data)])
        maxy <- max(d$y)
        maxx <- d$x[d$y==maxy][1]
    }, silent=T)
    return(maxx)
}

snp_allele_ratio <- FALSE
input <- opt$snp_vcf
if (!is.null(input)) {
    snp_allele_ratio <- TRUE

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
    rr <- rowRanges(vcf)
    snp_table <- as.data.table(rr)[,1:3][,chromosome:=str_remove(seqnames,'^chr')][,c(4,2,3)]
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

}




## LogR correction ------------------------------------------------------------
#save.image('ws.Rdata')

targets <- targets[bin %in% bins_for_mediancorrect$bin]


min1 <- function(data) {
    data[data<1] <- 1
    data[is.na(data)] <- 1
    return(data)
}

# Basic logR, targets (with correction for SNPs)
targets[,rawLR:=log2(min1(count))]
targets[,rawLR_short:=log2(min1(count_short))]

# median correct to backbone, separately for targets/background
targets[,rawLR:=rawLR-median(rawLR[is_backbone]),by='is_target']
targets[,rawLR_short:=rawLR_short-median(rawLR_short[is_backbone]),by='is_target']

# correct by reference median
targets[,rawLR:=rawLR-bins_for_mediancorrect$refmedian]
targets[,rawLR_short:=rawLR_short-bins_for_mediancorrect$refmedian_short]



## Retain all bins in output ------------------------------------------------------------
alltargets <- copy(targets)



# Correct using reference data ------------------------------------------------------------

## Remove query sample from reference --------------------
if (length(unique(reference_targets$sample))>1) for (s in unique(reference_targets$sample)){
    if (all(reference_targets[sample==s]$count_short==targets$count_short)) {
        reference_targets <- reference_targets[sample!=s]
        print('Removed this sample from reference')
    }
}



## References on matrix form ------------------------------------------------------------
mat <- dcast(data = reference_targets[,.(bin,sample,rawLR)],formula = bin ~ sample, value.var = 'rawLR')
mat_short <- dcast(data = reference_targets[,.(bin,sample,rawLR_short)],formula = bin ~ sample, value.var = 'rawLR_short')



## PCA 1: outliers ------------------------------------------------------------

targetbins <- target_template[is_target==T]$bin # the ontarget

set.seed(25)
tpca <- as.data.table(prcomp(mat[bin %in% targetbins,-1],center = F,scale. = F)$x)
tpca$keep <- T

for (pc in colnames(tpca)[-ncol(tpca)][1:min(100,ncol(tpca))]) {
    fact <- ifelse(pc %in% paste0('PC',1:3),3,4)
    sd <- sd(tpca[[pc]])
    tpca[tpca[[pc]] < -sd*fact, keep:=F]
    tpca[tpca[[pc]] > sd*fact, keep:=F]
}

remove <- mat[bin %in% targetbins]$bin[tpca$keep==F]


if (!this_is_wgs) if (any(targets$is_target==F)) {
    backgroundbins <- target_template[is_target==F]$bin # the offtarget
    bgpca <- as.data.table(prcomp(mat[bin %in% backgroundbins,-1],center = F,scale. = F)$x)
    bgpca$keep <- T
    
    for (pc in colnames(bgpca)[-ncol(bgpca)][1:min(100,ncol(bgpca))]) {
        fact <- ifelse(pc %in% paste0('PC',1:3),3,4)
        sd <- sd(bgpca[[pc]])
        bgpca[bgpca[[pc]] < -sd*fact, keep:=F]
        bgpca[bgpca[[pc]] > sd*fact, keep:=F]
    }
    
    remove <- c(remove,mat[bin %in% backgroundbins]$bin[bgpca$keep==F])
}

mat <- mat[!bin %in% remove]
mat_short <- mat_short[!bin %in% remove]
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


#save.image('ws.Rdata')


## Correct using reference PCA  ---------------------------

targets <- targets[bin %in% bins_outliers_removed]


jcorrect <- function(temp,train_ix=NULL, mult=T) {
    
    if (is.logical(train_ix[1])) train_ix <- which(train_ix)
    if (is.null(train_ix)) train_ix <- which(rep(TRUE,nrow(temp)))
    
    if (length(train_ix)>20e3) train_ix <- sample(train_ix,size = 20e3,replace = F)
    
    pcs <- sum(str_detect(colnames(temp),'^PC'))
    
    if (pcs > 50) pcs <- 50
    
    
    bins <- nrow(temp)
    
    
    # Dynamically construct formula
    formula_str <- paste("lr ~ ", paste0("PC", 1:min(50,pcs), collapse = "+"))
    formula <- as.formula(formula_str)
    rlm_temp <- rlm(formula, data=temp,
                    subset = train_ix)
    temp[,lr:=lr-predict(rlm_temp,temp)]
    
    # Sequential
    if (T) for (i in 1:(pcs)) {
        temp$thispc=temp[[paste0('PC',i)]]
        rlm_temp <- rlm(lr ~ poly(thispc,degree = 1), data=temp,
                        subset = train_ix)
        temp[,lr:=lr-predict(rlm_temp,temp)]
    }
    
    # GC correct
    if (T) if (all(!is.na(temp$gc))) {
        loess_temp=loess(lr ~ gc, data = temp,
                         subset = train_ix,
                         family="symmetric", control = loess.control(surface = "interpolate"))
        pred <- predict(loess_temp,temp)
        temp$lr <- temp$lr-pred
    }
    
    return(temp$lr)
}



targets[,log2:=rawLR]
targets[,log2_short:=rawLR_short]

# standard
ix <- targets$is_target# %in% c(T,F)
temp <- cbind(data.table(
    lr=targets[ix]$rawLR),
    gc=targets[ix]$gc,
    loc=targets[ix,paste(chromosome,round(mid/20e6),sep=':')],
    gene=targets[ix]$gene,
    type=targets[ix]$type,
    tpca)
targets[ix,log2:=jcorrect(temp,targets[ix]$is_backbone)]


# short
temp <- cbind(data.table(
    lr=targets[ix]$rawLR_short),
    gc=targets[ix]$gc,
    loc=targets[ix,paste(chromosome,round(mid/20e6),sep=':')],
    gene=targets[ix]$gene,
    type=targets[ix]$type,
    tpca_short)
targets[ix,log2_short:=jcorrect(temp,targets[ix]$is_backbone)]


if (!this_is_wgs) {

    if (any(targets$is_target==F)) {
        ix <- targets$is_target
        # standard bg
        temp <- cbind(data.table(
            lr=targets[!ix]$rawLR),
            gc=targets[!ix]$gc,
            loc=targets[!ix,paste(chromosome,round(mid/100e6),sep=':')],
            gene=targets[!ix]$gene,
            bgpca)
        targets[!ix,log2:=jcorrect(temp,targets[!ix]$is_backbone)]

        # bg short
        temp <- cbind(data.table(
            lr=targets[!ix]$rawLR_short),
            gc=targets[!ix]$gc,
            loc=targets[!ix,paste(chromosome,round(mid/100e6),sep=':')],
            gene=targets[!ix]$gene,
            bgpca_short)
        targets[!ix,log2_short:=jcorrect(temp,targets[!ix]$is_backbone)]

    }
}



## Set min/max ------------------------------------------------------------

targets[log2 < -5, log2:=-5]
targets[log2 > 7, log2:=7]
targets[log2_short < -4, log2_short:=-4]
targets[log2_short > 7, log2_short:=7]
targets[chromosome=='Y',log2:=NA]
targets[chromosome=='Y',log2_short:=NA]


#save.image('ws.Rdata')
## Adjust X to background ------------------------------------------------------------

if (!this_is_wgs) try(
    {
        # X-chromosome correction factor to targeted bins, based on background bins
        temp <- targets[chromosome=='X' & !is.na(log2)]
        temp[,pos_1M:=1*round(start / 1e6)]
        temp[,bg:=as.numeric(NA)][,bg:=median(log2[is_target==F],na.rm = T),by=pos_1M]
        temp[,tg:=as.numeric(NA)][,tg:=median(log2[is_target==T & is_tiled==F],na.rm = T),by=pos_1M]

        temp <- unique(temp[,.(pos_1M,bg,tg)])

        temp[,bg_median:=runmed(bg,11,na.action = 'na.omit')]
        temp[,tg_median:=runmed(tg,11,na.action = 'na.omit')]
        temp[,dif:=bg_median-tg_median]
        x_correct <- median(temp$dif,na.rm = T)
        targets[chromosome=='X' & is_target==T, log2:=log2+x_correct]
    }, silent=T
)

## Segmentation ------------------------------------------------------------

print('Segmentation')

#save.image('ws.Rdata')

getsegs <- function(targets, logratio) {

    alpha <- .02
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


        # adjust segment start and end to bin number, rather than bin order
        ix <- ceiling(segments[i]$start)
        segments[i,start:=targets[ix]$bin]
        ix <- floor(segments[i]$end)
        segments[i,end:=targets[ix]$bin]
    }


    return(segments)
}



# Do the segmentation
if (T) {

    targets[,chromosome:=str_replace(chromosome,'Y','24')][,chromosome:=str_replace(chromosome,'X','23')][,chromosome:=as.numeric(chromosome)]
    segments <- getsegs(targets, targets$log2)
    targets[,chromosome:=as.character(chromosome)][chromosome=='23',chromosome:='X'][chromosome=='24',chromosome:='Y']
    segments[,chromosome:=as.character(chromosome)][chromosome=='23',chromosome:='X'][chromosome=='24',chromosome:='Y']

    # Add cancer gene related information to segments

    segranges <- makeGRangesFromDataFrame(segments[,.(chromosome,start=start_pos,end=end_pos)])
    generanges <- makeGRangesFromDataFrame(cancergenes)
    exonranges <- makeGRangesFromDataFrame(cancerexons)

    gene_overlap <- as.data.table(findOverlaps(segranges,generanges))
    exon_overlap <- as.data.table(findOverlaps(segranges,exonranges))

    for (i in 1:nrow(segments)) {
        gene_ix <- gene_overlap[queryHits==i]$subjectHits

        if (length(gene_ix)>0) {
            genetable <- cancergenes[gene_ix,.(id,symbol,type,exonic=F,label='')]
            label <- paste(genetable$symbol,collapse=',')
            segments[i,genes:=label]

            exon_ix <- exon_overlap[queryHits==i]$subjectHits

            if (length(exon_ix)>0) {
                exontable <- cancerexons[exon_ix]
                genetable[id %in% exontable$id,exonic:=T]
            }

            genetable[exonic==T,label:=symbol]
            genetable[exonic==T & type!='',label:=paste0(symbol,'|',type)]

            label <- paste(genetable[exonic==T]$label,collapse=',')
            segments[i,relevance:=label]
        }
    }
}



# Output ------------------------------------------------------------

# Reintroduce the removed bins
targets <- merge(alltargets,targets,by=colnames(alltargets),all=T)[order(bin)]

# Jumble targets and background
saveRDS(targets,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble.RDS'))
if (snp_allele_ratio) saveRDS(snp_table,file = paste0(opt$output_dir,'/',clinbarcode,'.jumble_snps.RDS'))

if (this_is_wgs) targets[,type:='bin']


if (T) {

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



    noise <- function(targets) {
        targets[,loc:=paste(chromosome,round(start/10e6))]
        targets[,n:=.N,by=loc]
        std <- targets[n>=10,sd(log2),by=loc]
        round(2^median(std$V1,na.rm=T)*100-100,1)
    }

    if (!this_is_wgs) {
        
        stats <- paste0('Fragments per target: ',
                        paste(round(quantile(alltargetcount,c(.1,.9))),collapse = '-'),
                        ', Noise: ',
                        noise(targets),'%'
        )
    }
    
    
    
    if (this_is_wgs) stats <- paste0('Fragments per target: ',
                                     paste(round(quantile(alltargetcount,c(.1,.9))),collapse = '-'),
                                     ', Noise: ',
                                     noise(targets),'%'
    )
    



    # Plot ------------------------------------------------------------

    if (T) suppressWarnings( {
        
        
        
        p <- NULL
        targets[,smooth_log2:=runmed(log2,k=7),by=chromosome]
        ylims <- c(
            min(.4,min(2^targets[chromosome!='Y']$smooth_log2,na.rm = T)),
            max(3,max(2^targets$smooth_log2,na.rm = T))
        )
        ybreaks <- c(.5,.75,1,1.5,2,3,4,6,8)
        ybreaks <- ybreaks[ybreaks >= ylims[1] & ybreaks<=ylims[2]]
        yminorbreaks <- c(1.25,1.75)
        
        colorvalues <- c('background'='#909090','dense'='#F8766D',
                         'exonic'='#61ACFF','sparse'='#252525', 'bin'='#252525')
        
        size <- 1#; if (this_is_wgs) size <- 2
        
        pointcolor <- '#50505060'
        
        
        targets[,label:=type][label=='tiled',label:='dense'][label=='target',label:='sparse']
        
        
        alpha <- .5
        if (this_is_wgs) alpha <- .2
        
        
        ## Grid ------------------------------------------------------------
        
        if (snp_allele_ratio) {
            
            # defaults to using raw allele ratio:
            snp_table[,allele_ratio_use:=allele_ratio]
            # but if there is a corrected allele ratio, use it:
            if (!is.null(snp_table$allele_ratio_corrected2)) snp_table[,allele_ratio_use:=allele_ratio_corrected2]
            
            snp_table <- snp_table[type!='other'][allele_ratio_use < .99][allele_ratio_use > .01][RD>2][AD>2]
            snp_table <- snp_table[DP > median(DP)/5][DP < median(DP)*10]
            
            targets[,allele_ratio:=as.double(NA)][match(snp_table$bin,bin),allele_ratio:=snp_table$allele_ratio_use]
            targets[,maf:=abs(allele_ratio-.5)+.5]
            targets[!is.na(maf),maf:=runmed(maf,9)] # smoothing!
            # snp (grid) smooth-to-allele-ratio plot
            p$grid <- ggplot(targets) + xlim(c(.2,1.8)) + scale_y_continuous(limits=c(.5,1)) +
                xlab('Corrected depth (smooth)') + ylab('Major Allele Ratio') +
                geom_point(data=targets[,.(smooth_log2,maf)],aes(x=2^smooth_log2,y=maf),col='lightgrey',alpha=.2) +
                #geom_point(data=targets[,.(gc,log2)],aes(x=gc,y=2^smooth_log2),col='lightgrey',alpha=.2) +
                
                geom_text(data=unique(targets[,.(chromosome)]),mapping = aes(x=.25,y=.95,label=chromosome)) +
                
                #geom_point(aes(x=gc,y=2^smooth_log2),fill='#60606090',col='#20202090',shape=21) +
                geom_point(data=targets[label!='background'],aes(x=2^smooth_log2,y=maf,fill=label,col=label),shape=21,size=size,alpha=alpha) +
                facet_wrap(facets = vars(factor(chromosome,levels=unique(chromosome),ordered=T)),ncol = 8) +
                scale_fill_manual(values = colorvalues) +
                scale_color_manual(values = colorvalues) +
                #theme(panel.spacing = unit(0, "lines"),strip.text.x = element_text(size = 8)) +
                theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.text.x = element_blank())
            
            # snp (all) smooth-to-allele-ratio plot
            #temp <- targets[!is.na(label),median(log2),by=label]
            p$nogrid <- ggplot(targets) + xlim(c(0.2,1.8)) + ylim(c(.5,1)) + xlab('Corrected depth (smooth)') + ylab('Major allele ratio (smooth)') +
                #geom_point(data=targets[,.(smooth_log2,maf)],aes(x=2^smooth_log2,y=maf),fill='#60606050',col='#20202050',shape=21,alpha=alpha) +
                scale_fill_manual(values = colorvalues) +
                scale_color_manual(values = colorvalues) +
                geom_point(data=targets[label!=''],aes(x=2^smooth_log2,y=maf,fill=label,col=label),shape=21,size=size,alpha=alpha) #+
            #geom_point(data=temp,mapping=aes(x=2^V1,y=1,fill=label),size=2,shape=25,show.legend=F)
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
            #geom_point(data=targets[is.na(label)],mapping = aes(x=gpos,y=2^log2),fill='#60606070',col='#20202010',shape=21,size=1,alpha=alpha) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=2^log2,fill=label,col=label),shape=21,size=size,alpha=alpha) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=ylims, breaks=ybreaks, minor_breaks=yminorbreaks) +
            geom_segment(data=segments,col='green',linewidth=1,
                         mapping = aes(x=gstart,xend=gstop,y=2^mean,yend=2^mean)) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line())
        
        # depth by pos
        limits <- c(1e1,1e5)
        limits_labels <- c('10','30','100','300','1k','3k','10k','30k','100k')
        limits_breaks <- c(1e1,3e1,1e2,3e2,1e3,3e3,1e4,3e4,1e5)
        
        p$pos_rawdepth <- ggplot(targets) + xlab('Genomic position') + ylab('Fragments') +
            #geom_point(data=targets[is_target==T],mapping = aes(x=gpos,y=count),fill='#60606050',col='#20202010',size=1,shape=21,alpha=alpha) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=count,fill=label,col=label),shape=21,size=size,alpha=alpha) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=limits,breaks=limits_breaks,minor_breaks=yminorbreaks,labels=limits_labels) +
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
                #geom_point(data=targets[is.na(label)],mapping = aes(x=gpos,y=allele_ratio),fill='#60606080',col='#20202010',shape=21,size=1,alpha=alpha) +
                geom_point(data=targets[!is.na(label)],mapping = aes(x=gpos,y=allele_ratio,fill=label,col=label),shape=21,size=1,alpha=alpha/2) +
                scale_fill_manual(values = colorvalues) + ylim(0:1) +
                scale_color_manual(values = colorvalues) +
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
            geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2,fill=label,col=label),shape=21,size=size,alpha=alpha) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=ylims, breaks=ybreaks, minor_breaks=yminorbreaks) +
            geom_segment(data=segments,col='green',linewidth=1,
                         mapping = aes(x=start,xend=end,y=2^mean,yend=2^mean)) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line())
        # logR by gc
        p$gc_log2 <- ggplot(targets) + xlab('Target GC content') + ylab('Corrected depth') + xlim(c(.16,.84)) +
            geom_point(data=targets,mapping = aes(x=gc,y=2^log2,fill=label,col=label),
                       shape=21,size=1,alpha=alpha/2,show.legend = F) + # fill='#60606040'
            facet_wrap(facets = vars(label),nrow = 1) +
            theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.text.x = element_blank()) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=ylims, breaks=ybreaks, minor_breaks=yminorbreaks)
        
        # short by order
        p$order_log2_short <- ggplot(targets) + xlab('Order of genomic position') + ylab('Corrected depth') +
            #geom_point(data=targets[is.na(label) & is_target==T],mapping = aes(x=bin,y=2^log2),fill='#60606070',col='#20202010',size=1,alpha=alpha) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=2^log2_short,fill=label,col=label),shape=21,size=size,alpha=alpha) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=ylims, breaks=ybreaks, minor_breaks=yminorbreaks) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line())
        # short by gc
        p$gc_log2_short <- ggplot(targets) + xlab('Target GC content') + ylab('Corrected depth') + xlim(c(.16,.84)) +
            geom_point(data=targets,mapping = aes(x=gc,y=2^log2_short,fill=label,col=label),
                       shape=21,size=1,alpha=alpha/2,show.legend = F) + # fill='#60606040'
            facet_wrap(facets = vars(label),nrow = 1) +
            theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.text.x = element_blank()) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=ylims, breaks=ybreaks, minor_breaks=yminorbreaks)
        
        
        ## Allele ratio ------------------------------------------------------------
        
        m <- quantile(targets[is_target==T]$count,c(.01, .50, .99))
        if (m[1] > m[2]/3) m[1] <- m[2]/3
        if (m[3] < m[2]*3) m[3] <- m[2]*3
        m[1] <- max(1,m[1])
        if (snp_allele_ratio) {
            # allele ratio by order
            p$order_alleleratio <- ggplot(targets) + xlab('Order of genomic position') + ylab('Allele ratio') +
                #geom_point(data=targets[is.na(label)],mapping = aes(x=bin,y=allele_ratio),fill='#60606050',col='#20202010',shape=21,size=1,alpha=alpha) +
                geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=allele_ratio,fill=label,col=label),shape=21,size=size,alpha=alpha) +
                scale_fill_manual(values = colorvalues) +
                scale_color_manual(values = colorvalues) +
                scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                                   expand = c(.01,.01),labels = chroms$chromosome) +
                scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1)) +
                theme(panel.grid.major.x = element_blank(),
                      panel.grid.minor.y = element_line(),
                      panel.grid.minor.x = element_line(color = 'black'),
                      axis.line = element_line(),
                      axis.ticks = element_line())
            # allele ratio by depth
            p$depth_alleleratio <- ggplot(targets[type!='background']) + xlab('Fragments') + ylab('Allele ratio') +
                #geom_point(data=targets[is.na(label)],mapping = aes(x=count,y=allele_ratio),fill='#60606070',col='#20202070',shape=21,size=1) +
                geom_point(mapping = aes(x=count,y=allele_ratio,fill=label,col=label),
                           shape=21,size=1,alpha=alpha/2,show.legend = F) +
                #facet_wrap(facets = vars(label),nrow = 1) +
                theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.text.x = element_blank()) +
                scale_fill_manual(values = colorvalues) +
                scale_color_manual(values = colorvalues) +
                scale_x_log10(limits=c(m[1],m[3])) +
                scale_y_continuous(limits=c(-.05,1.05),breaks=c(0,.25,.5,.75,1))
        }
        
        # Depth ------------------------------------------------------------
        
        # depth by order
        p$order_rawdepth <- ggplot(targets) + xlab('Order of genomic position') + ylab('Fragments') +
            #geom_point(data=targets[is_target==T],mapping = aes(x=bin,y=count),fill='#60606050',col='#20202010',size=1,shape=21,alpha=alpha) +
            geom_point(data=targets[!is.na(label)],mapping = aes(x=bin,y=count,fill=label,col=label),shape=21,size=size,alpha=alpha) +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=limits,breaks=limits_breaks,minor_breaks=yminorbreaks,labels=limits_labels) +
            scale_x_continuous(breaks = chroms$mid,minor_breaks = chroms$start[-1],
                               expand = c(.01,.01),labels = chroms$chromosome) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_line(),
                  panel.grid.minor.x = element_line(color = 'black'),
                  axis.line = element_line(),
                  axis.ticks = element_line())
        # depth by GC
        p$gc_rawdepth <- ggplot(targets) + xlab('Target GC content') + ylab('Fragments') + xlim(c(.16,.84)) +
            geom_point(data=targets,mapping = aes(x=gc,y=count,fill=label,col=label),
                       shape=21,size=1,alpha=alpha/2,show.legend = F) + #
            facet_wrap(facets = vars(label),nrow = 1) +
            theme(panel.spacing = unit(0, "lines"),strip.background = element_blank(),strip.text.x = element_blank()) +
            geom_text(data=unique(targets[,.(label)])[label=='target',label:='sparse'],
                      mapping=aes(label = label, x = 0.16, y = Inf), hjust = 0, vjust = 1.5) +
            # geom_smooth(data=targets[!is.na(label)],
            #             mapping = aes(x=gc,y=count,col=label),size=.5,se=F,show.legend = F,method = 'loess') +
            scale_fill_manual(values = colorvalues) +
            scale_color_manual(values = colorvalues) +
            scale_y_log10(limits=limits,breaks=limits_breaks,minor_breaks=yminorbreaks,labels=limits_labels)
        
        
        for (i in 1:length(p)) p[[i]] <- p[[i]] + guides(fill=guide_legend(override.aes=list(shape=21,size=3))) + labs(fill = NULL)
        
        
        pa <- plot_annotation(
            title = paste(clinbarcode,'         ',stats),
            caption = paste('Jumble',jumble_version,'on',format(Sys.time(), "%a %b %e %Y, %H:%M"))
        )
        
        # Layout ------------------------------------------------------------
        
        
        if (snp_allele_ratio & !this_is_wgs) {
            
            layout <-  "ABBB
                CDDD
                EFFF
                GGGG
                HHHH
                IJJJ
                IJJJ"
            fig <-
                p$gc_rawdepth+p$order_rawdepth+
                p$gc_log2+p$order_log2+
                p$depth_alleleratio+p$order_alleleratio+
                p$pos_log2+
                p$pos_alleleratio+
                p$nogrid+p$grid+
                plot_layout(design = layout,guides = 'collect') &
                theme(legend.position='none')
        }
        if (!snp_allele_ratio & !this_is_wgs) {
            
            layout <-  "ABBB
                CDDD
                EEEE
                "
            fig <-
                p$gc_rawdepth+p$order_rawdepth+
                p$gc_log2+p$order_log2+
                p$pos_log2+
                plot_layout(design = layout,guides = 'collect') &
                theme(legend.position='none')
        }
        
        
        if (this_is_wgs) {
            
            layout <-  "ABBB
                CDDD
                "
            fig <-
                p$gc_rawdepth+p$pos_rawdepth+
                p$gc_log2+p$pos_log2+
                plot_layout(design = layout,guides = 'collect') &
                theme(legend.position='none')
        }
        
        png(file = paste0(opt$output_dir,'/',clinbarcode,'.png'),width = 1600,height=1300,res=100)
        print(fig+pa)
        
        # Close image ------------------------------------------------------------
        
        
        dev.off()

        
    })
    
}

