#!/usr/bin/env Rscript

jumble_version <- '0.3'

# Markus Mayrhofer 2022-2023

# Dependencies ------------------------------------------------------------
{
    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(GenomicRanges))
    #suppressPackageStartupMessages(library(biomaRt))
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    suppressPackageStartupMessages(library(BSgenome))
    suppressPackageStartupMessages(library(Repitools))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(VariantAnnotation))
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(MASS))
    

}
# Options ------------------------------------------------------------

# Define options
option_list <- list(
    make_option(c("-i", "--input_folder"), action = "store", type = "character",default = '.', 
                help = "folder of count RDS files"),
    make_option(c("-o", "--output_folder"), action = "store", type = "character",default = '.', 
                help = "folder to put reference file"),
    make_option(c("-c", "--cores"), action = "store", type = "character",default = '4', 
                help = "cores to use")
)
opt <- parse_args(OptionParser(option_list = option_list))


if (!is.null(opt$cores)) registerDoParallel(cores=as.integer(opt$cores))

# opt <- list(input_folder='.',output_folder='.',cores=4)


# Read count files --------------------------------------------------------

files=dir(path = opt$input_folder,pattern = 'counts.RDS$',full.names = T)

# read files
{
    allcounts <- NULL
    ntargets <- 0
    bed_files <- ''
    for (i in 1:length(files)) {
        counts <- readRDS(files[i])
        if (is.null(counts$input_bam_file))
            counts$input_bam_file <- files[i]
        ntargets[i] <- length(counts$count)
        if (!is.null(counts$target_bed_file)) 
            bed_files[i] <- str_remove(string = counts$target_bed_file,pattern = '^.*/') else bed_files[i] <- 'wgs'
        allcounts[[i]] <- counts
    }
    
    if (length(table(ntargets))>1) stop('Number of bins differs between samples.')
    if (length(table(bed_files))>1) stop('BED file differs between samples.')
}

wgs <- F
if (is.null(allcounts[[1]]$target_bed_file)) wgs <- T

# Make targets template ------------------------------------------------------------

counts <- allcounts[[1]]

targets <- as.data.table(counts$ranges)

targets <- targets[,.(sample='',
                      bin=1:.N,
                      is_target=T,
                      type='target',
                      is_tiled=F,
                      chromosome=as.character(seqnames),
                      start,end,
                      mid=round((end+start)/2),
                      width,
                      gene='')]
targets[width!=min(width),is_target:=F]
targets[is_target==F,type:='background']


# GC content ------------------------------------------------------------

ucsc_ranges <- counts$ranges
seqlevelsStyle(ucsc_ranges) <- "UCSC"
targets[,gc:=as.double(NA)]
targets[is_target %in% c(T,F)]$gc <- gcContentCalc(ucsc_ranges , organism=Hsapiens)

# Mappability ------------------------------------------------------------

targets[,map:=as.double(NA)]
targets[is_target %in% c(T,F)]$map <- mappabilityCalc(ucsc_ranges , organism=Hsapiens)



# Gene tables ------------------------------------------------------------
cancergenes_clinseq <- fread('~/Analysis/genes_and_exons/ANNOT_FOR_CURATOR.txt')

# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# 
# allgenes <- as.data.table(getBM(attributes=c('ensembl_gene_id','ensembl_exon_id','hgnc_symbol','chromosome_name','start_position','end_position'),
#                mart = ensembl))
#cancergenes <- allgenes[ensembl_gene_id %in% cancergenes_clinseq$ensembl_gene_id_version]

allgenes <- fread('~/Analysis/genes_and_exons/mart_export.txt')[`Chromosome/scaffold name` %in% c(1:22,'X','Y')]
allexons <- 
    fread('~/Analysis/genes_and_exons/mart_export_exons.txt')[`Chromosome/scaffold name` %in% c(1:22,'X','Y')]




# Read VCF files -------------------

vcf_files=dir(path = opt$input_folder,pattern = '.vcf',full.names = T)

if (length(vcf_files) > 0) {

    if (length(vcf_files) != length(files))
        warning('Not one VCF file per sample?')

    snp_tables <- foreach(i=1:length(vcf_files)) %dopar% try({
        vcf <- readVcf(vcf_files[i])

        # remove SNPs that do not have 2 alleles
        alleles <- as.data.table(table(as.data.table(alt(vcf))$group))$N
        vcf <- vcf[alleles==1 & sapply(geno(vcf)$AD,length)==2]

        this_list <- NULL
        for (j in 1:ncol(vcf)) {
            name <- colnames(vcf)[j]

            g <- geno(vcf[,j])

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
            snp_table$DP <- snp_table$AD+snp_table$RD
            snp_table[,logDP:=log2(DP)]

            # Compute raw allele ratio
            raw_allele_ratio <- unname(round(snp_table$AD/snp_table$DP,4))
            raw_allele_ratio[is.nan(raw_allele_ratio)] <- 0
            snp_table$allele_ratio <- raw_allele_ratio
            snp_table[allele_ratio==0,type:='none']

            # Add to list of tables
            #fwrite(snp_table,paste0(name,'.snptable.csv'))
            this_list[[j]] <- snp_table
        }
        this_set <- rbindlist(this_list)
        return(this_set)
    }, silent=T)

    snp_table <- unique(rbindlist(snp_tables)) # in case one sample appears more than once

    
    snp_table[,snp:=paste(id,ref_allele,alt_allele)]

    # map snps to bins
    unique_table <- unique(snp_table[,.(id,chromosome,start,end)])
    overlap <- findOverlaps(makeGRangesFromDataFrame(unique_table),counts$ranges)
    unique_table[queryHits(overlap),bin:=subjectHits(overlap)]
    snp_table <- merge(snp_table,unique_table,by=c('id','chromosome','start','end'))

    # Drop SNPs not mapping to bins
    snp_table <- snp_table[bin %in% target_template[is_target==T]$bin]

    setkey(snp_table,'sample')
    rm(snp_tables)



    # SNP GC


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

}
    


# Reference object ------------------------------------------------------------
if (wgs) allcounts[[1]]$target_bed_file <- 'wgs'
reference <- allcounts[[1]][c("target_bed_file","chromlength","ranges")]
reference$date <- date()
reference$samples <- unique(targets$sample)
reference$target_template <- targets

reference$allcounts <- allcounts

reference$allgenes <- allgenes
reference$allexons <- allexons
reference$cancergenes_clinseq <- cancergenes_clinseq

if (exists('snp_table')) {
    reference$snp_table <- snp_table
}


# Save ------------------------------------------------------------

name <- reference$target_bed_file
if (is.null(name)) name <- 'jumble.WGS'
saveRDS(reference,paste0(opt$output_folder,'/',
                         str_remove(name,'.*/'),
                         '.reference.RDS'))

print(paste(name,'done.'))
