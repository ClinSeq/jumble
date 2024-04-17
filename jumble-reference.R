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
    'GenomicRanges',
    #'biomaRt',
    'BSgenome',
    'BSgenome.Hsapiens.UCSC.hg19',
    'Repitools',
    'ggplot2',
    'VariantAnnotation',
    'doParallel',
    'MASS'
))

# Options ------------------------------------------------------------

# Define options
option_list <- list(
    make_option(c("-i", "--input_folder"), action = "store", type = "character",default = '.',
                help = "folder of count RDS files"),
    make_option(c("-a", "--annotation_folder"), action = "store", type = "character",default = '~/jumble/jumble_annotation',
                help = "folder with gene annotation files"),
    make_option(c("-o", "--output_dir"), action = "store", type = "character",default = '.',
                help = "folder to write reference file"),
    make_option(c("-c", "--cores"), action = "store", type = "character",default = '4',
                help = "cores to use")
)
opt <- parse_args(OptionParser(option_list = option_list))


if (!is.null(opt$cores)) registerDoParallel(cores=as.integer(opt$cores))


# Gene tables ------------------------------------------------------------

# # If BED file present in count table
# target_annotation <- allcounts[[1]]$bed # (NULL if absent)
# 
# # If BED file supplied to this script
# if (!is.null(opt$target_annotation)) target_annotation <- fread(opt$target_annotation)

# Read gene and exon tables
if (dir.exists(opt$annotation_folder)) {
    cancergenes_clinseq <- fread(file = paste(opt$annotation_folder,'cancergenes.txt',sep=.Platform$file.sep))
    allgenes <- fread(file = paste(opt$annotation_folder,'allgenes.txt',sep=.Platform$file.sep))
    allexons <- fread(file = paste(opt$annotation_folder,'allexons.txt',sep=.Platform$file.sep))
} else {
    stop('Annotation folder not found.')
}



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




# Reference object ------------------------------------------------------------
if (wgs) allcounts[[1]]$target_bed_file <- 'wgs'
reference <- allcounts[[1]][c("target_bed_file","chromlength","ranges")]

if (!is.null(allcounts[[1]]$flag)) reference$flag <- allcounts[[1]]$flag

reference$date <- date()
reference$samples <- files
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
saveRDS(reference,paste0(opt$output_dir,'/',
                         str_remove(name,'.*/'),
                         '.reference.RDS'))

print(paste(name,'done.'))
