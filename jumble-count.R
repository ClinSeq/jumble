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
    'stringr',
    'doParallel',
    'data.table',
    'bamsignals',
    'GenomicRanges',
    'Rsamtools',
    'csaw'
))

# Options ------------------------------------------------------------
option_list <- list(
    make_option(c("-t", "--target_bed_file"), action = "store", type = "character",default = NULL,
                help = "target bed file, unless WGS"),
    make_option(c("-b", "--input_bam_file"), action = "store", type = "character",default = NULL,
                help = "input BAM file, else all in current directory"),
    make_option(c("-o", "--output_dir"), action = "store", type = "character",default = '.',
                help = "output directory, else in current directory"),
    make_option(c("-s", "--wgs_bin_size"), action = "store", type = "integer",default = 10000,
                help = "bin size (WGS only)"),
    make_option(c("-c", "--cores"), action = "store", type = "integer",default = 1,
                help = "threads to use")
)

opt <- parse_args(OptionParser(option_list = option_list))
# opt <- list(target_bed_file='tests/testthat/count/testbed.bed',
#             input_bam_file='tests/testthat/count/sample1_1e6.bam')

registerDoParallel(cores = opt$cores)


if (!is.null(opt$input_bam_file)) {
    bamfiles <- opt$input_bam_file
} else {
    bamfiles <- dir(pattern = '.bam$',ignore.case = T)
}

foreach(i=1:length(bamfiles)) %dopar% {

    opt$input_bam_file <- bamfiles[i]
    bampath <- opt$input_bam_file
    bamfile <- BamFile(bampath)

    bam_is_ucsc <- any(str_detect(seqnames(seqinfo(bamfile)),'^chr'))
    if (bam_is_ucsc) {
        chromlength <- seqlengths(seqinfo(bamfile))[paste0('chr',c(1:22,'X','Y'))]

        # Here, if UCSC style, chromosome names are set to 1,2,3,...
        names(chromlength) <- str_remove(names(chromlength),'^chr')

    } else chromlength <- seqlengths(seqinfo(bamfile))[c(1:22,'X','Y')]


    # BED file ------------------------------------------------------------

    wgs <- F
    if (is.null(opt$target_bed_file)) wgs <- T

    bed <- NULL
    if (!wgs) {

        bed <- fread(opt$target_bed_file)

        # Let the first column be the target number
        targets_original <- cbind(1:nrow(bed),bed[,1:3])

        # Set useful column names
        colnames(targets_original) <- c('target','chromosome','start','end')

        # Make sure targets are EMBL style
        targets_original$chromosome <- str_remove(targets_original$chromosome,'^chr')

    }


    # WGS bins ------------------------------------------------------------


    if (wgs) {

        binsize <- opt$wgs_bin_size

        bins <- NULL
        for (i in 1:length(chromlength)) {
            bins[[i]] <- data.table(chromosome=names(chromlength[i]),
                                    start=seq(0,chromlength[i]-binsize*2,binsize)+1)
        }
        bins <- rbindlist(bins)
        bins[,end:=start+binsize-1]

    }



    # Target bins ------------------------------------------------------------

    if (!wgs) {

        binsize <- 200

        ranges <- makeGRangesFromDataFrame(targets_original)

        # Extend targets by 50
        start(ranges) <- start(ranges)-50
        end(ranges) <- end(ranges)+50

        # Extend targets by 50 again if any is below 200.
        width <- end(ranges)-start(ranges)
        if (any(width<200)) {
            start(ranges) <- start(ranges)-50
            end(ranges) <- end(ranges)+50
        }

        # Merge
        ranges <- reduce(ranges)

        # Shrink to multiple of binsize
        width <- end(ranges)-start(ranges)
        mid <- start(ranges)+round(width/2)
        new_width <- width - width %% binsize
        new_start <- mid - new_width/2
        new_end <- mid + new_width/2
        chrom <- as.character(seqnames(ranges))
        n_bins <- new_width/binsize


        # split longer targets
        targets_list <- NULL
        for (i in 1:length(chrom)) {
            bins <- data.table(chromosome=chrom[i],
                               start=seq(new_start[i],new_end[i]-binsize,binsize)+1,
                               end=seq(new_start[i]+binsize,new_end[i],binsize))
            bins[,mid:=ceiling((start+end)/2)]
            targets_list[[i]] <- bins
        }
        targets <- rbindlist(targets_list)
        targets[,length:=end-start]

        target_ranges <- makeGRangesFromDataFrame(targets)

    }

    # Background bins ------------------------------------------------------------

    if (!wgs) {

        # Add a dummy target at start and end of each chromosome to make sure it is included
        dummy_targets <- rbind(data.table(target=0,chromosome=c(1:22,'X','Y'),start=1,end=10),
                               data.table(target=0,chromosome=c(1:22,'X','Y'),start=chromlength-10,end=chromlength)
        )


        # Make background ranges object
        binsize <- 1e6
        minsize <- 3e5
        ranges <- makeGRangesFromDataFrame(rbind(targets_original[,chromosome:=as.character(chromosome)],dummy_targets)[order(chromosome,start)])

        # Extend original targets by 1k for good margin
        start(ranges) <- start(ranges)-1000
        end(ranges) <- end(ranges)+1000

        # Merge
        ranges <- reduce(ranges)

        # Invert
        ranges <- gaps(ranges)

        # Drop too short
        ranges <- ranges[width(ranges)>minsize]

        # Shrink to multiple of binsize
        width <- end(ranges)-start(ranges)
        mid <- start(ranges)+round(width/2)
        new_width <- width
        new_width[width > binsize] <- width[width > binsize] - width[width > binsize] %% binsize
        new_start <- mid - round(new_width/2)
        new_end <- mid + round(new_width/2)
        chrom <- as.character(seqnames(ranges))
        n_bins <- ceiling(new_width/binsize)

        # split longer antitargets
        antitargets_list <- NULL
        for (i in 1:length(chrom)) {
            bins <- data.table(chromosome=chrom[i],start=new_start[i],end=new_end[i])
            if (n_bins[i]>1) {
                bins <- data.table(chromosome=chrom[i],
                                   start=seq(new_start[i],new_end[i]-binsize,binsize)+1,
                                   end=seq(new_start[i]+binsize,new_end[i],binsize))
            } else {
                bins[,start:=start+1]
            }
            bins[,mid:=(start+end)/2]
            antitargets_list[[i]] <- bins
        }
        background <- rbindlist(antitargets_list)
        background[,length:=end-start]

        background_ranges <- makeGRangesFromDataFrame(background)

        overlap <- findOverlaps(target_ranges,background_ranges)
        if (length(overlap)>0) cat ('Warning! Unexpected overlap')


    }

    # Make counts list ------------------------------------------------------------

    if (wgs) ranges <- makeGRangesFromDataFrame(bins)

    if (!wgs) ranges <- sort(makeGRangesFromDataFrame(rbind(targets,background)))


    # If bam file is UCSC, make sure the ranges object is a match
    if (bam_is_ucsc)
        seqlevelsStyle(ranges) <- 'UCSC'

    counts <- opt
    counts$date_count <- date()
    counts$bed <- bed
    counts[['ranges']] <- ranges
    counts[['chromlength']] <- chromlength

    flag <- 2816
    counts$flag <- flag

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


    # Outfile ------------------------------------------------------------


    out <- paste0(basename(bampath),'.counts.RDS')

    saveRDS(counts,file = paste0(opt$output_dir,'/',out))

}
