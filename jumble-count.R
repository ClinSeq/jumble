# Markus Mayrhofer 2022
# Dependencies ------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bamsignals))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))


# Options ------------------------------------------------------------
option_list <- list(
  make_option(c("-t", "--target_bed_file"), action = "store", type = "character",default = NULL, 
              help = "target bed file, unless WGS"),
  make_option(c("-b", "--input_bam_file"), action = "store", type = "character",default = NULL, 
              help = "input BAM file"),
  make_option(c("-o", "--output_RDS_file"), action = "store", type = "character",default = NULL, 
              help = "output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))
# opt <- list(target_bed_file='/Users/markus/Data/CLINSEQ/targets.bed',
#             input_bam_file='/Users/markus/Data/CLINSEQ/sample1.bam',
#             output_RDS_file=NULL)


bampath <- opt$input_bam_file
bamfile <- BamFile(bampath)
chromlength <- seqlengths(seqinfo(bamfile))[c(1:22,'X','Y')]


# BED file ------------------------------------------------------------

wgs <- F
if (is.null(opt$target_bed_file)) wgs <- T

if (!wgs) {   
    
    targets_original <- fread(opt$target_bed_file)
    
    # Let the first column be the target number
    targets_original <- cbind(1:nrow(targets_original),targets_original)
    
    # Set useful column names
    colnames(targets_original) <- c('target','chromosome','start','end')
    
    
}


# WGS bins ------------------------------------------------------------


if (wgs) {
    
    binsize <- 1000
    
    bins <- NULL
    for (i in 1:length(chromlength)) {
        bins[[i]] <- data.table(chromosome=names(chromlength[i]),
                                start=seq(0,chromlength[i],binsize)+1)
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

counts <- opt
counts$date_count <- date()
counts[['ranges']] <- ranges
counts[['chromlength']] <- chromlength

counts[['count']] <- bamCount(bampath, ranges, paired.end="midpoint", 
                                     mapq=20, filteredFlag=1024, verbose=F)
counts[['count_short']] <- bamCount(bampath, ranges, paired.end="midpoint", 
                                          mapq=20, filteredFlag=1024, tlenFilter=c(0,150), verbose=F)


# Outfile ------------------------------------------------------------

out <- opt$output_RDS_file
if (is.null(out)) out <- paste0(bampath,'.counts.RDS')
if (!str_detect(out,'.[Rr][Dd][Ss]$')) out <- paste0(out,'.RDS')

saveRDS(counts,file = out)

