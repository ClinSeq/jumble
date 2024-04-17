# Jumble Quick-Start Guide

Jumble is a method and R-script(s) for copy number analysis of short read sequencing data, developed by Markus Mayrhofer and Johan Lindberg, Karolinska Institutet. 

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Questions and support: markus.mayrhofer@ki.se.


## Dependencies

To run Jumble you need R ≥4.2, with the Jumble R scripts available and these dependencies installed:

* bamsignals*
* BSgenome*
* BSgenome.Hsapiens.UCSC.hg19*
* data.table
* doParallel
* GenomicRanges*
* ggplot2
* MASS
* patchwork
* optparse
* PSCBS
* Repitools*
* Rsamtools*
* stringr
* VariantAnnotation*

*=Bioconductor

This R code should install all the dependencies. A container is recommended for running Jumble in a pipeline and ours will be made available shortly.

```
# Set a default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# List of CRAN packages
cran_packages <- c("MASS", "ggplot2", 
    "doParallel", "optparse", "patchwork")

# List of Bioconductor packages
bioconductor_packages <- c("BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", 
    "PSCBS", "Repitools", "Rsamtools", "VariantAnnotation", "bamsignals")

# Function to check and install missing packages and report failures
install_if_missing <- function(packages, source) {
    installed <- installed.packages()[, "Package"]
    to_install <- setdiff(packages, installed)
    failed_to_install <- character()

    for (pkg in to_install) {
        result <- tryCatch({
            if (source == "CRAN") {
                install.packages(pkg, dependencies = TRUE, type = "binary", ask = FALSE)
            } else {
                BiocManager::install(pkg, ask = FALSE)
            }
            TRUE
        }, error = function(e) FALSE)

        if (!result) {
            failed_to_install <- c(failed_to_install, pkg)
        }
    }

    failed_to_install
}

# Install missing CRAN packages and capture failures
failed_cran <- install_if_missing(cran_packages, "CRAN")

# Install missing Bioconductor packages and capture failures
failed_bioconductor <- install_if_missing(bioconductor_packages, "Bioconductor")

# Print summary of installation failures
if (length(c(failed_cran, failed_bioconductor)) > 0) {
    cat("Failed to install the following packages:\n")
    cat(paste(c(failed_cran, failed_bioconductor), collapse = "\n"))
} else {
    cat("All packages were successfully installed or were already present.")
}
```

If some dependency fails to install, install them manually. 

## Components

Jumble has three components: 
* **jumble_count.R** 
* **jumble_reference.R**
* **jumble_run.R** (for use in production/pipeline)


Save the scripts wherever convenient.


### jumble_count.R

This script is used with a bam file to create a small RDS file of read counts per target bin. It is not meant for use in "production" (where you would just run query sample bam files), but for speeding up development work, testing, and reference file generation. To run (in the terminal, not in R):

```
Rscript <path>/jumble-count.R -t <target_BED_file> -b <input_bam_file>
```

A "bai" index is also required, and will be expected in the same location as the bam file. The output file name will be based on the input bam file: **<input_bam_file>.counts.RDS**.

Alternatively, if you have all bam (and bai index) files available in a single folder, you can omit the `-b` argument and just run:
```
Rscript <path>/jumble-count.R -t <target_BED_file> -c <cores>
```
This will run the script on all bam files in the folder, and the `-c` option will multithread the work.


If you are just getting started with Jumble and your test dataset, go ahead and run the count script on all your bamfiles, including both the intended reference and query samples.

> Note on input data: Deduplication of aligned sequence data may result in signal saturation with high depth of coverage. Jumble works with deduplicated data, but it is likely better to retain all sequence reads. Marking duplicates and/or UMI processing is typically ok, but the latter has been found to introduce some issues.

### jumble_reference.R

A reference file is required with jumble_run. If you need to create one, collect appropriate <input_bam_file>counts.RDS files in a folder and run:

```
Rscript <path>/jumble-reference.R -i <path_to_counts_files> -a <path_to_annotation_folder> -o <output_path> -c cores
```

The gene annotation folder (-a) is available for download with the scripts and contain (hg19) gene and exon coordinates. You can replace the content with your preferred gene set.

Include "normal-like" samples without evidence of copy number alteration, that are otherwise as similar as possible to the query samples, largely spanning the technical and biological variation you might encounter among the query samples. 

For example, blood-based normal DNA is not the best match for use as reference samples with tissue-extracted or cell-free tumor DNA.

Do not include count files generated using different BED files.

Recommended number of references is 10-100, use at least 4.

The output reference file name will be based on the target BED file that was used for the included counts.RDS files: **<target_BED_file>.reference.RDS**


### jumble_run.R

To run Jumble, you need the query sample BAM file (sorted and indexed, or alternatively just a count file generated as described above) and a reference file based on samples that do not have (significant) copy number alterations. To run:

```
Rscript <path_to_script>/jumble-run.R -r <reference_file> -b <input_bam_file> -o <output_path>
```

Default output path (if omitted) is the current directory. Output file names will be based on the input BAM file: **<input_bam_file>.pdf** etc.

You can include a vcf file of SNPs for a better plot:

```
Rscript <path_to_script>/jumble-run.R -r <reference_file> -b <input_bam_file> -v <vcf_file>
```

> You can run all samples as query samples, including those that are also used as reference samples. This is a good way to verify that each reference sample is normal-like, and detect germline CNVs of interest. Where a query sample encounters itself in the reference set, it is omitted from the reference for that run.



## Results

jumble_run.R creates the following output. 

### Counts

Fragment counts per bin are saved to **<input_bam_file>.counts.RDS**. If there were no indication of copy number alterations in the sample, you can add this file to a future version of your reference pool to further reduce systematic noise. You can also use this file instead of the bam file if you rerun jumble_run.R, saving a few minutes.

### Figures

Figures are written to **<input_bam_file>.png**. 

### Segments

* Full segment details are written to **<input_bam_file>.segments.csv**.
* CNVkit-formatted segments are written to **<input_bam_file>.cns**. 
* DNAcopy-formatted segments are written to **<input_bam_file>.seg**. 

### Target and background bins

* CNVkit-formatted bins are written to **<input_bam_file>.cnr**.
* Full bin table is written to **<input_bam_file>.jumble.RDS**.
* SNPs are written to **<input_bam_file>.jumble.snps.RDS**.