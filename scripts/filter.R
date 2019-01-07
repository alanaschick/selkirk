## Alana Schick, August 2018
## This is a script to filter 16S microbiome data using dada2's function filterAndTrim
## Note that prior to running this script, primers have been trimmed using Cutadapt - resulting files called samplename_cut_r1.fastq.qc

## This is a version of the script to be run using snakemake

library(dada2)
packageVersion("dada2")

## Set variables
list_of_filenames <- snakemake@input$listfiles

## Set filtering parameters from config file
trimleft <- snakemake@config$trimleft
expected_errors <- c(snakemake@config$expected_errors_forward, snakemake@config$expected_errors_reverse)
truncate <- c(snakemake@config$truncate_forward, snakemake@config$truncate_reverse)


#########################################################


## Make a vector of sample names
samples <- scan(list_of_filenames, what = "character")

## file names of forward and reverse reads, before quality filtering
forward_reads <- file.path("data/cutdata", paste0(samples, "_r1_cut.fastq.gz"))
reverse_reads <- file.path("data/cutdata", paste0(samples, "_r2_cut.fastq.gz"))

####### Step 1: Quality filtering

## file names of forward and reverse reads, after quality filtering
filtered_forward_reads <- file.path("data", "filtered", paste0(samples, "_r1_filtered.fastq.gz"))
filtered_reverse_reads <- file.path("data", "filtered", paste0(samples, "_r2_filtered.fastq.gz"))

out <- filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads, filtered_reverse_reads, maxEE = expected_errors, multithread = TRUE, rm.phix=TRUE, truncLen=truncate, trimLeft = trimleft)

saveRDS(out, file.path("results", "filt_out.rds"))
