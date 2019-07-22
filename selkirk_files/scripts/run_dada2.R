## Alana Schick, August 2018
## This is a script to process 16S microbiome data using dada2
## Note that prior to running this script, primers have been trimmed using Cutadapt and reads have been filtered

## This is a version of the script to be run using snakemake

library(dada2)
packageVersion("dada2")

## Set variables
list_of_filenames <- snakemake@input$listfiles

## Make a vector of sample names
samples <- scan(list_of_filenames, what = "character")

out <- readRDS(snakemake@input$filt_out)


## file names of forward and reverse reads, after quality filtering
filtered_forward_reads <- file.path("data", "filtered", paste0(samples, "_r1_filtered.fastq.gz"))
filtered_reverse_reads <- file.path("data", "filtered", paste0(samples, "_r2_filtered.fastq.gz"))


####### Step 2: Generate error model of data

err_forward <- learnErrors(filtered_forward_reads, multithread = TRUE)
err_reverse <- learnErrors(filtered_reverse_reads, multithread = TRUE) 



####### Step 3: Derepliate sequences

derep_forward <- derepFastq(filtered_forward_reads, verbose = TRUE)
derep_reverse <- derepFastq(filtered_reverse_reads, verbose = TRUE)

names(derep_forward) <- samples
names(derep_reverse) <- samples



####### Step 4: Infer sequence variants

dadaF <- dada(derep_forward, err = err_forward, multithread = TRUE)
dadaR <- dada(derep_reverse, err = err_reverse, multithread = TRUE)



####### Step 5: Merge forward and reverse reads

merged <- mergePairs(dadaF, derep_forward, dadaR, derep_reverse, verbose = TRUE)



####### Step 6: Generate count table

seqtab <- makeSequenceTable(merged)




####### Step 7: Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)


## How are the reads concentrated in the merged sequence lengths
readsbyseqlen <- tapply(colSums(seqtab.nochim), nchar(colnames(seqtab.nochim)),sum)
pdf(file.path("results", "merged_sequence_lengths.pdf"))
plot(as.integer(names(readsbyseqlen)),readsbyseqlen, xlab = "Merged length", ylab = "Total reads")
dev.off()


####### Step 7b: Track reads throughout processing

getN <- function(x) sum (getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=out[,1], filtered=out[,2], dada_f=sapply(dadaF, getN), dada_r=sapply(dadaR, getN), merged=sapply(merged, getN), nonchim=rowSums(seqtab.nochim), total_perc_reads = round(rowSums(seqtab.nochim)/out[,1]*100,1))

## Write this table to output
write.table(summary_tab, file.path("results", "reads_tracked.txt"))


####### Step 8: Assign Taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/home/aschick/refs/dada2/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies (taxa, "/home/aschick/refs/dada2/silva_species_assignment_v132.fa.gz")




####### Step 9: Save output

saveRDS(seqtab.nochim, file.path("results", "seqtab_final.rds"))
saveRDS(taxa, file.path("results", "taxa_final.rds"))




