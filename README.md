# Selkirk

Bioinformatics pipeline for processing 16S rRNA amplicon sequence data.

## Pipeline

Custom parameters stored in `config.yaml`.

### Preprocessing

* Check quality with fastqc and multiqc to compile html report.
* Trim off adapters using cutadapt tool.
* Check quality post trimming with fastqc and multiqc again.

### Filtering

* Remove low quality reads using dada2::filterAndTrim function.
* Check quality post filtering with fastqc and multiqc again.

### Dada2 Workflow

* Generate error model using entire dataset.
* De-replicate and infer sequence variants.
* Remove bimeras, assign taxonomy.
* Track reads throughout processing, print results to table.
