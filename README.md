# Selkirk

Bioinformatics pipeline for processing 16S rRNA amplicon sequence data.

## Pipeline

Custom parameters stored in `config.yaml`.

### Preprocessing

* Check quality with fastqc and multiqc to compile html report.
* Trim off adapters using cutadapt tool.
* Check quality post trimming with fastqc and multiqc again. 

### Dada2 Workflow

* Quality filter using dada2::filterandTrim function
* Learn errors with whole dataset
* De-replicate and infer sequences
* Remove bimeras, assign taxonomy.

### Make tree

* Align sequences with `ssu-align`
* Make tree with `FastTree`
