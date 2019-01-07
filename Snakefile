# **********************************
# * Snakefile for 16S pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
        "results/multiqc_report_raw.html",
        "results/multiqc_report_cut.html",
        "results/multiqc_report_filt.html",
        "results/merged_sequence_lengths.pdf",
        "results/reads_tracked.txt",
        "results/seqtab_final.rds",
        "results/taxa_final.rds"

rule cutadapt:
    input:
        r1 = "data/rawdata/{sample}_R1_001.fastq.gz",
        r2 = "data/rawdata/{sample}_R2_001.fastq.gz"
    output:
        r1 = "data/cutdata/{sample}_r1_cut.fastq.gz",
        r2 = "data/cutdata/{sample}_r2_cut.fastq.gz"
    conda: "envs/cutadapt_env.yaml"
    shell:
            "cutadapt -e 0 -O 10 -m 50 -g {config[fwd_primer]} "
            "-G {config[rev_primer]} -a {config[rev_primer_rc]} "
            "-A {config[fwd_primer_rc]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2}"

rule fastqc_raw:
    input:
        r1 = "data/rawdata/{sample}_R1_001.fastq.gz",
        r2 = "data/rawdata/{sample}_R2_001.fastq.gz"
    output:
        r1 = "data/rawdata/fastqc/{sample}_R1_001_fastqc.html",
        r2 = "data/rawdata/fastqc/{sample}_R2_001_fastqc.html"
    conda: "envs/fastqc_env.yaml"
    shell: "fastqc -o data/rawdata/fastqc {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand("data/rawdata/fastqc/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        r2 = expand("data/rawdata/fastqc/{sample}_R2_001_fastqc.html", sample=SAMPLES),
    output: "results/multiqc_report_raw.html"
    conda: "envs/multiqc_env.yaml"
    shell: "multiqc -f data/rawdata/fastqc -o results -n multiqc_report_raw.html"

rule fastqc_cut:
    input:
        r1 = "data/cutdata/{sample}_r1_cut.fastq.gz",
        r2 = "data/cutdata/{sample}_r2_cut.fastq.gz"
    output:
        r1 = "data/cutdata/fastqc/{sample}_r1_cut_fastqc.html",
        r2 = "data/cutdata/fastqc/{sample}_r2_cut_fastqc.html"
    conda: "envs/fastqc_env.yaml"
    shell: "fastqc -o data/cutdata/fastqc {input.r1} {input.r2}"

rule multiqc_cut:
    input:
        r1 = expand("data/cutdata/fastqc/{sample}_r1_cut_fastqc.html", sample=SAMPLES),
        r2 = expand("data/cutdata/fastqc/{sample}_r2_cut_fastqc.html", sample=SAMPLES),
    output: "results/multiqc_report_cut.html"
    conda: "envs/multiqc_env.yaml"
    shell: "multiqc -f data/cutdata/fastqc -o results -n multiqc_report_cut.html"

rule fastqc_filt:
    input:
        r1 = "data/filtered/{sample}_r1_filtered.fastq.gz",
        r2 = "data/filtered/{sample}_r2_filtered.fastq.gz"
    output:
        r1 = "data/filtered/fastqc/{sample}_r1_filtered_fastqc.html",
        r2 = "data/filtered/fastqc/{sample}_r2_filtered_fastqc.html"
    conda: "envs/fastqc_env.yaml"
    shell: "fastqc -o data/filtered/fastqc {input.r1} {input.r2}"

rule multiqc_filt:
    input:
        r1 = expand("data/filtered/fastqc/{sample}_r1_filtered_fastqc.html", sample=SAMPLES),
        r2 = expand("data/filtered/fastqc/{sample}_r2_filtered_fastqc.html", sample=SAMPLES),
    output: "results/multiqc_report_filt.html"
    conda: "envs/multiqc_env.yaml"
    shell: "multiqc -f data/filtered/fastqc -o results -n multiqc_report_filt.html"

rule filter:
    input:
        listfiles = config["list_files"],
        r1 = expand("data/cutdata/{sample}_r1_cut.fastq.gz", sample=SAMPLES),
        r2 = expand("data/cutdata/{sample}_r2_cut.fastq.gz", sample=SAMPLES)
    output:
        r1 = expand("data/filtered/{sample}_r1_filtered.fastq.gz", sample=SAMPLES),
        r2 = expand("data/filtered/{sample}_r2_filtered.fastq.gz", sample=SAMPLES),
        filt_out = "results/filt_out.rds"
    conda: "envs/dada2_env.yaml"
    threads: config["num_threads"]
    script: "scripts/filter.R"

rule run_dada2:
    input:
        listfiles = config["list_files"],
        filt_out = "results/filt_out.rds"
    output:
        "results/merged_sequence_lengths.pdf",
        "results/reads_tracked.txt",
        "results/seqtab_final.rds",
        "results/taxa_final.rds"
    conda: "envs/dada2_env.yaml"
    threads: config["num_threads"]
    script: "scripts/run_dada2.R"
