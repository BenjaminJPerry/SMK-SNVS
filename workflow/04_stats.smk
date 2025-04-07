# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

#configfile: "config/config.yaml"


import os


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


SAMPLES, = glob_wildcards("results/01_mapping/{samples}.sorted.mkdups.bam")

rule all:
    input:
        expand("results/00_stats/samtools/{samples}.sorted.mkdups.bam.samtools_stats.txt", samples = SAMPLES),
        expand("results/00_stats/mosdepth/{samples}.mosdepth.summary.txt", samples = SAMPLES),
        "results/00_stats/bcftools/merged.MFF.chrom.norm.bcftools.vcf.gz.bcftools-stats.txt",
        "results/00_stats/bcftools/merged.MFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.bcftools-stats.txt",
        "results/00_stats/bcftools/merged.FFF.chrom.norm.bcftools.vcf.gz.bcftools-stats.txt",
        "results/00_stats/bcftools/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.bcftools-stats.txt",
        "results/00_stats/bcftools/OFF3.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz.bcftools-stats.txt",
        "results/00_stats/bcftools/1945.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz.bcftools-stats.txt",


rule samtools_stats_merged:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna"
    output:
        stats = "results/00_stats/samtools/{samples}.sorted.mkdups.bam.samtools_stats.txt"
    log:
        "logs/samtools_stats_merged.{samples}.log"
    benchmark:
        "benchmarks/samtools_stats_merged.{samples}.tsv"
    conda:
        "samtools-1.17"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "samtools stats --threads {threads} -r {input.referenceGenome} {input.bam} > {output.stats} "


rule mosdepth_stats:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
    output:
        stats = "results/00_stats/mosdepth/{samples}.mosdepth.summary.txt",
    log:
        "logs/mosdepth_stats_merged.{samples}.log"
    benchmark:
        "benchmarks/mosdepth_stats_merged.{samples}.tsv"
    conda:
        "mosdepth-0.3.6"
    threads: 12
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "mosdepth "
        "--use-median "
        "--fast-mode " # dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases).
        "--no-per-base " # dont output per-base depth.
        "--threads {threads} "
        "results/00_stats/mosdepth/{wildcards.samples} " # output prefix
        "{input.bam} "


rule bcftools_stats_merged_MFF:
    priority: 100
    input:
        vcf = "results/03_merged/merged.MFF.chrom.norm.bcftools.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/merged.MFF.chrom.norm.bcftools.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "


rule bcftools_stats_int_MFF:
    priority: 100
    input:
        vcf = "results/04_filtered/merged.MFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/merged.MFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "


rule bcftools_stats_merged_FFF:
    priority: 100
    input:
        vcf = "results/03_merged/merged.FFF.chrom.norm.bcftools.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/merged.FFF.chrom.norm.bcftools.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "


rule bcftools_stats_int_FFF:
    priority: 100
    input:
        vcf = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "


rule bcftools_stats_private_OFF3:
    priority: 100
    input:
        vcf = "results/05_private/OFF3.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/OFF3.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "


rule bcftools_stats_private_1945:
    priority: 100
    input:
        vcf = "results/05_private/1945.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz",
    output:
        stats = "results/00_stats/bcftools/1945.MFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.DP.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.vcf} > "
        "{output.stats} "

