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


SAMPLES, = glob_wildcards("results/01_mapping/{samples}.sorted.mkdups.merged.bam")


rule all:
    input:
        expand("results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.freebayes.vcf", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz", samples = SAMPLES),
        expand("results/02_snvs/{samples}.rawsnvs.haplotypeCaller.gvcf.gz", samples = SAMPLES),


rule gatk_HaplotypeCaller_gvcf:
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.gvcf.gz",
    log:
        "logs/gatk_HaplotypeCaller.gvcf.{samples}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller.gvcf.{samples}.tsv"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 128 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 7200 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.4.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'HaplotypeCaller '
        '--create-output-variant-index '
        '-I {input.bam} '
        '-R {input.referenceGenome} '
        '-O {output.gvcf} '
        '-ERC GVCF '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} '


rule gatk_HaplotypeCaller_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        gvcf = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz",
    log:
        "logs/gatk_HaplotypeCaller_vcf.{samples}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_vcf.{samples}.tsv"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 4320 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'module load GATK/4.4.0.0-gimkl-2022a ; '
        'gatk --java-options "-Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}" '
        'HaplotypeCaller '
        '--create-output-variant-index '
        '-I {input.bam} '
        '-R {input.referenceGenome} '
        '-O {output.gvcf} '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} '


rule bcftools_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcf = "results/02_snvs/{samples}.rawsnvs.bcftools.vcf.gz",
    log:
        "logs/bcftools_vcf.{samples}.log"
    benchmark:
        "benchmarks/bcftools_vcf.{samples}.tsv"
    threads: 24
    conda:
        "bcftools"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools mpileup --seed 1953 --threads {threads} --max-depth 500 -q 30 -Q 20 -m 10 -O u -f {input.referenceGenome} {input.bam} | bcftools call --threads {threads} -v -m -O z8 > {output.vcf}

        """


rule freebayes_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.merged.bam",
        referenceGenome = "/nesi/nobackup/agresearch03735/reference/ARS_lic_less_alts.male.pGL632_pX330_Slick_CRISPR_24.fa",
    output:
        vcf = "results/02_snvs/{samples}.rawsnvs.freebayes.vcf",
    log:
        "logs/freebayes_vcf.{samples}.log"
    benchmark:
        "benchmarks/freebayes_vcf.{samples}.tsv"
    threads: 2
    conda:
        "freebayes"
    resources:
        mem_gb = lambda wildcards, attempt: 64 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "milan",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        freebayes --standard-filters --trim-complex-tail --min-coverage 30 -C 10 -F 0.1 -f {input.referenceGenome} {input.bam} > {output.vcf} 

        """