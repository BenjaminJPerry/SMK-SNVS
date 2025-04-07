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


SAMPLES = ('OFF3DS', '1945')
CHROM = ('NC_056054.1', 'NC_056055.1', 'NC_056056.1', 'NC_056057.1', 'NC_056058.1', 'NC_056059.1', 'NC_056060.1', 'NC_056061.1', 'NC_056062.1', 'NC_056063.1', 'NC_056064.1', 'NC_056065.1', 'NC_056066.1', 'NC_056067.1', 'NC_056068.1', 'NC_056069.1', 'NC_056070.1', 'NC_056071.1', 'NC_056072.1', 'NC_056073.1', 'NC_056074.1', 'NC_056075.1', 'NC_056076.1', 'NC_056077.1', 'NC_056078.1', 'NC_056079.1', 'NC_056080.1')


wildcard_constraints:
    chromosome = "\w+.1"

rule all:
    input:
        "results/03_merged/merged.MFF.chrom.norm.haplotypeCaller.vcf.gz",
        "results/03_merged/merged.MFF.chrom.norm.haplotypeCaller.vcf.gz.csi"


rule gatk_HaplotypeCaller_vcf:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna",
    output:
        vcf_chrom = temp("results/02_snvs/{samples}.{chromosome}.rawsnvs.haplotypeCaller.vcf.gz"),
        csi = temp("results/02_snvs/{samples}.{chromosome}.rawsnvs.haplotypeCaller.vcf.gz.csi")
    params:
        chromosome = '{chromosome}'
    log:
        "logs/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.log"
    benchmark:
        "benchmarks/gatk_HaplotypeCaller_vcf.{samples}.{chromosome}.tsv"
    conda:
        "gatk-4.5.0.0"
    threads: 4
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 12),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        'export PATH=/home/agresearch.co.nz/perrybe/.conda/envs/gatk-4.5.0.0/bin:$PATH; '
        'gatk --java-options "-Xmx{resources.mem_gb}G -XX:ParallelGCThreads={threads}"  '
        'HaplotypeCaller '
        '--pileup-detection '
        '--min-base-quality-score 10 ' # default 10
        '--minimum-mapping-quality 20 ' # default 20
        '--standard-min-confidence-threshold-for-calling 10 ' # default 30; default 10 pre v4.1
        '--indel-size-to-eliminate-in-ref-model 20 ' #default 10
        '--create-output-variant-index '
        '-R {input.referenceGenome} '
        '-I {input.bam} '
        '-L {params.chromosome} '
        '-O {output.vcf_chrom} '
        '--tmp-dir {resources.DTMP} '
        '&> {log}.attempt.{resources.attempt} && '
        'rm results/02_snvs/{wildcards.samples}.{wildcards.chromosome}.rawsnvs.haplotypeCaller.vcf.gz.tbi; '
        'bcftools index --threads {threads} {output.vcf_chrom} -o {output.csi}'



rule concatenate_replicons_vcf: #TODO
    priority:1000
    input:
        vcfgz = expand("results/02_snvs/{samples}.{chromosome}.rawsnvs.haplotypeCaller.vcf.gz", chromosome = CHROM, allow_missing = True),
        csi = expand("results/02_snvs/{samples}.{chromosome}.rawsnvs.haplotypeCaller.vcf.gz", chromosome = CHROM, allow_missing = True),
    output:
        merged = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi",
    benchmark:
        "benchmarks/{samples}_concatenate_replicons_vcf.tsv"
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        
        bcftools concat --threads {threads} {input.vcfgz} | bcftools sort -O z8 -o {output.merged} - ;

        bcftools index --threads {threads} {output.merged} -o {output.csi}

        """


rule view_haplotype_chrom:
    priority:100
    input:
        vcf = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz",
        csi = "results/02_snvs/{samples}.rawsnvs.haplotypeCaller.vcf.gz.csi",
    output:
        filtered_vcf = temp("results/03_merged/{samples}.chrom.haplotypeCaller.vcf.gz",),
        filtered_vcf_csi = temp("results/03_merged/{samples}.chrom.haplotypeCaller.vcf.gz.csi"),
    params:
        chromosomes = "NC_056054.1,NC_056055.1,NC_056056.1,NC_056057.1,NC_056058.1,NC_056059.1,NC_056060.1,NC_056061.1,NC_056062.1,NC_056063.1,NC_056064.1,NC_056065.1,NC_056066.1,NC_056067.1,NC_056068.1,NC_056069.1,NC_056070.1,NC_056071.1,NC_056072.1,NC_056073.1,NC_056074.1,NC_056075.1,NC_056076.1,NC_056077.1,NC_056078.1,NC_056079.1,NC_056080.1"

    benchmark:
        "benchmarks/view_haplotypeCaller_chrom.{samples}.tsv"
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools view {input.vcf} -O z8 -o {output.filtered_vcf} --regions {params.chromosomes};

        bcftools index --threads {threads} {output.filtered_vcf} -o {output.filtered_vcf_csi}; 

        echo "Total snps in {output.filtered_vcf}: $(cat {output.filtered_vcf} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 
        
        """


rule bcftools_norm_samples:
    priority: 100
    input:
        unnormal = "results/03_merged/{samples}.chrom.haplotypeCaller.vcf.gz",
        filtered_vcf_csi = "results/03_merged/{samples}.chrom.haplotypeCaller.vcf.gz.csi"
    output:
        norm = temp("results/03_merged/{samples}.chrom.norm.haplotypeCaller.vcf.gz"),
        csi = temp("results/03_merged/{samples}.chrom.norm.haplotypeCaller.vcf.gz.csi"),
    threads:6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools norm --threads {threads} -O z8 -m+ -s -D -f resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -o {output.norm} {input.unnormal};
    
        echo "Total snps in {output.norm}: $(bcftools view --threads {threads} {output.norm} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        bcftools index --threads {threads} {output.norm};
    
        """


# rule filter_DP:
#     priority:100
#     input:
#         norm = "results/03_filtered/{samples}.chrom.norm.haplotypeCaller.vcf.gz",
#         csi = "results/03_filtered/{samples}.chrom.norm.haplotypeCaller.vcf.gz.csi",
#     output:
#         filtered = temp("results/03_filtered/{samples}.chrom.norm.DPFilt.haplotypeCaller.vcf.gz"),
#         csi = temp("results/03_filtered/{samples}.chrom.norm.DPFilt.haplotypeCaller.vcf.gz.csi"),
#     threads: 8
#     conda:
#         "bcftools-1.19"
#     resources:
#         mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
#         time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
#         partition = "compute",
#         DTMP = "tmp",
#         attempt = lambda wildcards, attempt: attempt,
#     shell:
#         """
#         # -e is 'exclude'

#         bcftools view --threads {threads} -e 'INFO/DP<2 || INFO/DP>2500' {input.norm} -O z8 -o {output.filtered};

#         bcftools index --threads {threads} {output.filtered} -o {output.csi};

#         echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

#         """


# rule filter_QUAL60: 
#     priority:100
#     input:
#         filtered = "results/03_filtered/{samples}.chrom.norm.DPFilt.haplotypeCaller.vcf.gz",
#         csi = "results/03_filtered/{samples}.chrom.norm.DPFilt.haplotypeCaller.vcf.gz.csi",
#     output:
#         filtered = temp("results/03_filtered/{samples}.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz"),
#         csi = temp("results/03_filtered/{samples}.chrom.norm.DPFilt.QUAL60.haplotypeCaller.vcf.gz.csi"),
#     threads:8
#     conda:
#         "bcftools-1.19"
#     resources:
#         mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
#         time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
#         partition = "compute",
#         DTMP = "tmp",
#         attempt = lambda wildcards, attempt: attempt,
#     shell:
#         '''
#         # -e is 'exclude'

#         bcftools view -e 'QUAL<60' {input.filtered} -O z8 -o {output.filtered};

#         bcftools index --threads {threads} {output.filtered} -o {output.csi};

#         echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

#         '''


rule merge_animals_vcf:
    priority:100
    input:
        vcfgz = expand("results/03_merged/{samples}.chrom.norm.haplotypeCaller.vcf.gz", samples = SAMPLES),
        csi = expand("results/03_merged/{samples}.chrom.norm.haplotypeCaller.vcf.gz.csi", samples = SAMPLES),
    output:
        merged = "results/03_merged/merged.MFF.chrom.norm.haplotypeCaller.vcf.gz",
        csi = "results/03_merged/merged.MFF.chrom.norm.haplotypeCaller.vcf.gz.csi",
    benchmark:
        "benchmarks/merge_animals_vcf.tsv"
    threads: 16
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools merge --threads {threads} {input.vcfgz} -O z8 -o {output.merged};

        bcftools index --threads {threads} {output.merged} -o {output.csi};

        echo "Total snps in {output.merged}: $(cat {output.merged} | gunzip | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt 

        """
