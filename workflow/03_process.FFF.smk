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


SAMPLES = ('1941', '1927', '1936', '1938')


# ENTRY POINTS
        # "results/03_merged/merged.FFF.chrom.norm.bcftools.vcf.gz",
        # "results/03_merged/merged.FFF.chrom.norm.freebayes.vcf.gz",
        # "results/03_merged/merged.FFF.chrom.norm.haplotypeCaller.vcf.gz",


rule all:
    input:
        "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
        "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.intersection.vcf.gz",
        "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.intersection.vcf.gz",
        expand("results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz", samples = SAMPLES),
        expand("results/00_stats/{samples}.mosdepth.summary.txt", samples = SAMPLES),
        expand("results/00_stats/{samples}.sorted.mkdups.merged.bam.samtools-stats.txt", samples = SAMPLES),
        expand("results/00_stats/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz.bcftools-stats.txt", samples = SAMPLES),



rule get_eva_snvs:
    output:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
    params:
        vcf = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/ovis_aries/ARSUI_Ramb_v2.0/9940_GCA_016772045.1_current_ids.vcf.gz",
        index = "https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/ovis_aries/ARSUI_Ramb_v2.0/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
    threads: 2
    shell:
        """
        
        wget -t 5 -O {output.vcf} {params.vcf} ;
        wget -t 5 -O {output.csi} {params.index} ; 
        
        """


rule rename_eva_snvs:
    input:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.vcf.gz.csi",
        sed_file = "resources/eva_rename.sed",
    output:
        vcf = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi",
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        zcat {input.vcf} | sed -f {input.sed_file} | bcftools view --threads {threads} -O z8 -o {output.vcf} -

        bcftools index --threads {threads} {output.vcf} -o {output.csi} ; 
        
        """


# BCFTOOLS
rule filter_monomorphic_bcftools: 
    priority:100
    input:
        norm = "results/03_merged/merged.FFF.chrom.norm.bcftools.vcf.gz",
        csi = "results/03_merged/merged.FFF.chrom.norm.bcftools.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.bcftools.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.bcftools.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''

        zcat {input.norm} | grep -P -v "0\/1\S+\t0\/1\S+\t0\/1\S+\t0\/1\S+" | grep -P -v "1\/1\S+\t1\/1\S+\t1\/1\S+\t1\/1\S+" |  bcftools view - -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        '''


rule filter_DP_bcftools:
    priority:100
    input:
        norm = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.bcftools.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.bcftools.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.bcftools.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.bcftools.vcf.gz.csi"),
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
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/DP<66 || INFO/DP>282' {input.norm} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule isec_bcftools_eva:
    priority:100
    input:
        filtered = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.bcftools.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.bcftools.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_filtered/isec_bcftools_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_filtered/isec_bcftools_eva/0001.vcf -o {output.filtered} ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi} ; 

        rm -r results/04_filtered/isec_bcftools_eva ;

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


# FREEBAYES
rule filter_monomorphic_freebayes: 
    priority:100
    input:
        norm = "results/03_merged/merged.FFF.chrom.norm.freebayes.vcf.gz",
        csi = "results/03_merged/merged.FFF.chrom.norm.freebayes.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.freebayes.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.freebayes.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''

        zcat {input.norm} | grep -P -v "0\/1\S+\t0\/1\S+\t0\/1\S+\t0\/1\S+" | grep -P -v "1\/1\S+\t1\/1\S+\t1\/1\S+\t1\/1\S+" |  bcftools view - -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        '''


rule filter_DP_freebayes:
    priority:100
    input:
        norm = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.freebayes.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.freebayes.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.freebayes.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.freebayes.vcf.gz.csi"),
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
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/DP<66 || INFO/DP>282' {input.norm} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule isec_freebayes_eva:
    priority:100
    input:
        filtered = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.freebayes.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.freebayes.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_filtered/isec_freebayes_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_filtered/isec_freebayes_eva/0001.vcf -o {output.filtered} ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi} ; 

        rm -r results/04_filtered/isec_freebayes_eva ;

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


# GATK
rule filter_monomorphic_haplotypeCaller: 
    priority:100
    input:
        norm = "results/03_merged/merged.FFF.chrom.norm.haplotypeCaller.vcf.gz",
        csi = "results/03_merged/merged.FFF.chrom.norm.haplotypeCaller.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.haplotypeCaller.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.haplotypeCaller.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 120),
        partition = "compute",
        DTMP = "tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        '''

        zcat {input.norm} | grep -P -v "0\/1\S+\t0\/1\S+\t0\/1\S+\t0\/1\S+" | grep -P -v "1\/1\S+\t1\/1\S+\t1\/1\S+\t1\/1\S+" |  bcftools view - -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        '''


rule filter_DP_haplotypeCaller:
    priority:100
    input:
        norm = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.haplotypeCaller.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.haplotypeCaller.vcf.gz.csi",
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.haplotypeCaller.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.haplotypeCaller.vcf.gz.csi"),
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
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/DP<66 || INFO/DP>282' {input.norm} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule isec_haplotypeCaller_eva:
    priority:100
    input:
        filtered = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.haplotypeCaller.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.haplotypeCaller.vcf.gz.csi",
        eva = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz",
        eva_csi = "resources/eva/9940_GCA_016772045.1_current_ids.sed.vcf.gz.csi"
    output:
        filtered = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.vcf.gz"),
        csi = temp("results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.vcf.gz.csi"),
    threads:8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        bcftools isec -O v -p results/04_filtered/isec_haplotypeCaller_eva {input.eva} {input.filtered} ;

        bcftools view --threads {threads} -O z8 results/04_filtered/isec_haplotypeCaller_eva/0001.vcf -o {output.filtered} ;

        bcftools index --threads {threads} {output.filtered} -o {output.csi} ; 

        rm -r results/04_filtered/isec_haplotypeCaller_eva ;

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt; 

        """


# INTERSECTION
rule ensemble_intersection:
    priority: 100
    input:
        bcftools = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.vcf.gz",
        bcftools_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.vcf.gz.csi",
        freebayes = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.vcf.gz",
        freebayes_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.vcf.gz.csi",
        haplotypeCaller = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.vcf.gz",
        haplotypeCaller_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.vcf.gz.csi",
    output:
        bcftools_common = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
        bcftools_common_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.csi",
        freebayes_common = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.intersection.vcf.gz",
        freebayes_common_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.freebayes.intersection.vcf.gz.csi",
        haplotypeCaller_common = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.intersection.vcf.gz",
        haplotypeCaller_common_csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.haplotypeCaller.intersection.vcf.gz.csi",
    threads:6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """

        bcftools isec -O z8 -p results/04_filtered -n~111 --threads {threads} {input.bcftools} {input.freebayes} {input.haplotypeCaller};

        mv results/04_filtered/0000.vcf.gz {output.bcftools_common};
        echo "Total snps in {output.bcftools_common}: $(bcftools view --threads {threads} {output.bcftools_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/04_filtered/0001.vcf.gz {output.freebayes_common};
        echo "Total snps in {output.freebayes_common}: $(bcftools view --threads {threads} {output.freebayes_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        mv results/04_filtered/0002.vcf.gz {output.haplotypeCaller_common};
        echo "Total snps in {output.haplotypeCaller_common}: $(bcftools view --threads {threads} {output.haplotypeCaller_common} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        bcftools index --threads {threads} {output.bcftools_common};

        bcftools index --threads {threads} {output.freebayes_common};

        bcftools index --threads {threads} {output.haplotypeCaller_common};

        rm results/04_filtered/0000.vcf.gz.tbi
        rm results/04_filtered/0001.vcf.gz.tbi
        rm results/04_filtered/0002.vcf.gz.tbi

        """


rule bcftools_private_snps:
    priority:100
    input:
        merged = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
        csi = "results/04_filtered/merged.FFF.chrom.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.csi",
    output:
        private = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
        csi = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.csi",
    threads:6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """ 
        sleep 5;

        bcftools view -O z8 --samples {wildcards.samples} --private --threads {threads} {input.merged} -o {output.private};

        bcftools index --threads {threads} {output.private} -o {output.csi};

        echo "Total snps in {output.private}: $(bcftools view --threads {threads} {output.private} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule filter_MQ60_bcftools:
    priority:100
    input:
        private = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz",
        csi = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.vcf.gz.csi",
    output:
        filtered = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz",
        csi = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz.csi",
    threads: 8
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 8 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 60 + ((attempt - 1) * 60),
        partition = "compute",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        """
        # -e is 'exclude'

        bcftools view --threads {threads} -e 'INFO/MQ<60' {input.private} -O z8 -o {output.filtered};

        bcftools index --threads {threads} {output.filtered} -o {output.csi};

        echo "Total snps in {output.filtered}: $(bcftools view --threads {threads} {output.filtered} | grep -v "#" | wc -l)" | tee -a snps.counts.summary.txt;

        """


rule bcftools_stats_private_MQ:
    priority: 100
    input:
        filtered = "results/05_private/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz",
    output:
        stats = "results/00_stats/{samples}.FFF.norm.monomorphic.DPFilt.eva.bcftools.intersection.MQ60.vcf.gz.bcftools-stats.txt",
    threads: 6
    conda:
        "bcftools-1.19"
    resources:
        mem_gb = lambda wildcards, attempt: 12 + ((attempt - 1) * 64),
        time = lambda wildcards, attempt: 120 + ((attempt - 1) * 60),
        partition = "compute",
        DTMP = "/nesi/nobackup/agresearch03735/SMK-SNVS/tmp",
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "bcftools stats "
        "--fasta-ref resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna "
        "--samples - "
        "--threads 6 "
        "{input.filtered} > "
        "{output.stats} "


rule samtools_stats_merged:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
        referenceGenome = "resources/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna"
    output:
        stats = "results/00_stats/{samples}.sorted.mkdups.merged.bam.samtools-stats.txt"
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
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "samtools stats --threads {threads} -r {input.referenceGenome} {input.bam} > {output.stats} "


rule mosdepth_stats_merged:
    priority: 100
    input:
        bam = "results/01_mapping/{samples}.sorted.mkdups.bam",
    output:
        stats = "results/00_stats/{samples}.mosdepth.summary.txt",
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
        attempt = lambda wildcards, attempt: attempt,
    shell:
        "mosdepth "
        "--use-median "
        "--fast-mode " # dont look at internal cigar operations or correct mate overlaps (recommended for most use-cases).
        "--no-per-base " # dont output per-base depth.
        "--threads {threads} "
        "results/00_stats/{wildcards.samples} " # output prefix
        "{input.bam} "

