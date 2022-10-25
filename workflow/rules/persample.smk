import pandas as pd

def sample2UUID(wildcards):
    units_df = pd.read_table(config["units"])
    UUID_path = units_df[units_df['Sample_ID'] == wildcards.sample].iloc[0]['Mapped_Path']
    return(UUID_path)

rule get_mapped:
    input:
        bam_path = sample2UUID
    output:
        protected("mapped/{sample}.bam")
    shell:
        'cp "{input.bam_path}" mapped/ ; '
        'fname=$(basename "{input.bam_path}") ; '
        'mv "mapped/$fname" "mapped/{wildcards.sample}.bam"'

rule deduplicate:
    input:
        "mapped/{sample}.bam"
    output:
        bam=temp("dedupe/{sample}.bam"),
        bai=temp("dedupe/{sample}.bai"), 
        metrics="dedupe/{sample}.metrics.txt",
    log:
        "logs/dedupe/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options ""-Xmx4g"" '
            'MarkDuplicates '
            '-I {input} '
            '-O {output.bam} '
            '-M {output.metrics} '
            '--CREATE_INDEX true &> {log}'

rule split_n_cigar_reads:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        bam="dedupe/{sample}.bam",
        bai="dedupe/{sample}.bai"
    output:
        bam=temp("cigar/{sample}.bam"),
        bai=temp("cigar/{sample}.bai")
    log:
        "logs/cigar/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk SplitNCigarReads '
            '-R {input.ref} '
            '-I {input.bam} '
            '-O {output.bam} &> {log}'

def getID(wildcards):
    units_df = pd.read_table(config["units"])
    uID = units_df[units_df['Sample_ID'] == wildcards.sample].iloc[0]['ID']
    return(str(uID))

rule fix_RG:
    input:
        bam="cigar/{sample}.bam",
        bai="cigar/{sample}.bai",
    output:
        # temp("fixRG/{sample}.bam")
        bam=temp("fixRG/{sample}.bam"),
        bai=temp("fixRG/{sample}.bai")
    params:
        uID=lambda wildcards: getID(wildcards)
    log:
        "logs/fixRG/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -SO "coordinate" \
            -LB "bar" \
            -SM "{wildcards.sample}" \
            -PL "illumina" \
            -PU "ID{params.uID}" \
            --CREATE_INDEX true  &> {log}
        """

rule fix_RO:
    input:
        ref_fasta="references/Homo_sapiens_assembly38.fasta",
        ref_dict="references/Homo_sapiens_assembly38.dict",
        bam="fixRG/{sample}.bam",
        bai="fixRG/{sample}.bai"
    output:
        # temp("fixRO/{sample}.bam")
        bam=temp("fixRO/{sample}.bam"),
        bai=temp("fixRO/{sample}.bai")
    log:
        "logs/fixRO/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk ReorderSam \
            -I {input.bam} \
            -R {input.ref_fasta} \
            -SD {input.ref_dict} \
            -O {output.bam} \
            --CREATE_INDEX true &> {log}
        """

rule haplotype_caller:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        bam="fixRO/{sample}.bam",
        bai="fixRO/{sample}.bai"
    output:
        gvcf="hcGVCF/{sample}.g.vcf.gz",
        gvcftbi="hcGVCF/{sample}.g.vcf.gz.tbi"
    log:
        "logs/hcGVCF/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller --java-options "-Xmx6g" \
            -I {input.bam} \
            -R {input.ref} \
            -O {output.gvcf} \
            -ERC "GVCF" \
            -OVI true &> {log}
        """

rule genotype_GVCFs:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        gvcf="hcGVCF/{sample}.g.vcf.gz",
        gvcftbi="hcGVCF/{sample}.g.vcf.gz.tbi"
    output:
        vcf="hcVCF/{sample}.vcf.gz",
        vcftbi="hcVCF/{sample}.vcf.gz.tbi"
    log:
        "logs/hcVCF/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk GenotypeGVCFs --java-options "-Xmx4g" \
            -V {input.gvcf} \
            -R {input.ref} \
            -O {output.vcf} \
            -OVI true &> {log}
        """
