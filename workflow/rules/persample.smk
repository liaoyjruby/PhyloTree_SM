import pandas as pd

configfile: "config/config.yml"
workdir: config['workdir']
wildcard_constraints:
    sample = '|'.join(config["samples"])

def sample2UUID(wildcards):
    units_df = pd.read_table(config["units"])
    UUID_path = units_df[units_df['Sample_ID'] == wildcards.sample].iloc[0]['Mapped']
    return(UUID_path)

rule get_mapped:
    input:
        bam_path = sample2UUID
    output:
        "mapped/{sample}.bam"
    shell:
        "cp {input.bam_path} mapped/{wildcards.sample}.bam"

rule deduplicate:
    input:
        "mapped/{sample}.bam"
    output:
        bam="dedupe/{sample}.bam",
        metrics="dedupe/{sample}.metrics.txt",
    log:
        "logs/dedupe/{sample}.log"
    shell:
        'gatk --java-options ""-Xmx4g"" '
            'MarkDuplicates '
            '-I {input} '
            '-O {output.bam} '
            '-M {output.metrics} '
            '--REMOVE_DUPLICATES true '
            '--CREATE_INDEX true'

rule split_n_cigar_reads:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        bam="dedupe/{sample}.bam"
    output:
        bam="cigar/{sample}.bam"
    log:
        "logs/cigar/{sample}.log"
    shell:
        'gatk SplitNCigarReads '
            '-R {input.ref} '
            '-I {input.bam} '
            '-O {output.bam} '

def getID(sample_ID):
    units_df = pd.read_table(config["units"])
    uID = units_df[units_df['Sample_ID'] == sample_ID].iloc[0]['ID']
    return(str(uID))

rule fix_RG:
    input:
        "cigar/{sample}.bam"
    output:
        temp("fixRG/{sample}.bam")
    log:
        "logs/fixRG/{sample}.log"
    run:
        print(sample)
        uID = getID(wildcards.sample_ID)
        print(uID)
        shell(
            """
            gatk AddOrReplaceReadGroups \
                -I {input} \
                -O {output} \
                -SO "coordinate" \
                -LB "bar" \
                -SM "{wildcards.sample_ID}" \
                -PL "illumina" \
                -PU "ID{uID}" \
                --CREATE_INDEX true
            """)

rule fix_RO:
    input:
        ref_fasta="references/Homo_sapiens_assembly38.fasta",
        ref_dict="references/Homo_sapiens_assembly38.dict",
        bam="fixRG/{sample}.bam"
    output:
        temp("fixRO/{sample}.bam")
    log:
        "logs/fixRO/{sample}.log"
    shell:
        """
        gatk ReorderSam \
            -I {input.bam} \
            -R {input.ref_fasta} \
            -SD {input.ref_dict} \
            -O {output} \
            --CREATE_INDEX true
        """

rule haplotype_caller:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        bam="fixRO/{sample}.bam"
    output:
        protected("hcGVCF/{sample}.g.vcf.gz")
    log:
        "logs/hcGVCF/{sample}.log"
    shell:
        """
        gatk HaplotypeCaller --java-options "-Xmx6g" \
            -I {input.bam} \
            -R {input.ref} \
            -O {output} \
            -ERC "GVCF" \
            -OVI true
        """

rule genotype_GVCFs:
    input:
        ref="references/Homo_sapiens_assembly38.fasta",
        gvcf="hcGVCF/{sample}.g.vcf.gz"
    output:
        "hcVCF/{sample}.vcf.gz"
    log:
        "logs/hcVCF/{sample}.log"
    shell:
        """
        gatk GenotypeGVCFs --java-options "-Xmx4g" \
            -V {input.gvcf} \
            -R {input.ref} \
            -O {output} \
            -OVI true
        """
