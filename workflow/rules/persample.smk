import pandas as pd

configfile: "config/config.yml"

def sample2UUID(wildcards):
    units_df = pd.read_table("config/units.tsv")
    UUID_path = units_df[units_df['Sample_ID'] == wildcards.sample_ID].iloc[0]['Mapped']
    return(UUID_path)

rule get_mapped:
    input:
        bam_path = sample2UUID
    output:
        "mapped/{sample_ID}.bam"
    run:
        sample_ID = UUID2sample({input})
        shell("echo {input.bam_path}; cp {input.bam_path} /mapped/{sample_ID}", sample_ID=sample_ID)

# rule deduplicate:
#     input:
#         bam=unit_paths
#     output:
#         bam="/dedupe/{sample}_DD.bam",
#         metrics="/dedupe/{sample}_DD.metrics.txt",
#     log:
#         "logs/dedupe/{sample}.log"
#     shell:
#         'gatk --java-options ""-Xmx4g"" '
#             'MarkDuplicates '
#             '-I {input.bam} '
#             '-O {output.bam} '
#             '-M {output.metrics} '
#             '--CREATE_INDEX true'


rule split_n_cigar_reads:
    input:
        
    output:

rule fix_RG:
    input:
    output:

rule fix_RO:
    input:
    output:

rule validate_BAM:
    input:
    output:

rule haplotype_caller:
    input:
    output:

rule genotype_GVCFs:
    input:
    output:

rule validate_VCF:
    input:
    output: