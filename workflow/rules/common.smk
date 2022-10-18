import pandas as pd

samples = list(pd.read_table(config["units"])['Sample_ID'])

# Run up until before HC step with validation
rule preHC:
    input:
        expand("fixRO/{sample}.bam", sample=samples),
        expand("validateBAM/{sample}_fixRO_SUMMARY.txt", sample=samples)

rule dag:
    output:
        "figures/dag.png"
    conda:
        "../envs/cmd.yaml"
    shell:
        "snakemake -n --dag | dot -Tpng > figures/dag.png"