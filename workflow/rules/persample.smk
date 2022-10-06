import pandas as pd

samples = pd.read_table("config/units.tsv")

rule mark_duplicates:
    input:
        bams="mapped/{sample}.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="/dedup/{sample}_DD.bam",
        metrics="dedup/{sample}_DD.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v1.14.0/bio/picard/markduplicates"