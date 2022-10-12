workdir: config['dir']

rule validate_BAM:
    input:
        bam = "{step}/{sample}.bam",
        ref = "references/Homo_sapiens_assembly38.fasta",
    output:
        summary = "validateBAM/{sample}_{step}_summary.txt",
        errors = "validateBAM/{sample}_{step}_errors.txt",
        warnings = "validateBAM/{sample}_{step}_warnings.txt"
    log:
        "logs/validateBAM/{sample}_{step}.log"
    shell:
        """
        gatk ValidateSamFile \
            -I {input.bam} \
            -R {input.ref} \
            -M SUMMARY \
            -O {output.summary}

        gatk ValidateSamFile \
            -I {input.bam} \
            -R {input.ref} \
            -M VERBOSE \
            --IGNORE_WARNINGS true \
            -O {output.errors}
    
        gatk ValidateSamFile \
            -I {input.bam} \
            -R {input.ref} \
            -M VERBOSE \
            --IGNORE MISSING_TAG_NM \
            -O {output.warnings}
        """

rule validate_GVCF:
    input:
        gvcf = "{step}/{sample}.g.vcf.gz",
        ref_fasta = "references/Homo_sapiens_assembly38.fasta",
        ref_dbsnp = "references/Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        "validateGVCF/{sample}_{step}.txt",
    log:
        "logs/validateGVCF/{sample}_{step}.log"
    shell:
        """
        gatk ValidateVariants \
        -V {input.gvcf} \
        -R {input.ref_fasta} \
        -D {input.ref_dbsnp} \
        -Xtype ALLELES \
        -verbosity ERROR \
        -gvcf > {output}
        """

rule validate_VCF:
    input:
        gvcf = "{step}/{sample}.vcf.gz",
        ref_fasta = "references/Homo_sapiens_assembly38.fasta",
        ref_dbsnp = "references/Homo_sapiens_assembly38.dbsnp138.vcf"
    output:
        "validateVCF/{sample}_{step}.txt",
    log:
        "logs/validateVCF/{sample}_{step}.log"
    shell:
        """
        gatk ValidateVariants \
        -V {input.gvcf} \
        -R {input.ref_fasta} \
        -D {input.ref_dbsnp} \
        -Xtype ALLELES \
        -verbosity ERROR > {output}
        """