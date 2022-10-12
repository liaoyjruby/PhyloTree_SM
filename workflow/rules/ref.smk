# configfile: "config/config.yaml"

rule get_refs2:
    output:
        fasta_ref = "references/Homo_sapiens_assembly38.fasta",
        fasta_ref_idx = "references/Homo_sapiens_assembly38.fasta.fai",
        dict_ref = "references/Homo_sapiens_assembly38.dict",
        dbsnp_ref = "references/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbsnp_ref_idx = "references/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
    shell:
        """
        gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta references/
        gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai references/
        gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict references/
        gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf references/
        gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx references/
        """

# rule get_ref:
#     output:
#         "references/Homo_sapiens_assembly38.{type}"
#     shell:
#         """
#         gsutil -o \"GSUtil:parallel_process_count=1\" -m cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{wildcards.type} references/
#         """
