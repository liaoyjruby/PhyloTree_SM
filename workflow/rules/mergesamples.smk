import pandas as pd

samples = list(pd.read_table(config["units"])['Sample_ID'])

rule build_VCF_list:
    input:
        expand("hcVCF/{sample}.vcf.gz", sample=samples)
    output:
        "merged/VCFList.txt"
    run:
        print({input})
        with open("merged/VCFList.txt", "w") as f:
            f.write('\n'.join(input))

rule merge_VCFs:
    input:
        vcflist = "merged/VCFList.txt",
        vcfs = expand("hcVCF/{sample}.vcf.gz", sample=samples),
        vcfidxs = expand("hcVCF/{sample}.vcf.gz.tbi", sample=samples)
    output:
        vcf="merged/VCFMerged.vcf.gz",
        vcfidx="merged/VCFMerged.vcf.gz.tbi"
    log:
        "logs/merged/VCFMerged.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        'bcftools merge '
            '-l {input.vcflist} '
            '-o {output.vcf} '
            '-O z &> {log} ; '
        'gatk IndexFeatureFile -I {output.vcf} &> {log}'

rule get_vcf2phylip_script:
    output:
        "workflow/rules/scripts/vcf2phylip.py"
    log:
        "logs/scripts/get_vcf2phylip_script.log"
    conda:
        "../envs/wget.yaml"
    shell:
        'wget -O {output} https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py  &> {log}'

rule vcf2phylip:
    input:
        vcf2phy=rules.get_vcf2phylip_script.output,
        vcf="merged/VCFMerged.vcf.gz",
        vcfidx="merged/VCFMerged.vcf.gz.tbi"
    output:
        "merged/VCFMerged.min4.phy"
    log:
        "logs/tree/vcf2phylip.log"
    conda:
        "../envs/cmd.yaml"
    shell:
        'python {input.vcf2phy} '
            '-i {input.vcf} '
            '--output-prefix merged/VCFMerged '
            '-m 4 &> {log}'

rule iqtree:
    input:
        "merged/VCFMerged.min4.phy"
    output:
        "merged/tree/{treename}.treefile"
    log:
        "logs/tree/iqtree.log"
    conda:
        "../envs/iqtree.yaml"
    shell:
        'iqtree2 -s {input} '
            '-pre merged/tree/{wildcards.treename} '
            '-T AUTO '
            '-st DNA '
            '-redo '
            '-m GTR &> {log}'

rule plot_tree:
    input:
        "merged/tree/VCFMerged.treefile"
    output:
        "figures/tree_rect.png",
        "figures/tree_circ.png",
        "figures/tree_eqan.png",
        "figures/tree_daylight.png",
        "figures/heatmap_distances.png"
    conda:
        "../envs/r.yaml"
    script:
        "scripts/plot_tree.R"