# PhyloTree_SM

RNA-seq SNP based Phylogenetic Tree Pipeline, implemented as a Snakemake workflow

## Links

[09/14 Graeber Lab Presentation](https://docs.google.com/presentation/d/13BMLpEajmpYnvCELF9d1-8kaRoVbmbUg/edit?usp=sharing&ouid=117414729187227390562&rtpof=true&sd=true)

[Original Non-Snakemake Pipeline Notes](https://docs.google.com/document/d/1BgOVv_zX04O1sf_FatdRbjEjrhTxYfw8i2QTu3F1fp4/edit#heading=h.epnsuaukicdh)

[Original Pipeline Implementation using Shell](https://github.com/liaoyjruby/glab_R/tree/main/OV/PhyloTree)

## Usage

### 1. Install Snakemake
Using either the [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) (recommended) or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager, install Snakemake & Snakedeploy in an isolated environment:
```shell
mamba create -c conda-forge -c bioconda -n snakemake snakemake snakedeploy
```
Ensure that the newly created environment is activated for all following steps:
```shell
mamba activate snakemake
```

### 2. Deploy workflow

Create a working directory for this project and enter it for all following steps:
```shell
mkdir -p path/to/workdir
cd path/to/workdir
```

If you want to run the pipeline according to the main branch of this repository, run:
```shell
snakedeploy deploy-workflow https://github.com/liaoyjruby/PhyloTree_SM . --branch main
```

If you want to have all files locally, clone this repository into the working directory:
```shell
git clone https://github.com/liaoyjruby/PhyloTree_SM.git .
```

There are two main folders:
- `workflow`: contains the Snakemake rule that implement the workflow
- `config`: contains configuration files that should be edited according to needs



### 3. Configure workflow

**General Settings:**

Modify `config.yaml` as needed according to comments in the file.

**Units & Samples Sheets:**

- `units.tsv`: Required columns `Sample_ID` and `ID`. Add column `Mapped_Path` with absolute paths if aligned BAM files are elsewhere.
- `samples.tsv`: Sample annotation sheet with required column `Sample_ID`. Add columns with information about conditions of interest as desired.

The pipeline will include all samples listed in the `units.tsv` sheet in the final phylogenetic tree output.

If you have the aligned BAM files already and do not want to copy them into the working directory, place them into subdirectory `mapped/` with name `<Sample_ID>.bam`.

If you have the VCF file + index, place them into subdirectory `hcVCF/` with name `<Sample_ID>.vcf.gz` & `<Sample_ID>.vcf.gz.tbi`.

### 5. Run workflow

See DAG of pipeline jobs by running:
```shell
snakemake -c 1 dag --use-conda
```

After configuration, run the Snakemake workflow while deploying any necessary software in the process with:
```shell
snakemake -c all --use-conda
```
The main script `Snakefile` in the `workflow` subfolder will automatically be detected and executed.

Change `all` to desired number of cores to use to run the pipeline.