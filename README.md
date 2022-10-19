# PhyloTree_SM

RNA-seq SNP based Phylogenetic Tree Pipeline, implemented as a Snakemake workflow

[Original Shell-based Pipeline Notes](https://docs.google.com/document/d/1BgOVv_zX04O1sf_FatdRbjEjrhTxYfw8i2QTu3F1fp4/edit#heading=h.epnsuaukicdh)

## Usage

### 1. Install Snakemake
Using either the [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) (recommended) or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager, install Snakemake & Snakedeploy in an isolated environment:
```shell
mamba create -c conda-forge -c bioconda -n snakemake snakemake snakedeploy
```
Ensure that the newly created environment is activated for all following steps:
```shell
conda activate snakemake
```

### 2. Deploy workflow

Create a working directory for this project and enter it for all following steps:
```shell
mkdir -p path/to/workdir
cd path/to/workdir
```
Then run:
```
snakedeploy deploy-workflow https://github.com/liaoyjruby/PhyloTree_SM . --branch main
```

Snakedeploy creates two folders:
- `workflow`: contains the Snakemake rules that implement the workflow
- `config`: contains configuration files that should be edited according to needs

### 3. Configure workflow

General Settings:

Modify `config.yaml` as needed according to comments in the file.

Units & Samples Sheets:

- `units.tsv`: Required columns `Sample_ID` and `ID`. Add column `Mapped_Path` with absolute paths if aligned BAM files are elsewhere.
- `samples.tsv`: Sample annotation sheet with required column `Sample_ID`. Add columns with information about conditions of interest as desired.

The pipeline will include all samples listed in the `units.tsv` sheet in the final phylogenetic tree output.'

### 5. Run workflow

After configuration, run the Snakemake workflow while deploying any necessary software in the process with:
```shell
snakemake --cores all --use-conda
```
The main driver script `Snakefile` in the `workflow` subfolder will automatically be detected and executed.