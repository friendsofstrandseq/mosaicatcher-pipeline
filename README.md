![MosaiCatcher](docs/mosaic_logo.png)


Structural variant calling from single-cell Strand-seq data - summarized in a single [Snakemake](https://github.com/snakemake/snakemake) pipeline.


##  Overview of this workflow

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) to
execute all steps of MosaiCatcher in order. The starting point are single-cell
BAM files from Strand-seq experiments and the final output are SV predictions in
a tabular format as well as in a graphical representation. To get to this point,
the workflow goes through the following steps:

  1. Binning of sequencing reads in genomic windows of 100kb via [mosaic](https://github.com/friendsofstrandseq/mosaicatcher)
  2. Strand state detection
  3. [Optional]Normalization of coverage with respect to a reference sample
  4. Multi-variate segmentation of cells ([mosaic](https://github.com/friendsofstrandseq/mosaicatcher))
  5. Haplotype resolution via [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR)
  6. Bayesian classification of segmentation to find SVs using MosaiClassifier
  7. Visualization of results using custom R plots



## System requirements

This workflow is meant to be run in a Unix-based operating system (tested on Ubuntu 18.04 & CentOS 7). 

Minimum system requirements vary based on the use case. We highly recommend running it in a server environment with 32+GB RAM and 24 cores.


- [Conda install instructions](https://conda.io/miniconda.html)
- [Singularity install instructions](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps)

## Installation

### 🐍 1. Mosaicatcher basic conda environment install

MosaiCatcher leverages snakemake built-in features such as execution within container and conda predefined modular environments. That's why it is only necessary to create an environment that relies on [snakemake](https://github.com/snakemake/snakemake) (to execute the pipeline) and [pandas](https://github.com/pandas-dev/pandas) (to handle basic configuration). If you plan to generate HTML Web report including plots, it is also necessary to install [imagemagick](https://github.com/ImageMagick/ImageMagick).

If possible, it is also highly recommanded to install and use mamba package manager instead of conda, which is much more efficient.

```
conda install -c conda-forge mamba
mamba create -n mosaicatcher_env -c conda-forge -c bioconda snakemake pandas imagemagick
conda activate mosaicatcher_env
```



### ⤵️ 2. Clone repository & go into workflow directory

After cloning the repo, go into the `workflow` directory which correspond to the pipeline entry point. 

<<<<<<< HEAD



```
GIT_LFS_SKIP_SMUDGE=1 git clone https://git.embl.de/tweber/mosaicatcher-update.git
cd mosaicatcher-update/workflow/
```

---

**ℹ️ Note**

`GIT_LFS_SKIP_SMUDGE=1` is currently added to `git clone` command to prevent cloning of the LFS files tracked by git (35GB of data).
*(Will be corrected in next release)*

---


=======
```
git clone https://git.embl.de/tweber/mosaicatcher-update.git
cd mosaicatcher-update/workflow/
```

>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
### ⚙️ 3. MosaiCatcher config and execution

MosaiCatcher takes different arguments to run. Default configuration (`worfklow/config/config.yaml`) looks like the following. 

```yaml
# Required arguments

## Modes ["count", "segmentation", "mosaiclassifier"]
mode: "count"
## Plot enabled [True] or disabled [False]
plot: False
## Enable / Disable comparison for each BAM file between folder name & SM tag
check_sm_tag: False

# Chromosomes list to process
chromosomes: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX]


# I/O path

## Input BAM location
input_bam_location: "TEST_EXAMPLE_DATA/bam/"
## Output location
output_location: "TEST_OUTPUT"


# External files 

## 1000G SNV sites to genotype : https://sandbox.zenodo.org/record/1060653/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz
snv_sites_to_genotype: "/path/to/SNV_sites"

## Reference genome : https://sandbox.zenodo.org/record/1060653/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
reference: "/path/to/ref"
```
<<<<<<< HEAD

You can either change it or override YAML file by using snakemake CLI arguments as the following : 

```
--config mode=segmentation plot=False input_bam_location=/HELLO_WORLD output_location=/AU_REVOIR
```

The `--config` argument will here overrides value of each of the keys present in the YAML file.

Following commands allow to retrieve 1000G VCF file (+ .tbi index) as well as Fasta reference genome file (+ .fai index).

```
wget https://sandbox.zenodo.org/record/1060653/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz
wget https://sandbox.zenodo.org/record/1060653/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz.tbi
wget https://sandbox.zenodo.org/record/1060653/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
wget https://sandbox.zenodo.org/record/1060653/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
```
=======
You can either change it or override YAML file by using snakemake CLI arguments as the following : 
```
--config mode=segmentation plot=False input_bam_location=/HELLO_WORLD output_location=/AU_REVOIR
```
The `--config` argument will here overrides value of each of the keys present in the YAML file.

>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4


#### 3A. Download example data automatically with snakemake [Optional] 

```
snakemake -c1 --config mode=download_example_data input_bam_location=/path/to/INPUT
```
<<<<<<< HEAD
**Warning:** Download example data currently requires 35GB of free space disk. 
=======
**Warning:** Download example data currently requires 65GB of free space disk. 
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4


#### 3B. Prepare input data 

In its current flavour, MosaiCatcher requires that input data must be formatted the following way :


```bash
Parent_folder
|-- Sample_1
|   |-- all
|   |   |-- Cell_01.sort.mdup.bam
|   |   |-- Cell_02.sort.mdup.bam
|   |   |-- Cell_03.sort.mdup.bam
|   |   `-- Cell_04.sort.mdup.bam
|   `-- selected
|       |-- Cell_01.sort.mdup.bam
|       |-- Cell_02.sort.mdup.bam
|       `-- Cell_04.sort.mdup.bam
`-- Sample_2
    |-- all
    |   |-- Cell_21.sort.mdup.bam
    |   |-- Cell_22.sort.mdup.bam
    |   |-- Cell_23.sort.mdup.bam
    |   `-- Cell_24.sort.mdup.bam
    `-- selected
        |-- Cell_22.sort.mdup.bam
        |-- Cell_23.sort.mdup.bam
        `-- Cell_24.sort.mdup.bam
```

In a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Your Strand-seq BAM files of this sample go into two folders:

* `all/` for the total set of BAM files
* `selected/` for the subset of successful Strand-seq libraries (possibly hard-linked to `all/`)

It is important to follow these rules for single-cell data

* One BAM file per cell
* Sorted and indexed
  * If BAM files are not indexed, please use a writable folder in order that the pipeline generate itself the index `.bai` files
* Timestamp of index files must be newer than of the BAM files
* Each BAM file must contain a read group (`@RG`) with a common sample name (`SM`), which must match the folder name (`sampleName` above)


### ⚡️ 4. Run the pipeline

After defining your configuration, you can launch the pipeline the following way:


<<<<<<< HEAD



```bash
snakemake \
    --use-conda  \
    --cores 20  \
=======
```bash
snakemake \
    --use-conda  \
    --cores 40  \
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
    --config \
        plot=True \
        mode=mosaiclassifier \ 
        output_location=/path/to/OUTPUT_FOLDER/ \
        input_bam_location=/path/to_/INPUT_FOLDER/  \
<<<<<<< HEAD
    -p \
    --conda-frontend mamba \
    --use-singularity \
    --singularity-args "-B /:/" \
    --dry-run
```

---
**⚠️ Warning**

If you are experiencing any issues with conda-frontend snakemake option, please use `--conda-frontend conda` instead of `mamba` 

---


=======
    --printshellcmds \
    --conda-frontend mamba \
    --use-singularity \
    --singularity-args "-B /:/"
```

>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4

#### Snakemake & Singularity arguments

```
--cores 1
```
Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the number of total cores used over all jobs (made available to rules via workflow.cores).

```
--printshellcmds, -p
```
Recommended to print out the shell commands that will be executed.

```
--use-conda
```
If defined in the rule, run job in a conda environment. If this flag is not set, the conda directive is ignored and use the current environment (and path system) to execute the command.

```
--conda-frontend mamba|conda 
```
Choose the conda frontend for installing environments. Mamba is much faster and highly recommended. Default: “mamba”



```
--use-singularity 
```
If defined in the rule, run job within a singularity container. If this flag is not set, the singularity directive is ignored.

```
--singularity-args "-B /:/"
```
Pass additional args to singularity. `-B` stands for binding point between the host and the container.

---
**ℹ️ Note**

Currently, raise the following WARNING message : 
```
WARNING: Skipping mount /etc/localtime [binds]: /etc/localtime doesn't exist in container
```
---

Obviously, all other snakemake CLI options can also be used. 



#### MosaiCatcher arguments

!# TODO : simplified pipeline image  

MosaiCatcher currently supports three different modes of execution : `count`, `segmentation` and `mosaiclassifier`.
- `count` (selected by default) will only performs `Mosaic count` binning and count reads for each bin produced
- `segmentation` will run the pipeline until the `Mosaic segmentation` and selection of the correct segments
- `mosaiclassifier` will run the complete pipeline until the detection of SV in each selected cell of the samples

To select your mode of execution, use the following argument `--config mode=[count|segmentation|mosaiclassifier]`

For each of these modes, you can *enable* or *disable* the plots generation by using `--config plot=[True|False]`


###  📊 5. Generate report  Optional]

Optionally, you can also MosaiCatcher rules that produce plots 

```bash
snakemake \
    --use-conda  \
<<<<<<< HEAD
    --cores 20  \
=======
    --cores 40  \
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
    --config \
        plot=True \
        mode=mosaiclassifier \ 
        output_location=OUTPUT_FOLDER \
        input_bam_location=INPUT_FOLDER  \
    --report report.zip
```

<<<<<<< HEAD
## 📆 Roadmap 

- [ ] Zenodo automatic download of external files + indexes
- [ ] Change of reference genome (currently only GRCh38)
- [ ] Plotting options (enable/disable segmentation back colors)
- [ ] Full singularity image with preinstalled conda envs
- [ ] Automatic testing of BAM SM tag compared to sample folder name
- [ ] On-error/success e-mail
- [ ] Upstream QC pipeline and FastQ handle  

=======
>>>>>>> 726eb7567c393d423926ee719d112336c279d4a4
## 📕 References


> Strand-seq publication: Falconer, E., Hills, M., Naumann, U. et al. DNA template strand sequencing of single-cells maps genomic rearrangements at high resolution. Nat Methods 9, 1107–1112 (2012). https://doi.org/10.1038/nmeth.2206

> scTRIP/MosaiCatcher original publication: Sanders, A.D., Meiers, S., Ghareghani, M. et al. Single-cell analysis of structural variations and complex rearrangements with tri-channel processing. Nat Biotechnol 38, 343–354 (2020). https://doi.org/10.1038/s41587-019-0366-x


