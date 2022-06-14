
# Usage 

## Quick Start

1. Install [Singularity](https://www.sylabs.io/guides/3.0/user-guide/) 
2. Create a dedicated conda environment 
```
conda create -n mosaicatcher_env -c conda-forge -c bioconda snakemake pandas pysam imagemagick && conda activate mosaicatcher_env
```
3. Clone the repository and go in the `workflow` directory
``` 
git clone https://git.embl.de/tweber/mosaicatcher-update.git && cd mosaicatcher-update/workflow/
```
4. Download test and reference data 
```
snakemake -c1 --config mode=download_data dl_external_files=True dl_bam_example=True input_bam_location=TEST_EXAMPLE_DATA/
```
5. Run on test data
```
snakemake --cores 12 input_bam_location=TEST_EXAMPLE_DATA/ output_location=TEST_OUTPUT/
```

6. Start running your own analysis
```
snakemake --cores 12 --config input_bam_location=<INPUT_DATA_FOLDER> output_location=<OUTPUT_DATA_FOLDER>

```
7. Generate report 
```
snakemake --cores 12 --config input_bam_location=<INPUT_DATA_FOLDER> output_location=<OUTPUT_DATA_FOLDER> --report <REPORT.zip>
```

## System requirements

This workflow is meant to be run in a Unix-based operating system (tested on Ubuntu 18.04 & CentOS 7). 

Minimum system requirements vary based on the use case. We highly recommend running it in a server environment with 32+GB RAM and 12+ cores.


- [Conda install instructions](https://conda.io/miniconda.html)
- [Singularity install instructions](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps)

## Detailed usage

### 🐍 1. Mosaicatcher basic conda environment install

MosaiCatcher leverages snakemake built-in features such as execution within container and conda predefined modular environments. That's why it is only necessary to create an environment that relies on [snakemake](https://github.com/snakemake/snakemake) (to execute the pipeline) and [pandas](https://github.com/pandas-dev/pandas) (to handle basic configuration). If you plan to generate HTML Web report including plots, it is also necessary to install [imagemagick](https://github.com/ImageMagick/ImageMagick). Finally, [pysam](https://pysam.readthedocs.io/en/latest/api.html) is currently required to enable `check_sm_tag` feature that compare BAM SM tag to folder name. 


If possible, it is also highly recommended to install and use `mamba` package manager instead of `conda`, which is much more efficient.

```
conda install -c conda-forge mamba
mamba create -n mosaicatcher_env -c conda-forge -c bioconda snakemake pandas pysam imagemagick
conda activate mosaicatcher_env
```



### ⤵️ 2. Clone repository & go into workflow directory

After cloning the repo, go into the `workflow` directory which correspond to the pipeline entry point. 




```
git clone https://git.embl.de/tweber/mosaicatcher-update.git
cd mosaicatcher-update/workflow/
```


### ⚙️ 3. MosaiCatcher config and execution

MosaiCatcher takes different arguments to run. Default configuration (`worfklow/config/config.yaml`) looks like the following. 

```yaml
# Required arguments

## Modes ["count", "segmentation", "mosaiclassifier"]
mode: "count"
## Plot enabled [True] or disabled [False]
plot: False
## Enable / Disable comparison for each BAM file between folder name & SM tag
check_sm_tag: True
## Enable / Disable download of BAM examples (RPE-BM510)
dl_bam_example: False
## Enable / Disable download of external files (1000G SNV & Fasta ref genome)
dl_external_files: False
## Input BAM location
input_bam_location: "TEST_EXAMPLE_DATA/"
## Output location
output_location: "TEST_OUTPUT/"

# External files
## 1000G SNV sites to genotype : https://sandbox.zenodo.org/record/1060653/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz
snv_sites_to_genotype: "sandbox.zenodo.org/record/1062182/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz"
# Reference genome : https://sandbox.zenodo.org/record/1060653/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
reference: "sandbox.zenodo.org/record/1062182/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Chromosomes list to process
chromosomes: [chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX]
```

You can either change it or override YAML file by using snakemake CLI arguments as the following : 

```
--config mode=segmentation plot=False input_bam_location=/HELLO_WORLD output_location=/AU_REVOIR
```

The `--config` argument will here overrides value of each of the keys present in the YAML file.




#### 3A. Download external data automatically with snakemake [Optional] 

```
snakemake -c1 --config mode=download_data dl_external_files=True
```

#### 3B. Strand-Seq BAM input data

##### (i) Download example data automatically with snakemake [Optional] 

```
snakemake -c1 --config mode=download_data dl_bam_example=True input_bam_location=TEST_EXAMPLE_DATA/
```
**Warning:** Download example data currently requires 35GB of free space disk. 


##### (ii) Use your own data

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

#### Local execution (without batch scheduler)

After defining your configuration, you can launch the pipeline the following way if you downloaded BAM example data using 3A:


```bash
snakemake \
    --cores 20 \
    --config \
        plot=True \
        mode=mosaiclassifier
```

Otherwise, you must specify your input and output folder like the following:

```bash
snakemake \
    --cores 20 \
    --config \
        plot=True \
        mode=mosaiclassifier \
        output_location=OUTPUT_FOLDER \
        input_bam_location=INPUT_FOLDER
```

---
**ℹ️ Note**

It is recommended to first run the command and check if there are any anomalies with the `--dry-run` option

---

---
**⚠️ Warning**

If you are experiencing any issues with conda-frontend snakemake option, please use `--conda-frontend conda` instead of `mamba` 

---

#### HPC execution

MosaiCatcher can be executed on HPC using [Slurm](https://slurm.schedmd.com/documentation.html) by leveraging snakemake profile feature. Current Slurm profile [`workflow/profiles/slurm/config.yaml`] was defined and tested on EMBL HPC cluster but can be modified, especially regarding **partition** setting. 

##### Current strategy to solve HPC job OOM 

Workflow HPC execution usually needs to deal with out of memory (OOM) errors, out of disk space, abnormal paths or missing parameters for the scheduler. To deal with OOM, we are currently using snakemake restart feature (thanks [@Pablo Moreno](https://github.com/pcm32)) in order to automatically double allocated memory to the job at each attempt (limited to 8 for the moment). Then, if a job fails to run with the default 1GB of memory allocated, it will be automatically restarted tith 2GB at the 2nd attempt, 4GB at the 3rd, etc. 

To execute MosaiCatcher on HPC, use the following command. 

##### Command 

```bash
snakemake \
    --profile profiles/slurm/ \
    --config \
        plot=True \
        mode=mosaiclassifier \
        output_location=OUTPUT_FOLDER \
        input_bam_location=INPUT_FOLDER
```

The `logs` and `errors` directory will be automatically created in the current directory, corresponding respectively to the `output` and `error` parameter of the `sbatch` command. 


###  📊 5. Generate report  [Optional]

Optionally, you can also MosaiCatcher rules that produce plots 

```bash
snakemake \
    --cores 20  \
    --config \
        plot=True \
        mode=mosaiclassifier \ 
        output_location=OUTPUT_FOLDER \
        input_bam_location=INPUT_FOLDER  \
    --report report.zip
```



---
**ℹ️ Note**

The zip file produced can be heavy (~1GB for 24 HGSVC samples ; 2000 cells) if multiple samples are processed in parallel in the same output folder.

---