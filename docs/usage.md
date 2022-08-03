
# Usage 

# Quick Start

0. [Optional] Install [Singularity](https://www.sylabs.io/guides/3.0/user-guide/) 

1. Create a dedicated conda environment 
```bash
conda create -n snakemake -c bioconda snakemake && conda activate snakemake
```

2. Clone the repository 
```bash
git clone https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git && cd mosaicatcher-pipeline
```

3. Download reference data 
```bash
snakemake -c1 --config dl_external_files=True 
```

4. Run on example data on only one small chromosome (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)
```bash
# Snakemake Profile: if singularity installed: workflow/profiles/local/conda_singularity/
# Snakemake Profile: if singularity NOT installed: workflow/profiles/local/conda/
snakemake --cores 6 --configfile .tests/config/simple_config.yaml --profile workflow/profiles/local/conda_singularity/ 
```

5. Generate report on example data
```bash
snakemake --cores 6 --configfile .tests/config/simple_config.yaml --profile workflow/profiles/local/conda_singularity/ --report report.zip
```

6. Start running your own analysis
```bash
snakemake \
    --cores <N> --config input_bam_location=<INPUT_DATA_FOLDER> output_location=<OUTPUT_DATA_FOLDER> \
    --profile workflow/profiles/local/conda_singularity/ 

```
7. Generate report 
```bash
snakemake \
    --cores <N> --config input_bam_location=<INPUT_DATA_FOLDER> output_location=<OUTPUT_DATA_FOLDER> \
    --profile workflow/profiles/local/conda_singularity/ \
    --report <REPORT.zip>
```

## System requirements

This workflow is meant to be run in a Unix-based operating system (tested on Ubuntu 18.04 & CentOS 7). 

Minimum system requirements vary based on the use case. We highly recommend running it in a server environment with 32+GB RAM and 12+ cores.


- [Conda install instructions](https://conda.io/miniconda.html)
- [Singularity install instructions](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps)

## Detailed usage

### 🐍 1. Mosaicatcher basic conda environment install

MosaiCatcher leverages snakemake built-in features such as execution within container and conda predefined modular environments. That's why it is only necessary to create an environment that relies on [snakemake](https://github.com/snakemake/snakemake) (to execute the pipeline) and [pandas](https://github.com/pandas-dev/pandas) (to handle basic configuration). If you plan to generate HTML Web report including plots, it is also necessary to install [imagemagick](https://github.com/ImageMagick/ImageMagick). Finally, [pysam](https://pysam.readthedocs.io/en/latest/api.html) and [tqdm](https://github.com/tqdm/tqdm) are currently required to enable `check_sm_tag` feature that compare BAM SM tag to folder name. 


If possible, it is also highly recommended to install and use `mamba` package manager instead of `conda`, which is much more efficient.

```
conda install -c conda-forge mamba
mamba create -n snakemake -c bioconda snakemake
conda activate mosaicatcher_env
```



### ⤵️ 2. Clone repository & go into workflow directory

After cloning the repo, go into the `workflow` directory which correspond to the pipeline entry point. 




```
git clone https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git
cd mosaicatcher-pipeline/workflow/
```


### ⚙️ 3. MosaiCatcher execution


#### 3A. Download external data automatically with snakemake [Optional] 

```
snakemake -c1 --config dl_external_files=True
```

#### 3B. Strand-Seq BAM input data

##### (i) Download large example data automatically with snakemake [Optional] 

```
snakemake -c1 --config dl_bam_example=True input_bam_location=TEST_EXAMPLE_DATA/
```
**Warning:** Download example data currently requires 3GB of free space disk. 


##### (ii) Use your own data

In its current flavour, MosaiCatcher requires that input data must be formatted the following way :


```bash
Parent_folder
|-- Sample_1
|   `-- all
|       |-- Cell_01.sort.mdup.bam
|       |-- Cell_02.sort.mdup.bam
|       |-- Cell_03.sort.mdup.bam
|       `-- Cell_04.sort.mdup.bam
| 
`-- Sample_2
    `-- all
        |-- Cell_21.sort.mdup.bam
        |-- Cell_22.sort.mdup.bam
        |-- Cell_23.sort.mdup.bam
        `-- Cell_24.sort.mdup.bam
```

In a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Your Strand-seq BAM files of this sample go into the following folder:

* `all` for the total set of BAM files

It is important to follow these rules for single-cell data

* One BAM file per cell
* Sorted and indexed
  * If BAM files are not indexed, please use a writable folder in order that the pipeline generate itself the index `.bai` files
* BAM file name ending by suffix: `.sort.mdup.bam`
* Timestamp of index files must be newer than of the BAM files
* Each BAM file must contain a read group (`@RG`) with a common sample name (`SM`), which must match the folder name (`sampleName` above)



### ⚡️ 4. Run the pipeline

#### Local execution (without batch scheduler)

After defining your configuration, you can launch the pipeline the following way if you downloaded BAM example data using 3A:


```bash
snakemake \
    --cores <N> \
    --profile workflow/profiles/local/conda_singularity/
 
```

Otherwise, you must specify your input and output folder like the following:

```bash
snakemake \
    --cores <N> \
    --config \
        input_bam_location=<INPUT_FOLDER> \
        output_location=<OUTPUT_FOLDER> \
    --profile workflow/profiles/local/conda_singularity/ --singularity-args "-B /<mounting_point>:/<mounting_point>"
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
    --config \
        input_bam_location=<INPUT_FOLDER> \
        output_location=<OUTPUT_FOLDER> \
    --singularity-args "-B /<mounting_point>:/<mounting_point>" \
    --profile workflow/profiles/slurm/
```

The `logs` and `errors` directory will be automatically created in the current directory, corresponding respectively to the `output` and `error` parameter of the `sbatch` command. 


###  📊 5. Generate report  [Optional]

Optionally, you can also MosaiCatcher rules that produce plots 

```bash
snakemake \
    --config \
        input_bam_location=<INPUT_FOLDER> \
        output_location=<OUTPUT_FOLDER> \
    --singularity-args "-B /<mounting_point>:/<mounting_point>" \
    --profile workflow/profiles/slurm/ \
    --report <report>.zip
```



---
**ℹ️ Note**

The zip file produced can be heavy (~1GB for 24 HGSVC samples ; 2000 cells) if multiple samples are processed in parallel in the same output folder.

---