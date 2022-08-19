# Usage

## Quick Start

0. [Optional] Install [Singularity](https://www.sylabs.io/guides/3.0/user-guide/)

1. Create a dedicated conda environment

```bash
conda create -n snakemake -c bioconda "snakemake=7.4.1" && conda activate snakemake
```

2. Clone the repository

```bash
git clone https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git && cd mosaicatcher-pipeline
```

3. Run on example data on only one small chromosome (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)

```bash
# Snakemake Profile: if singularity installed: workflow/profiles/local/conda_singularity/
# Snakemake Profile: if singularity NOT installed: workflow/profiles/local/conda/
snakemake --cores 6 --configfile .tests/config/simple_config.yaml --profile workflow/profiles/local/conda_singularity/
```

4. Generate report on example data

```bash
snakemake --cores 6 --configfile .tests/config/simple_config.yaml --profile workflow/profiles/local/conda_singularity/ --report report.zip
```

## 🔬​ Start running your own analysis

Following commands show you an example using local execution (not HPC or cloud)

1. Chose your reference (hg38 or hg19 ; T2T is coming soon) and download related data

```bash
snakemake -c1 --config dl_external_files=True reference=hg38
```

2. Start running your own Strand-Seq analysis

```bash
snakemake \
    --cores <N> --config input_bam_location=<INPUT_DATA_FOLDER> output_location=<OUTPUT_DATA_FOLDER> \
    --profile workflow/profiles/local/conda_singularity/

```

3. Generate report

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

```bash
conda install -c conda-forge mamba
mamba create -n snakemake -c bioconda "snakemake=7.4.1"
conda activate mosaicatcher_env
```

### ⤵️ 2. Clone repository & go into workflow directory

After cloning the repo, go into the `workflow` directory which correspond to the pipeline entry point.

```bash
git clone https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git
cd mosaicatcher-pipeline/workflow/
```

### ⚙️ 3. MosaiCatcher execution (without preprocessing)

#### 3A. Download external data automatically with snakemake [Optional]

```bash
snakemake -c1 --config dl_external_files=True reference=[hg38|hg19] # T2T coming in the next release, mm10 soon
```

#### 3B. Strand-Seq BAM input data

#### (i) Download & use example data automatically with snakemake [Optional]

In a first step, `dl_bam_example=True` will allow to retrieve data stored on Zenodo registry.

```bash
snakemake -c1 --config dl_bam_example=True input_bam_location=TEST_EXAMPLE_DATA/
```

Then, the following command will process the 18 cells (full BAM with all chromosomes) present in this example.

```bash
snakemake -c6 --config input_bam_location=TEST_EXAMPLE_DATA output_location=TEST_OUTPUT \
    --profile workflow/profiles/local/conda_singularity/
```

**Warning:** Download example data currently requires 3GB of free space disk.

#### (ii) Use your own data

#### BAM input requirements

It is important to follow these rules for Strand-Seq single-cell BAM data

- **BAM file name ending by suffix: `.sort.mdup.bam`**
- One BAM file per cell
- Sorted and indexed
  - If BAM files are not indexed, please use a writable folder in order that the pipeline generate itself the index `.bai` files
- Timestamp of index files must be newer than of the BAM files
- Each BAM file must contain a read group (`@RG`) with a common sample name (`SM`), which must match the folder name (`sampleName` above)

#### Note: filtering of low-quality cells impossible to process

From Mosaicatcher version ≥ 1.6.1, a snakemake checkpoint rule was defined ([filter_bad_cells](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/blob/master/workflow/rules/count.smk#L91)) to automatically remove from the analysis cells flagged as low-quality.

The ground state to define this status (usable/not usable) is determined from column pass1 (enough coverage to process the cell) in mosaic count info file (see example below). The value of pass1 is not static and can change according the window size used (the larger the window, the higher the number of reads to retrieve in a given bin).

Cells flagged as low-quality cells are listed in the following TSV file: `<OUTPUT>/cell_selection/<SAMPLE>/labels.tsv`.

#### Classic behavior

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

- `all` for the total set of BAM files

> Using the classic behavior, cells flagged as low-quality will only be determined based on coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process).

#### Old behavior

Previous version of MosaiCatcher (version ≤ 1.5) needed not only a `all` directory as described above, but also a `selected` folder, presenting only high-quality selected libraries wished to be processed for the rest of the analysis.

You can still use this behavior by enabling the config parameter either by the command line: `input_old_behavior=True` or by modifying the corresponding entry in the config/config.yaml file.

Thus, in a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Your Strand-seq BAM files of this sample go into the following folder:

- `all` for the total set of BAM files
- `selected` for the subset of BAM files that were identified as high-quality (either by copy or symlink)

Your `<INPUT>` directory should look like this:

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

> Using the `old behavior`, cells flagged as low-quality will be determined both based on their presence in the `selected` folder presented above and on coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process).

---

**⚠️ Warning**

Using the `old behavior`, only **intersection** between cells present in the selected folder and with enough coverage will be kept. Example: if a library is present in the selected folder but present a low coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process), this will not be processed.

---

#### 3C. FASTQ input & Preprocessing module\*

From Mosaicatcher version ≥ 1.6.1, it is possible to use [ashleys-qc-pipeline preprocessing module](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) as part of MosaiCatcher. To do so, the user need to enable the config parameter `ashleys_pipeline=True` and create a directory respecting the structure below (based on the model used for BAM inputs):

```bash
Parent_folder
|-- Sample_1
|   `-- all
|       |-- Cell_01.1.fastq.gz
|       |-- Cell_01.2.fastq.gz
|       |-- Cell_02.1.fastq.gz
|       |-- Cell_02.2.fastq.gz
|       |-- Cell_03.1.fastq.gz
|       |-- Cell_03.2.fastq.gz
|       |-- Cell_04.1.fastq.gz
|       `-- Cell_04.2.fastq.gz
|
`-- Sample_2
    `-- all
        |-- Cell_21.1.fastq.gz
        |-- Cell_21.2.fastq.gz
        |-- Cell_22.1.fastq.gz
        |-- Cell_22.2.fastq.gz
        |-- Cell_23.1.fastq.gz
        |-- Cell_23.2.fastq.gz
        |-- Cell_24.1.fastq.gz
        `-- Cell_24.2.fastq.gz
```

Thus, in a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Each Strand-seq FASTQ files of this sample need to go into the `fastq` folder and respect the following syntax: `<CELL>.<1|2>.fastq.gz`, `1|2` corresponding to the pair identifier.

Informations and modes of execution can be found on the [ashleys-qc-pipeline documentation](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/blob/main/README.md).

---

**⚠️ Warnings**

1. When using ashleyq-qc-pipeline preprocessing module, Singularity execution is not possible at the moment (will be fixed soon).
2. `ashleys_pipeline=True` and `input_old_behavior=True` are mutually exclusive

---

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

### 📊 5. Generate report [Optional]

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
