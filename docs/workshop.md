# Mosaicatcher workshop

## AFAC

Context: You are working on a fancy project and Jan is suggesting at some point to generate some Strand-seq data

That's the first time you are working with Strand-seq and you are starting to panic

But You remember that you heard that a complex tool was developed in the lab in order to process in a systematic way Strand-seq data: tadam MOSAICATCHER

Prerequisite: I asked you to select a sample to process during today's workshop
TD: feedback on how they traced back the name of the sample, the associated run/flowcell, the date when it was sequenced ...

So here's the plan for today:

Small intro (~20/30 min) about Mosaicatcher, the different steps, options, branches, possibilities

Outputs examples
SV trustfullness

Web report analysis of RPE-MIXTURE

Hands on: pipeline install, module load, test data execution

    Data preparation, parent folder structure


    ashleys_pipeline_only=True

    vim cell_selection/labels.tsv

    ashleys_pipeline_only=False

    mkdir -p scNOVA_input_user

    # TODO: script to precreate input_subclonality.txt




    vim scNOVA_input_user/input_subclonality.txt

Then trigger the pipeline on YOUR data

Once this is running, web report analysis together with questions

Then, Strand-scape
Still in beta, some microservices instable, main application for QC and web report consultation
Remove MC trigger, too complex in the backend
Cell selection with username

    cp --preserve=timestamps FROM_ TO_

    snakemake ...

## Technical prerequisites

- SSHFS/SFTP connection to visualise/download/access files created (WinSCP/FileZilla/Cyberduck)
- Functional terminal connected to the EMBL cluster (if not follow SSH key configuration here: https://www.embl.org/internal-information/it-services/hpc-resources-heidelberg/)
- Have a workspace on /g/korbel

## Workshop prerequisites

- Pick a sample name to be processed
- Download this MosaiCatcher report: https://oc.embl.de/index.php/s/WBgrzBjyzdYdVJA/download

## EMBL cheatsheet

### connect to seneca

ssh USERNAME@seneca.embl.de

### connect to login nodes

ssh USERNAME@login0[1,2,3,4].embl.de (login01 to login04)

**ℹ️ Important Note**

From 2.2.0, you don't need to clone both [ashleys-qc-pipeline preprocessing module](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) and [mosaicatcher-pipeline](https://github.com/friendsofstrandseq/mosaicatcher-pipeline). By using `ashleys_pipeline_only=True` combined with `ashleys_pipeline=True` in the configuration of MosaiCatcher, this will stop the execution after the generation of the files related to ashleys-qc-pipeline. This allow you to use a single pipeline, repository and container and limit the potential conflicts by processing the same folder (data_location) by different repositories and set of files (including the workflow/data/ref_genomes references files).

---

## Snakemake cheat sheet & important things

Snakemake is a workflow manager, will handle execution and distribute computing based on configuration (local execution: use local CPUs, cluster execution, will generate the sbatch commands for each job, 1 single instance locally = snakemake, other jobs computed on the cluster, except when specifically defined)

To run the pipeline, go into the repository and run: snakemake

Snakemake important arguments/options

--dry-run
--profile
--config
--rerun-triggers
--touch

## MosaiCatcher important files

- Counts: PARENT_FOLDER/SAMPLE_NAME/counts/SAMPLE_NAME.txt.raw.gz
- Counts statistics: PARENT_FOLDER/SAMPLE_NAME/counts/SAMPLE_NAME.info_raw
- Ashleys predictions: PARENT_FOLDER/SAMPLE_NAME/cell_selection/labels.tsv
- Counts plot: PARENT_FOLDER/SAMPLE_NAME/plots/CountComplete.raw.pdf
- Count normalied plot: PARENT_FOLDER/SAMPLE_NAME/plots/CountComplete.normalised.pdf
- Phased W/C regions: PARENT_FOLDER/SAMPLE_NAME/strandphaser/strandphaser_phased_haps_merged.txt
- SV calls (stringent): PARENT_FOLDER/SAMPLE_NAME/mosaiclassifier/sv_calls/stringent_filterTRUE.tsv
- SV calls (lenient): PARENT_FOLDER/SAMPLE_NAME/mosaiclassifier/sv_calls/lenient_filterFALSE.tsv
- Plots folder: PARENT_FOLDER/SAMPLE_NAME/plots/
- scNOVA outputs:

## CLI usage of the pipeline

### Quick Start

1. Clone the repository

```bash
git clone --recurse-submodules https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git && cd mosaicatcher-pipeline
```

Notes

- MC single repository, ashleys do not need to be downloaded anymore, automatically handled
- Config definition is crucial / via command line or via YAML file, will define where to stop, which mode, which branch, which options to be used
- Profile

2. Load snakemake

A. Use module load OR create a dedicated conda environment

```bash
module load snakemake/7.32.4-foss-2022b
```

Give a look at the folder structure:

```bash
tree -h .tests/data_CHR17
```

Similar to this

Parent_folder
|-- Sample_1
| `-- fastq
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
    `-- fastq
|-- Cell_21.1.fastq.gz
|-- Cell_21.2.fastq.gz
|-- Cell_22.1.fastq.gz
|-- Cell_22.2.fastq.gz
|-- Cell_23.1.fastq.gz
|-- Cell_23.2.fastq.gz
|-- Cell_24.1.fastq.gz
`-- Cell_24.2.fastq.gz

````



1. Run on example data on only one small chromosome (`<disk>` must be replaced by your disk letter/name)

First using the `--dry-run` option of snakemake to make sure the Graph of Execution is properly connected. (In combination with `--dry-run`, we use the `local/conda` profile as snakemake still present a bug by looking for the singularity container).



```bash
snakemake \
    --cores 6 \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \ # DATA LOCATION
        ashleys_pipeline=True \ # DOWNLOAD & TRIGGER ASHLEYS QC UPSTREAM MODULE
        ashleys_pipeline_only=True \ # STOP AFTER ASHLEYS QC - VALIDATION PURPOSE
        multistep_normalisation=True \ # TRIGGER MARCO'S MULTISTEP NORMALISATION
        MultiQC=True \ # TRIGGER samtools stats, FastQC & MultiQC reporting
    --profile workflow/snakemake_profiles/local/conda/ \ # EXECUTION PROFILE TO BE USED
    --dry-run # ONLY CHECK IF EVERYTHING CONNECTS WELL AND READY FOR COMPUTING
````

If no error message, you are good to go!

```bash
snakemake \
    --cores 6 \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \ # DATA LOCATION
        ashleys_pipeline=True \ # DOWNLOAD & TRIGGER ASHLEYS QC UPSTREAM MODULE
        ashleys_pipeline_only=True \ # STOP AFTER ASHLEYS QC - VALIDATION PURPOSE
        multistep_normalisation=True \ # TRIGGER MARCO'S MULTISTEP NORMALISATION
        MultiQC=True \ # TRIGGER samtools stats, FastQC & MultiQC reporting
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --cores 24
```

```bash
cat .tests/data_CHR17/RPE-BM510/counts/RPE-BM510.info_raw
zcat .tests/data_CHR17/RPE-BM510/counts/RPE-BM510.txt.raw.gz | less
cat .tests/data_CHR17/RPE-BM510/cell_selection/labels.tsv
```

Look at the plots

.tests/data_CHR17/RPE-BM510/plots

REPORT

```bash
snakemake \
    --cores 6 \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \ # DATA LOCATION
        ashleys_pipeline=True \ # DOWNLOAD & TRIGGER ASHLEYS QC UPSTREAM MODULE
        ashleys_pipeline_only=False \ # STOP AFTER ASHLEYS QC - VALIDATION PURPOSE
        multistep_normalisation=True \ # TRIGGER MARCO'S MULTISTEP NORMALISATION
        MultiQC=True \ # TRIGGER samtools stats, FastQC & MultiQC reporting
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --cores 24 \
    --report TEST_DATA_REPORT.zip \
    --report-stylesheet workflow/report/custom-stylesheet.css
```

Questions???




SCNOVA

mkdir -p .tests/data_CHR17/RPE-BM510/scNOVA_input_user
awk 'BEGIN {FS=OFS="\t"} NR==1 {print "Filename", "Subclonality"} NR>1 && $2==1 {sub(/\.sort\.mdup\.bam/, "", $1); print $1, "clone"}' .tests/data_CHR17/RPE-BM510/cell_selection/labels.tsv > .tests/data_CHR17/RPE-BM510/scNOVA_input_user/input_subclonality.txt

```bash
snakemake \
    --cores 6 \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \ # DATA LOCATION
        ashleys_pipeline=True \ # DOWNLOAD & TRIGGER ASHLEYS QC UPSTREAM MODULE
        ashleys_pipeline_only=False \ # CONTINUES AFTER ASHLEYS QC - VALIDATION PURPOSE
        multistep_normalisation=True \ # TRIGGER MARCO'S MULTISTEP NORMALISATION
        MultiQC=True \ # TRIGGER samtools stats, FastQC & MultiQC reporting
        scNOVA=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --cores 24
```


########################################################################


### ⚙️ 3. MosaiCatcher execution (without preprocessing)

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
|   `-- bam
|       |-- Cell_01.sort.mdup.bam
|       |-- Cell_02.sort.mdup.bam
|       |-- Cell_03.sort.mdup.bam
|       `-- Cell_04.sort.mdup.bam
|
`-- Sample_2
    `-- bam
        |-- Cell_21.sort.mdup.bam
        |-- Cell_22.sort.mdup.bam
        |-- Cell_23.sort.mdup.bam
        `-- Cell_24.sort.mdup.bam
```

In a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Your Strand-seq BAM files of this sample go into the following folder:

- `bam` for the total set of BAM files

> Using the classic behavior, cells flagged as low-quality will only be determined based on coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process).

#### Old behavior

Previous version of MosaiCatcher (version ≤ 1.5) needed not only a `all` directory as described above, but also a `selected` folder (now renamed `bam`), presenting only high-quality selected libraries wished to be processed for the rest of the analysis.

You can still use this behavior by enabling the config parameter either by the command line: `=True` or by modifying the corresponding entry in the config/config.yaml file.

Thus, in a `Parent_Folder`, create a subdirectory `Parent_Folder/sampleName/` for each `sample`. Your Strand-seq BAM files of this sample go into the following folder:

- `all` for the total set of BAM files
- `selected` for the subset of BAM files that were identified as high-quality (either by copy or symlink)

Your `<INPUT>` directory should look like this:

```bash
Parent_folder
|-- Sample_1
|   |-- bam
|   |   |-- Cell_01.sort.mdup.bam
|   |   |-- Cell_02.sort.mdup.bam
|   |   |-- Cell_03.sort.mdup.bam
|   |   `-- Cell_04.sort.mdup.bam
|   `-- selected
|       |-- Cell_01.sort.mdup.bam
|       `-- Cell_04.sort.mdup.bam
|
`-- Sample_2
    |-- bam
    |   |-- Cell_01.sort.mdup.bam
    |   |-- Cell_02.sort.mdup.bam
    |   |-- Cell_03.sort.mdup.bam
    |   `-- Cell_04.sort.mdup.bam
    `-- selected
        |-- Cell_03.sort.mdup.bam
        `-- Cell_04.sort.mdup.bam
```

> Using the `input_bam_legacy` parameter, cells flagged as low-quality will be determined both based on their presence in the `selected` folder presented above and on coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process).

---

**⚠️ Warning**

Using the `input_bam_legacy` parameter, only **intersection** between cells present in the selected folder and with enough coverage will be kept. Example: if a library is present in the selected folder but present a low coverage [see Note here](#note:-filtering-of-low-quality-cells-impossible-to-process), this will not be processed.

---

#### 3C. FASTQ input & Preprocessing module\*

From Mosaicatcher version ≥ 1.6.1, it is possible to use [ashleys-qc-pipeline preprocessing module](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) as part of MosaiCatcher. To do so, the user need to enable the config parameter `ashleys_pipeline=True` and create a directory respecting the structure below (based on the model used for BAM inputs):

```bash
Parent_folder
|-- Sample_1
|   `-- fastq
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
    `-- fastq
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

- `ashleys_pipeline=True` and `input_bam_legacy=True` are mutually exclusive
- We advice to limit to the number of 2 the "_" characters in your file names (Ashleys-qc ML tool will react weidly with more than 3 "_")

---

### ⚡️ 4. Run the pipeline

#### Local execution (without batch scheduler) using conda only

After defining your configuration, you can launch the pipeline the following way:

```bash
snakemake \
    --cores <N> \
    --config \
        data_location=<INPUT_FOLDER> \
    --profile workflow/snakemake_profiles/local/conda/
```

#### Local execution (without batch scheduler) using singularity X conda (recommanded)

After defining your configuration, you can launch the pipeline the following way:

```bash
snakemake \
    --cores <N> \
    --config \
        data_location=<INPUT_FOLDER> \
    --profile workflow/snakemake_profiles/local/conda_singularity/ --singularity-args "-B /<mouting_point>:/<mounting_point>"
```

---

**ℹ️ Note**

It is possible to provide multiple mouting points between system and cointainer using as many `-B` as needed in the `singularity-args` command like the following: "-B /<mouting_point1>:/<mounting_point1> -B /<mouting_point2>:/<mounting_point2>"
For EMBL users, you don't need to specify this as this is already part of the execution profile (workflow/snakemake_profiles/HPC/slurm_EMBL)

---

---

**ℹ️ Note**

It is recommended to first run the command and check if there are any anomalies with the `--dry-run` option

---

---

**⚠️ Warning**

If you are experiencing any issues with conda-frontend snakemake option, please use `--conda-frontend conda` instead of `mamba`

---

#### HPC execution

MosaiCatcher can be executed on HPC using [Slurm](https://slurm.schedmd.com/documentation.html) by leveraging snakemake profile feature. Current Slurm profile [`workflow/snakemake_profiles/HPC/slurm_EMBL/config.yaml`] was defined and tested on EMBL HPC cluster but can be modified, especially regarding **partition** setting.

##### Current strategy to solve HPC job OOM

Workflow HPC execution usually needs to deal with out of memory (OOM) errors, out of disk space, abnormal paths or missing parameters for the scheduler. To deal with OOM, we are currently using snakemake restart feature that allows to automatically double allocated memory to the job at each attempt (limited to 6 for the moment). Then, if a job fails to run with the default 1GB of memory allocated, it will be automatically restarted tith 2GB at the 2nd attempt, 4GB at the 3rd, etc.

To execute MosaiCatcher on HPC, use the following command.

##### Command

```bash
snakemake \
    --config \
        data_location=<INPUT_FOLDER> \
    --singularity-args "-B /<mounting_point>:/<mounting_point>" \
    --profile workflow/snakemake_profiles/HPC/slurm_generic/
```

The `logs` and `errors` directory will be automatically created in the current directory, corresponding respectively to the `output` and `error` parameter of the `sbatch` command.

### 📊 5. Generate report [Optional]

Optionally, it's also possible to generate a static interactive HTML report to explore the different plots produced by MosaiCatcher, using the same command as before and just adding at the end the following:
`--report <report>.zip --report-stylesheet workflow/report/custom-stylesheet.css`, which correspond to the following complete command:

```bash
snakemake \
    --config \
        data_location=<INPUT_FOLDER> \
    --singularity-args "-B /<mounting_point>:/<mounting_point>" \
    --profile workflow/snakemake_profiles/slurm_generic/ \
    --report <report>.zip --report-stylesheet workflow/report/custom-stylesheet.css
```

---

**ℹ️ Note**

The zip file produced can be heavy (~1GB for 24 HGSVC samples ; 2000 cells) if multiple samples are processed in parallel in the same output folder.

---

## ArbiGent mode of execution

From 1.9.0, it's now possible to run MosaiCatcher in order to genotype a given list of positions specified in a bed file. To do so, `arbigent` config parameter need to be set to `True`.
Thus, an alternative branch of the pipeline will be executed instead of the classic one, targetting results produced by ArbiGent.

```bash
snakemake \
    --cores <N> \
    --config \
        data_location=<INPUT_DATA_FOLDER> \
        arbigent=True \
    --profile workflow/snakemake_profiles/local/conda_singularity/

```

A generic BED file is provided in `workflow/data/arbigent/manual_segmentation.bed`. A custom BED file can be specified by modifying `arbigent_bed_file` config parameter to the path of your choice.

---

**ℹ️ Note**

If you modify the chromosome list to be processed (remove chrX & chrY for instance), a dedicated rule will extract only from the BED file, the rows matching the list of chromosomes to be analysed.

---

## scNOVA mode of execution

From 1.9.2, it's now possible to run [scNOVA](https://github.com/jeongdo801/scNOVA/) directly from MosaiCatcher in order to determine the nucleosome occupancy associated to the SV calls provided during MosaiCatcher execution.

To do so, mosaicatcher must be executed during a first step in order to generate the SV calls and associated plots/tables required to determine subclonality.

Once the subclonality was determined, a table respecting the following template must be provided at this path `<DATA_LOCATION>/<SAMPLE>/scNOVA_input_user/input_subclonality.txt` where the different clones are named `clone<N>`:

| Filename             | Subclonality |
| -------------------- | ------------ |
| TALL3x1_DEA5_PE20406 | clone2       |
| TALL3x1_DEA5_PE20414 | clone2       |
| TALL3x1_DEA5_PE20415 | clone1       |
| TALL3x1_DEA5_PE20416 | clone1       |
| TALL3x1_DEA5_PE20417 | clone1       |
| TALL3x1_DEA5_PE20418 | clone1       |
| TALL3x1_DEA5_PE20419 | clone1       |
| TALL3x1_DEA5_PE20421 | clone1       |
| TALL3x1_DEA5_PE20422 | clone1       |

Once done, you can run the exact same command as previously and set `scNOVA` config parameter to `True` like the following:

```bash
snakemake \
    --cores <N> \
    --config \
        data_location=<INPUT_DATA_FOLDER> \
        scNOVA=True \
    --profile workflow/snakemake_profiles/local/conda_singularity/

```

---

**ℹ️ Note**

scNOVA related snakemake rules and scripts only support conda execution at the moment

---

## Multistep normalisation

From 2.1.0, it's now possible to use multistep normalisation (library size normalisation, GC correction, Variance Stabilising Transformation) resulting counts file, not only for visualisation purpose, but also to rely on it during the SV calling framework. To do so

```bash
snakemake \
    --cores <N> \
    --config \
        data_location=<INPUT_DATA_FOLDER> \
        multistep_normalisation=True \
        multistep_normalisation_for_SV_calling=True \
    --profile workflow/snakemake_profiles/local/conda_singularity/

```

## scTRIP multiplot

From 2.2.2, scTRIP multiplot (from Marco Cosenza) is now compatible with MosaiCatcher. The single requirement is to clone scTRIP multiplot repository (please reach out Marco if you want to access the repository, currently private) inside `workflow/scripts/plotting/scTRIP_multiplot`.

By default, scTRIP multiplot is set to false, to enable it, please update the `config/config.yaml` file or use `scTRIP_multiplot=True` in the config section of the command line. An example of a scTRIP multiplot is available [here](/docs/output.md#sctrip-multiplot-marco-cosenza)

---

**ℹ️ Note**

`multistep_normalisation_for_SV_calling` is mutually exclusive with `hgsvc_based_normalisation`

---

## Pipeline update procedure

If you already use a previous version of mosaicatcher-pipeline, here is a short procedure to update it:

- First, update all origin/<branch> refs to latest:

`git fetch --all`

- Jump to a new version (for example 2.2.2) & pull code:

`git checkout 2.2.2 && git pull`

Then, to initiate or update git snakemake_profiles submodule:

`git submodule update --init --recursive`
