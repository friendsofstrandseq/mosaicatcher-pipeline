# Mosaicatcher workshop

## Technical prerequisites

- SSHFS/SFTP connection to visualise/download/access files created (WinSCP/FileZilla/Cyberduck)
- Functional terminal connected to the EMBL cluster (if not follow SSH key configuration here: https://www.embl.org/internal-information/it-services/hpc-resources-heidelberg/)
- Have a workspace on /g/korbel

## Workshop prerequisites

- Select a sample you sequenced in the past to be processed using MosaiCatcher
- Download this MosaiCatcher report: https://oc.embl.de/index.php/s/WBgrzBjyzdYdVJA/download

## EMBL cheatsheet

### connect to seneca

```
ssh USERNAME@seneca.embl.de
```

### connect to login nodes

```
ssh USERNAME@login0[1,2,3,4].embl.de (login01 to login04)
```

## Snakemake cheat sheet & important things

Snakemake is a workflow manager, it will handle execution and distribute computing based on configuration (local execution: use local CPUs; cluster execution: will generate the sbatch commands for each job). Under cluster execution, snakemake is using 1 single instance locally, the other jobs computed on the cluster.

To run the pipeline, go into the repository and run: snakemake

### Snakemake important arguments/options

- --dry-run
- --profile
- --config
- --rerun-triggers
- --touch

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
- scNOVA outputs: PARENT_FOLDER/SAMPLE_NAME/scNOVA_result/

## CLI usage of the pipeline

### Quick Start

---

**ℹ️ Important Notes**

- From 2.2.0, you don't need to clone both [ashleys-qc-pipeline preprocessing module](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) and [mosaicatcher-pipeline](https://github.com/friendsofstrandseq/mosaicatcher-pipeline). By using `ashleys_pipeline_only=True` combined with `ashleys_pipeline=True` in the configuration of MosaiCatcher, this will stop the execution after the generation of the files related to ashleys-qc-pipeline. This allow you to use a single pipeline, repository and container and limit the potential conflicts by processing the same folder (data_location) by different repositories and set of files (including the workflow/data/ref_genomes references files).

- Configuration is crucial in the pipeline, via command line or via YAML file, it will define where to stop, which mode, which branch, which options to be used.

- The profile parameter will define if you execute locally/on the cluster the jobs and also defines the dependencies management system (conda + container / just conda)

---

1. Clone the repository & update to latest profiles

```bash
git clone --recurse-submodules https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git && cd mosaicatcher-pipeline && git submodule update --remote
```

For download time and space purpose in that workshop, let symlink the latest version of the pipeline container. If you're updating the pipeline or if there is any issue in the future, you can just delete the symbolic link and snakemake will download the real image in your local workspace.

```bash
ln -s /g/korbel2/weber/workspace/mosaicatcher-update/.snakemake/singularity/43cee69ec12532570abaa913132bfaa5.simg .snakemake/singularity/
```

Give a look at the folder structure:

```bash
tree -h .tests/data_CHR17
```

Similar to below (multiple samples can be processed in the same parent_folder)

```
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
```

2. Load snakemake

```bash
module load snakemake/7.32.4-foss-2022b
```

3. Run on example data on only one small chromosome.

List all the options available in MosiiCatcher.

```bash
snakemake --config list_commands=True --core 1 --dry-run
```

Then dry-run on the test dataset.

```bash
snakemake \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \
        ashleys_pipeline=True \
        ashleys_pipeline_only=True \
        multistep_normalisation=True \
        MultiQC=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 20 \
    --dry-run
```

If no error message, you are good to go!

```bash
snakemake \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \
        ashleys_pipeline=True \
        ashleys_pipeline_only=True \
        multistep_normalisation=True \
        MultiQC=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 20
```

4. Let's get a look at the important files

```bash
cat .tests/data_CHR17/RPE-BM510/counts/RPE-BM510.info_raw # Counts stats
zcat .tests/data_CHR17/RPE-BM510/counts/RPE-BM510.txt.raw.gz | less # Counts
cat .tests/data_CHR17/RPE-BM510/cell_selection/labels.tsv # Ashleys predictions
```

Now using the SSHFS/SFTP, get a look at the plots generated at that location:

```
.tests/data_CHR17/RPE-BM510/plots
```

5. Continue with SV calling using MosaiCatcher

We now set `ashleys_pipeline_only=False` and the pipeline will use the cells labeled as correct by the model (or by you if you changed some predictions) in the file (.tests/data_CHR17/RPE-BM510/cell_selection/labels.tsv).

```bash
snakemake \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \
        ashleys_pipeline=True \
        ashleys_pipeline_only=False \
        multistep_normalisation=True \
        MultiQC=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 20
```

---

Notes:

- The list of jobs in the jobs table doesn't mandatory represent the complete list, as some rule are "checkpoint" rules, the list of jobs to be executed is re-evaluated on the fly after the completion of those checkpoints.
- If you trust blindly ashleys-qc and just want to trigger the pipeline from the FASTQ processing to the SV calling, you can directly run the command below without doing the step 3 above.

---

6. Once done, generate a report using the **same command** and add the --report option at the end

```bash
snakemake \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17
        ashleys_pipeline=True
        ashleys_pipeline_only=False \
        multistep_normalisation=True \
        MultiQC=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 20 \
    --report TEST_DATA_REPORT.zip \
    --report-stylesheet workflow/report/custom-stylesheet.css
```

Using your SFTP tool, download the .zip file on your laptop, extract and check the content!

7. Bonus for test dataset: run scNOVA!

---

Note:

scNOVA mode in MosaiCatcher is considered as a wrapper from the original pipeline. If you are facing an issue, please contact Hyobin Jeong who developed the pipeline.

---

Let's now first prepare the subclonality input file:

```bash
mkdir -p .tests/data_CHR17/RPE-BM510/scNOVA_input_user # TO CREATE THE CORRESPONDING FOLDER

# THEN USE ASHLEYS LABELS AS A TEMPLATE TO GENERATE input_subclonality.txt
awk 'BEGIN {FS=OFS="\t"} NR==1 {print "Filename", "Subclonality"} NR>1 && $2==1 {sub(/\.sort\.mdup\.bam/, "", $1); print $1, "clone"}' .tests/data_CHR17/RPE-BM510/cell_selection/labels.tsv > .tests/data_CHR17/RPE-BM510/scNOVA_input_user/input_subclonality.txt
```

---

Note:

You can use the block of commands above as a template for your own analysis in the future if you want to use scNOVA. Thus, you will limit discrepancies and potential issues in the pipeline.

You can now trigger scNOVA mode by adding `scNOVA=True` to the config section:

---

```bash
snakemake \
    --configfile .tests/config/simple_config.yaml \
    --config \
        data_location=.tests/data_CHR17 \
        ashleys_pipeline=True \
        ashleys_pipeline_only=False \
        multistep_normalisation=True \
        MultiQC=True \
        scNOVA=True \
        chromosomes_to_exclude="[chrY]"
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 20
```

---

Note:

chrY is not handled yet by scNOVA and need to discarded

---

8. Let's now use a real sample

Using the sample name (in 2023) you picked, you can process the corresponding data quite simply with MosaiCatcher.

genecore: enable genecore mode that will symlink & reformat the data
genecore_prefix: define location where the date_folder will be picked up
genecore_date_folder: folder corresponding to the date + flowcell used
genecore_regex_element: index system used (PE20 or iTRU)

```bash
snakemake \
    --cores 6 \
    --config \
        genecore=True \
        genecore_prefix=/g/korbel/STOCKS/Data/Assay/sequencing/2023 \
        genecore_date_folder=2023-XX-XX-XXXXXX \
        samples_to_process="[YOUR_SAMPLE]" \
        data_location=DATA_LOCATION_XXX \
        ashleys_pipeline=True \
        ashleys_pipeline_only=True \
        multistep_normalisation=True \
        MultiQC=True \
    --profile workflow/snakemake_profiles/HPC/slurm_EMBL/ \
    --jobs 100
```

If you don't have pick a sample yet, you can select one in the list below:

- 2023-03-24-HCNJ5AFX5 -- TAllPDX6340p6RELs1p1x
- 2023-03-24-HCNJ5AFX5 -- IMR90E6E7PD106s1p3x01
- 2023-03-24-HCNJ5AFX5 -- IMR90E6E7PD103s1p2x01
- 2023-03-24-HCNGCAFX5 -- BAB3161
- 2023-03-24-HCNGCAFX5 -- BAB3114
- 2023-03-24-HCNGCAFX5 -- BAB14547
- 2023-03-08-HCNGHAFX5 -- HGSVCpool2inWell5ul
- 2023-03-08-HCNGHAFX5 -- HGSVCpool2inWell2ul
- 2023-03-08-HCNGHAFX5 -- HGSVCpool2OPSfromFrozen2ul
- 2023-02-08-HCN3VAFX5 -- HGSVCpool2
- 2023-01-24-H75WVAFX5 -- TAllPDX4973p8RELx
- 2023-01-24-H75WVAFX5 -- RPEH2BdendraMicroNx01
- 2023-01-24-H75WVAFX5 -- IMR90E6E7PD106x01

Tip tool:

```
# make sure your sample name is exact & without typo
SAMPLE="YOUR_SAMPLE" && find /g/korbel/STOCKS/Data/Assay/sequencing/2023 -name "*$SAMPLE*" | head -n 1 | tr "/" "\t" | cut -f 9
```

## Strand-Scape usage

Strand-Scape current adress: http://seneca.embl.de:8060

If you want to use Strand-Scape and further reprocess data available there, you can pick the location of raw data on /scratch. If you do so, please use the following while copying over your workspace:

```
cp --preserve=timestamps STRAND_SCAPE_DATA YOUR_LOCATION
```

Snakemake engine is file-based and its default behavior is to look at various elements in order provide a reproducible analysis, meaning that if: snakemake version, container, software dependencies, timestamps changed, the system will want to recompute everything in order to provide an up to date result. To bypass this, we can require snakemake to run with the following parameter `--rerun-triggers mtime`, meaning that snakemake will only look at timestamps and will expect that the files timestamps chronological order is matching the graph computed by the pipeline (input output rationale). By using the command above, you should preserve the timestamps used during the original execution, and thus, limit any further issues.

Strand-scape provides a feature (cell selection menu) to allow you to save your corrected predictions, by looking simultaneously at the report generated by the pipeline. If you want to retrigger a second pipeline execution yourself after using your expertise, you can use the parameter `--strandscape_labels_path` and specify the location of your version of the cell_selection/labels.tsv file which should be named the following (strandscape_labels-USERNAME.tsv). Thus, MosaiCatcher will bypass the default cell_selection/labels.tsv file and use your instead.

## Other elements to keep in consideration

- Please verify you consulted the [documentation](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/blob/master/docs/usage.md) that contains additional information.
- Snakemake logs are colored in yellow/green/red. If you notice that some jobs are failing while using the cluster profile (slurm_EMBL), that might not indicate directly that the pipeline is likely to crash. It can be to Out Of Memory (OOM) issue where a job is missing of memory. MosaiCatcher provides a system that restart automatically jobs by doubling the memory at each new attempt. Also, another issue related to our BeeGFS shared file system is well-known to lead to asynchronicity problems where the next job is expecting an input to be available while the output of the previous job is not compltely written yet on the disk. This is currently fixed by the community.
- [Pipeline update procedure](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/blob/master/docs/usage.md#pipeline-update-procedure)
