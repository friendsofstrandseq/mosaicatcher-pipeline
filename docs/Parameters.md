
# Parameters

## MosaiCatcher arguments




```
mode
```

MosaiCatcher currently supports three different modes of execution : `count`, `segmentation` and `mosaiclassifier`.
- `count` (selected by default) will only performs `Mosaic count` binning and count reads for each bin produced
- `segmentation` will run the pipeline until the `Mosaic segmentation` and selection of the correct segments
- `mosaiclassifier` will run the complete pipeline until the detection of SV in each selected cell of the samples

To select your mode of execution, use the following argument `--config mode=[count|segmentation|mosaiclassifier]`


```
plot
```

For each of these modes, you can *enable* or *disable* the plots generation by using `--config plot=[True|False]`


```
check_sm_tag
```
Based on pysam, will compare for each BAM file, if the header SM tag is identical to the folder name in order to prevent further issues.

```
dl_bam_example
```
Allow to retrieve automatically BAM example data to run the pipeline.

```
dl_external_files
```
Allow to retrieve automatically external files (GRCh38 reference genome + 1000G SNV VCF file) required to run the pipeline.



## Snakemake arguments

Here are presented some essential snakemake options that could help you. 

```
--cores, -c
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
--conda-frontend [mamba|conda] 
```
Choose the conda frontend for installing environments. Mamba is much faster and highly recommended. Default: “mamba”


```
--use-singularity 
```
If defined in the rule, run job within a singularity container. If this flag is not set, the singularity directive is ignored.

```
--singularity-args "-B /mounting_point:/mounting_point"
```
Pass additional args to singularity. `-B` stands for binding point between the host and the container.

```
--dryrun, -n 
```
Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.

```
--rerun-incomplete, --ri
```
Re-run all jobs the output of which is recognized as incomplete.


```
--keep-going, -k
```
Go on with independent jobs if a job fails.

```
-T, --retries, --restart-times
```
Number of times to restart failing jobs (defaults to 0).

```
--forceall, -F
```
Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.

---
**ℹ️ Note**

Currently, the binding command needs to correspond to the mounting point of your system (i.e: "/tmp:/tmp").
On seneca for example (EMBL), use `"/g:/g"` if you are working on `/g/korbel[2]` or `"/scratch:/scratch"` if you plan to work on `scratch`.

---

Obviously, all other [snakemake CLI options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) can also be used. 


