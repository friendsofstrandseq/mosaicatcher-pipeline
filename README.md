# Strand-seq pipeline

> **Ongoing work**

Preliminary SV calling using Strand-seq data - summarized in a [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


### Bioconda environment
To install the correct environment, you can use Bioconda.

  1. **Install MiniConda:**
    In case you do not have Conda yet, it is easiest to just install [MiniConda](https://conda.io/miniconda.html).

  2. **Create environment:**
    ```
    conda env create -n strandseqnation -f conda-environment.yml
    source activate strandseqnation
    ```
    That's it, you are ready to go.

### How to use it

  1. **Install required software:**

    * Install [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher) (*currently you will need the `develop` branch*)
    * Get the R-scripts from [strandsequtils](https://github.com/friendsofstrandseq/strandsequtils)
    * Install BSgenome.Hsapiens.UCSC.hg38 (can be skipped of you use the Bioconda environment, see above):
      ```
      source("https://bioconductor.org/biocLite.R")
      biocLite('BSgenome.Hsapiens.UCSC.hg38')
      ```
    * [Strand-Phaser](https://github.com/daewoooo/StrandPhaseR) is installed automatically

  2. **Set up the configuration of the snakemake pipeline**

    * Open `Snake.config.json` and specify the path to the executatables
      (such as Mosaicatcher) and to the R scripts.
    * Create a subdirectory `bam/` and another subdirectory per sample (e.g.
      `bam/NA12878/`). **Multiple samples can be run together not**.
      Then copy (or soft-link) the Strand-seq single-cell libraries (one BAM
      file per cell) in there. Note that bam files need to be sorted and indexed,
      contain a read group and should have duplicates marked.

  3. **Run Snakemake**

    * run `snakemake` to compute all tasks locally
    * Alternatively, you can ask Snakemake to submit your jobs to a HPC cluster. To this end edit the `Snake.cluster.json` file according to your available HPC environment and call

      ```
      snakemake -j 100 \
        --cluster-config Snake.cluster.json \
        --cluster "???"
      ```
      
### SNV calls

  The pipeline will run simple SNV calling using [samtools](https://github.com/samtools/samtools)
  and [bcftools](https://github.com/samtools/bcftools). If you **already have
  SNV calls**, you can avoid that by entering your VCF files into the pipeline.
  To so, make sure the files are [tabix](https://github.com/samtools/tabix)-indexed
  and specifigy them inside the `Snake.config.json` file:
  ```
  "snv_calls"     : {
        "NA12878" : "path/to/snp/calls.vcf.gz"
    },
  ```
