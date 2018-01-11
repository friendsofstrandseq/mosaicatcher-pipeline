# Strand-seq pipeline

> **Ongoing work**

Preliminary SV calling using Strand-seq data - summarized in a [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.

### How to use it

  1. **Install required software:**
  
    * Install [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher) (*currently you will need the `develop` branch*)
    * Get the R-scripts from [strandsequtils](https://github.com/friendsofstrandseq/strandsequtils)
    * Install BSgenome.Hsapiens.UCSC.hg38:
      ```
      source("https://bioconductor.org/biocLite.R")
      biocLite('BSgenome.Hsapiens.UCSC.hg38')
      ```
    * [Strand-Phaser](https://github.com/daewoooo/StrandPhaseR) is installed automagically

  2. **Set up the configuration of the smakemake pipeline**
  
    * Open `Snake.config.json` and specify the path to the executatables (such as Mosaicatcher) and to the R scripts.
    * Create a subdirectory `bam` and copy (or soft-link) the Strand-seq single-cell libraries in there. Note that bam files need to contain a read group and should have duplicates marked.
   
  3. **Run Snakemake**

    * run `snakemake` to compute all tasks locally
    * Alternatively, you can ask Snakemake to submit your jobs to a HPC cluster. To this end edit the `Snake.cluster.json` file according to your available HPC environment and call 
      
      ```
      snakemake -j 100 \
        --cluster-config Snake.cluster.json \
        --cluster "???"
      ```
