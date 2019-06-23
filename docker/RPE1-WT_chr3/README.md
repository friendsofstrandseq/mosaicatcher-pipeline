# Docker image including an example data set

This is the same as `smei/mosaicatcher-pipeline` but including the reference 
genome *hg38*, a single-chromosome data set (chr3 from RPE1-WT data) and
`snakemake`.

**This can be used to run a complete example pipeline directly out of Docker**.


### Availability

* On [dockerhub](https://cloud.docker.com/repository/docker/smei/mosaicatcher-pipeline-rpe1-chr3)

### Updates

This image is *nor* configured for auto-build and likely won't be updated in the future. 

### Access to RPE1-WT data

* The chr3-only dataset is baked into this image, so you can extract it from here.

* The whole-genome data set can be downloaded from ENA. The file names are listed
  in `ena_files.txt`. Here is an example command of how they can be downloaded:

  ```
  mkdir example_data \
    && cd example_data \
    && for x in $(grep RPE1WT ena_files.txt); \
      do \
        echo "downloading ${x}"; \
        wget -q ${x%.bai}; \
        wget -q ${x}; \
      done
  ```
  The files should look like this:
  ```
  $ ls 
  RPE1WTPE20430.sort.mdup.bam.bai  RPE1WTPE20464.sort.mdup.bam.bai
  RPE1WTPE20401.sort.mdup.bam      RPE1WTPE20431.sort.mdup.bam      RPE1WTPE20465.sort.mdup.bam
  RPE1WTPE20401.sort.mdup.bam.bai  RPE1WTPE20431.sort.mdup.bam.bai  RPE1WTPE20465.sort.mdup.bam.bai
  (...)
  ```
  **Note** The dataset on ENA only contains the 80/96 high-quality cells.