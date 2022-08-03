# General settings
To configure this workflow, modify config/config.yaml according to your needs, following the explanations provided in the file.

# Samples

At this moment, mosaicatcher-pipeline is designed to automatically produce a configuration file saved in the input folder selected based on the parameter <input_bam_location>.
This configuration file looks like the following when running the pipeline with [ashleys-qc-pipeline](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) preprocessing module **disabled** (BAM input):


| File           | Folder            | Sample    | Cell           | Full_path                                                    |
| -------------- | ----------------- | --------- | -------------- | ------------------------------------------------------------ |
| BM510x3PE20401 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20401 | .tests/data_CHR21/RPE-BM510/all/BM510x3PE20401.sort.mdup.bam |
| BM510x3PE20402 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20402 | .tests/data_CHR21/RPE-BM510/all/BM510x3PE20402.sort.mdup.bam |


If the [ashleys-qc-pipeline](https://github.com/friendsofstrandseq/ashleys-qc-pipeline) preprocessing module is **enabled** (FASTQ input), thus the configuration file looks like the following:

| File             | Folder            | Sample    | Cell           | Full_path                                                   |
| ---------------- | ----------------- | --------- | -------------- | ----------------------------------------------------------- |
| BM510x3PE20401.1 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20401 | .tests/data_CHR21/RPE-BM510/fastq/BM510x3PE20401.1.fastq.gz |
| BM510x3PE20401.2 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20401 | .tests/data_CHR21/RPE-BM510/fastq/BM510x3PE20401.2.fastq.gz |
| BM510x3PE20402.1 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20402 | .tests/data_CHR21/RPE-BM510/fastq/BM510x3PE20402.1.fastq.gz |
| BM510x3PE20402.2 | .tests/data_CHR21 | RPE-BM510 | BM510x3PE20402 | .tests/data_CHR21/RPE-BM510/fastq/BM510x3PE20402.2.fastq.gz |

This configuration file is then used to create the dictionnary and lists needed to flag wildcards and identify the samples, the different cells for each of them and the associated files.  