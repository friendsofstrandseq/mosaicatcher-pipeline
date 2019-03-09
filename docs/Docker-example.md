# Docker workflow (including example data)

We provide a Docker image ([mosaicatcher-pipeline-rpe-1](https://hub.docker.com/r/smei/mosaicatcher-pipeline-rpe-1)) that was made specifically to re-run a complete analysis on the epethelial cell line data set RPE-1.

This image includes all required software + data and is ~30GB large.


## How to re-run RPE-1 analysis

Run Docker as shown below. The image will be fetched automatically from Dockerhub.

```
sudo docker run -v $(pwd):/host -ti smei/mosaicatcher-pipeline-rpe-1 bash
root@70768001ace0:/pipeline#
```

Now, from within the Docker container, run snakemake with the `Snake.config-singularity-rpe1.json` config file (it differs from `Snake.config-singularity.json` only by specifying the file path of SNV sites).

```
snakemake --configfile Snake.config-singularity-rpe1.json
```

This should reproduce the whole analysis. Note that certain steps, for example *FreeBayes* may take a long time.


### Extracting results from within Docker

Data within Docker containers is not persistent. To extract the results of running this pipeline, you can bind mount a folder of the host system into the container, as shown above using `-v`.

In this case, just copy the relevant files to the host directory, for example like this

```
root@70768001ace0:/pipeline# cp -r postprocessing /host/
```

An alternative to bind mounts are [Docker volumes](https://docs.docker.com/storage/volumes/).


## Explanation of how this container was created

The RUN commands within the [Dockerfile](../RPE-1/Dockerfile) describe exaclty how data is added to the workflow. Here is a step-by-step explanation:

```
FROM smei/mosaicatcher-pipeline
```

This container is based on the image used for Singularity (explained [here](Singularity.md)) - it hence already contains all required software.


**Add SNV calls**

At first, a set of SNV calls is prepared, which will be re-genotyped within the given RPE-1 sample. To this end 1000 Genomes SNV calls are downloaded per chromosome using `wget`:

```
RUN mkdir snv_sites_to_genotype \
    && mkdir tmp \
    && cd tmp \
    && for x in $(seq 1 22) X; \
	do \
	    echo "${x} chr${x}" >> rename-chrs; \
	    echo "downloading chromosome $x"; \
	    wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${x}_GRCh38_sites.20170504.vcf.gz; \
	    echo "ALL.chr${x}_GRCh38_sites.20170504.vcf.gz" >> file-list; \
	done \

```

Then, VCF files are concatenated into a single file...

```
    && bcftools concat -Oz -o tmp.vcf.gz -f file-list \
```

...and chromosome names are changed to match the names of the reference genome.

```
    && bcftools annotate --rename-chrs rename-chrs tmp.vcf.gz \
        | bgzip \
        > ../snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz \
```

At last, the VCF is indexed and remaining files are removed.

```
    && tabix ../snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz \
    && cd .. \
    && rm -rf tmp
```


**Add single-cell BAM files**

```
RUN mkdir -p bam/RPE1-WT/all \
    && cd bam/RPE1-WT/all \
    && for x in $(grep RPE1WT ../../../docker/ena_files.txt); \
        do \
        echo "downloading ${x}"; \
        wget -q ${x}; \
        wget -q ${x%.bai}; \
        done
```

Here, the single-cell BAM files are fetched from the European Nucleotide Archive.

```
RUN mkdir -p bam/RPE1-WT/selected \
    && cd bam/RPE1-WT/selected \
    && for x in 01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 19 21 22 23 \
                24 25 26 28 29 30 31 32 34 35 36 37 38 39 40 41 42 45 47 48 49 \
                50 51 52 53 54 55 56 57 59 62 63 64 65 66 67 68 69 70 71 72 73 \
                74 76 77 78 80 81 83 84 85 86 89 90 91 92 93 94 95; \
        do \
        ln /pipeline/bam/RPE1-WT/all/RPE1WTPE204${x}*.bam; \
        ln /pipeline/bam/RPE1-WT/all/RPE1WTPE204${x}*.bam.bai; \
        done
```

Then, the set of "good" cells is selected and placed within the `selected` folder. Choosing good cells requires knowledge of Strand-seq and is described in [Sanders *et. al.*, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28492527).

**Add the reference genome**

```
RUN cd /\
    && wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    && gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    && mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference.fa \
    && samtools faidx reference.fa
```

At last, the refernce genome is taken from NCBI and added as `/reference.fa`.

```
RUN cd /tmp \
    && wget -q https://bioconductor.org/packages/3.6/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.1.tar.gz \
    && tar -xf BSgenome.Hsapiens.UCSC.hg38_1.4.1.tar.gz  \
    && mv BSgenome.Hsapiens.UCSC.hg38/inst/extdata/single_sequences.2bit /usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit \
    && rm -rf BSgenome.Hsapiens.UCSC.hg38*
```

Also, the R-version of the reference genome (*BSgenome.Hsapiens.UCSC.hg38*) is provided. This was previously stripped off the image to save disk space.
