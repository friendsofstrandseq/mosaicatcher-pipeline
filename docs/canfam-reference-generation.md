# CanFam Reference Files Generation Guide

## Required Files

### 1. Reference Genome (`canfam.fa`)
Download and filter to main chromosomes only (chr1-38, chrX, chrY), excluding scaffolds/unplaced contigs and chrM:

```bash
module load SAMtools/1.21-GCC-13.3.0
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz  # Or canFam4/canFam6
gzip -d canFam3.fa.gz
grep "^>" canFam3.fa | grep -v "[_M]" | sed 's/>//' | xargs samtools faidx canFam3.fa > canfam.fa
bwa index canfam.fa
```

### 2. GC Content Matrix (`canfam.GC_matrix.txt.gz`)
Generate 200 kb windows and compute %GC per bin:

```bash
module load BEDTools/2.30.0-GCC-12.2.0
module load SAMtools/1.21-GCC-13.3.0

# Create 200kb windows
samtools faidx canfam.fa
cut -f 1,2 canfam.fa.fai | grep -v "[_M]" > canfam.txt
bedtools makewindows -g canfam.txt -w 200000 > canfam.bed
bedtools getfasta -fi canfam.fa -bed canfam.bed > canfam.win.fa

# Count nucleotides and compute %GC
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faCount
chmod +x faCount
./faCount canfam.win.fa | grep -v "^total" > canfam_nt_content.txt

# Reformat and compute %GC
(grep "^#" canfam_nt_content.txt | sed 's/#seq\tlen/chrom\tstart\tend/;s/$/\t%GC/';
grep -v "^#" canfam_nt_content.txt | awk '
BEGIN { FS=OFS="\t" }
{ GC = $4 + $5; tot = $3 + $4 + $5 + $6; perc = tot ? 100*GC/tot : "NA"; split($1, arr, "[:-]") }
{ print arr[1], arr[2], arr[3], $3, $4, $5, $6, $7, $8, perc }') | bgzip > canfam.GC_matrix.txt.gz
```

### 3. Normalization/Blacklist File (`HGSVC.200000.txt`)
**Requires** 2 high-quality Strand-seq samples:
- **Control**: Diploid, no clonal SVs, high coverage (use XX karyotype to avoid chrX issues)
- **Test**: Any quality sample

```bash
# 1. Run mosaicatcher count on both samples (200kb windows)
mosaicatcher count --verbose --do-not-blacklist-hmm \
  -o control_sample.txt.raw.gz -i control_sample.info_raw \
  -x chroms_to_exclude.txt -w 200000 *.bam

mosaicatcher count --verbose --do-not-blacklist-hmm \
  -o test_sample.txt.raw.gz -i test_sample.info_raw \
  -x chroms_to_exclude.txt -w 200000 *.bam

# 2. Generate normalization/blacklist file
wget https://raw.githubusercontent.com/friendsofstrandseq/mosaicatcher/master/R/makeNorm.R
Rscript makeNorm.R control_sample.txt.raw.gz test_sample.txt.raw.gz HGSVC.200000.txt

# 3. Manually curate centromeric regions if needed
```

**Critical**: Use XX/diploid control to avoid chrX normalization issues (see `docs/mm39-notes.md`)

### 4. Binning BED File (`canfam.bin_200kb_all.bed`)
```bash
# Reuse the bed file from GC matrix generation
cp canfam.bed canfam.bin_200kb_all.bed
```

### 5. Segmental Duplications (Optional - Not Required)
Can be skipped for initial implementation. If needed later, download from UCSC database.

## File Locations in Repository

Add files to `mosaicatcher-referencedata` submodule:
```
workflow/data/
├── ref_genomes/canfam.fa{,.amb,.ann,.bwt,.pac,.sa,.fai}
├── GC/canfam.GC_matrix.txt.gz
└── normalization/canfam/HGSVC.200000.txt
```

## Testing (Can Be Done Later)

Lightweight test dataset can be created after reference files are ready:
1. Extract single chromosome (e.g., chr1) from sample BAMs/FASTQs
2. Test with config override: `chromosomes=[chr1]`
3. See `.tests/data_chr19_mm/` for example structure

## References

- mm39 generation example: `/g/korbel/shared/adelgado/mosaicatcher_reffiles_mm39/`
- mm39 caveats (chrX issues): `docs/mm39-notes.md`
