# Whoeps 16th Feb 2021
# Snakemake for bringing bams into mosaicatcher-complying format.

# This snakefile contains rules to:
# 1) exchange the SM and ID tag,
# 2) Index the resulting bam files
# 3) Make automatically symlinks from the /all/ to /selected/ folders.

# According to my labbook (check key [BAM sample prep]), here is what bams need to fulfill:

##############################
# 'ID' has to be the same in
# 1) The bam name
# 2) The RG ID thing in the header
# 3) In the end of each individual read.
#
# 'SM' has to be the same in
# 1) The folder name in /bam/
# 2) The header RG SM thing
# 3) The bed file
##############################

from collections import defaultdict

Testmode = False

path_to_orig_samples = "/scratch/tweber/SCO_COURSE/HJ_MIXTURE_RPE1_Mix_renamed2"


SAMPLE, BAM = glob_wildcards(path_to_orig_samples + "/{sm}/raw/{id}.sort.mdup.bam")
SAMPLES = sorted(set(SAMPLE))
print(SAMPLE)
print(BAM)
# print(ONEKG)

### THIS PART IS STOLEN FROM MOSAICATCHER ###
CELL_PER_SAMPLE = defaultdict(list)
BAM_PER_SAMPLE = defaultdict(list)
for sample, bam in zip(SAMPLE, BAM):
    BAM_PER_SAMPLE[sample].append(bam)
    CELL_PER_SAMPLE[sample].append(bam.replace("_sorted", ""))

ALLBAMS_PER_SAMPLE = BAM_PER_SAMPLE
print(BAM_PER_SAMPLE)
print(ALLBAMS_PER_SAMPLE)

print("Detected {} samples:".format(len(SAMPLES)))
for s in SAMPLES:
    print(
        "  {}:\t{} cells\t {} selected cells".format(
            s, len(ALLBAMS_PER_SAMPLE[s]), len(BAM_PER_SAMPLE[s])
        )
    )
#################################33

# Targets #
bams_all = []
bais_all = []
bams_select = []
bais_select = []
for s in SAMPLES:
    bams_all.extend(
        expand(
            "{path}/{SM}/bam/{ID}.sort.mdup.bam",
            path=path_to_orig_samples,
            SM=s,
            ID=sorted(ALLBAMS_PER_SAMPLE[s]),
        )
    )
# bais_all.append(
#     expand(
#         "{path}/{SM}/all/{ID}.sort.mdup.bam.bai",
#         path=path_to_orig_samples,
#         SM=s,
#         ID=ALLBAMS_PER_SAMPLE[s],
#     )
# )
# bams_select.append(
#     expand(
#         "{path}/{SM}/selected/{ID}.bam",
#         path=path_to_orig_samples,
#         SM=s,
#         ID=ALLBAMS_PER_SAMPLE[s],
#     )
# )
# bais_select.append(
#     expand(
#         "{path}/{SM}/selected/{ID}.bam",
#         path=path_to_orig_samples,
#         SM=s,
#         ID=ALLBAMS_PER_SAMPLE[s],
#     )
# )

# bams_all = ['HG00513/all/HG00513_IV_045.bam']
rule all:
    input:
        bams_all,
        # bais_all,
        # bams_select,
        # bais_select,


rule change_id_and_sam:
    input:
        # bam_orig=expand(
        #     "{path}/{SM}/raw/{ID}.bam",
        #     zip,
        #     path=path_to_orig_samples,
        #     SM=SAMPLE,
        #     ID=BAM,
        # ),
        bam_orig="{path}/{SM}/raw/{ID}.sort.mdup.bam",
    output:
        bam_out="{path}/{SM}/bam/{ID}.sort.mdup.bam",
    envmodules:
        "SAMtools/1.14-GCC-11.2.0"
    resources:
        mem_mb="16000",
        time="10:00:00",
    shell:
        """
        # old_id=$(samtools view -H {input.bam_orig} | grep "^@RG" | grep -P -o "\tID:[A-Za-z0-9_\-]*" | sed 's/ID://g')
        # old_sm=$(samtools view -H {input.bam_orig} | grep "^@RG" | grep -P -o "\tSM:[A-Za-z0-9_\-]*" | sed 's/SM://g')
        # echo "{wildcards.ID} {wildcards.SM} $old_id $old_sm"
        # First, the 'ID' tag
        samtools view -H {input.bam_orig} | sed "s/ID:[A-Za-z0-9_-]*/ID:{wildcards.ID}/g;s/SM:[A-Za-z0-9_-]*/SM:{wildcards.SM}/g;s/RG:Z:[A-Za-z0-9_-]*/RG:Z:{wildcards.ID}/g" > {output.bam_out}.header
        samtools view -h {input.bam_orig} | sed "s/ID:[A-Za-z0-9_-]*/ID:{wildcards.ID}/g;s/SM:[A-Za-z0-9_-]*/SM:{wildcards.SM}/g;s/RG:Z:[A-Za-z0-9_-]*/RG:Z:{wildcards.ID}/g" | samtools view -bS > {output.bam_out}.core
        samtools reheader -P {output.bam_out}.header {output.bam_out}.core > {output.bam_out}
        # Remove intermediate files
        rm {output.bam_out}.header {output.bam_out}.core
        """


# rule add_idx:
#     input:
#         bam="{path}/{SM}/all/{ID}.sort.mdup.bam",
#     output:
#         bai="{path}/{SM}/all/{ID}.sort.mdup.bam.bai",
#     shell:
#         """
#         samtools index {input.bam}
#         """
