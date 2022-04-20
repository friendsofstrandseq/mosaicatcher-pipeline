import pysam
import sys, os
import pandas as pd

bam_file_path = "/g/korbel2/weber/MosaiCatcher_files/bam_KG_full/RPE1-WT/all/RPE1WTPE20401.sort.mdup.bam"
# bam_file_path = sys.argv[1]


# def check_bam_header(bam_file_path):
#     # GET BAM FILE HEADER WITH PYSAM
#     h = pysam.view("-H", bam_file_path)
#     h = [e.split("\t") for e in h.split("\n")]
#     sm_tag_list = list(set([sub_e.replace("SM:", "") for e in h for sub_e in e if "SM:" in sub_e]))

#     # FOLDER NAME BASED ON PATH
#     folder_name = bam_file_path.split("/")[-3]

#     # CHECKS
#     assert len(sm_tag_list) == 1, "Two different SM tags in the header of BAM file {}".format(bam_file_path)
#     assert sm_tag_list[0] == folder_name, "SM tag in BAM file {} do not correspond to folder name {}".format(bam_file_path, folder_name)


# check_bam_header(bam_file_path=bam_file_path)



print(df_h)
