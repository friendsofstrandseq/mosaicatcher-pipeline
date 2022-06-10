import os

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))

import os, sys
from pprint import pprint 
from tqdm import tqdm 
import shutil 

scratch_dir = "/scratch/jeong/pipeline_mosaicatcher_backup13/pipeline_20190625_HGSVC24_rerun/"

map_dict = {
	"pipeline_20190625_HGSVC24_rerun_GM12329" : "/g/korbel2/StrandSeq/20200109_U24/20191023_GM12329A/bam/",
	"pipeline_20190625_HGSVC24_rerun_GM18534" : "/g/korbel2/StrandSeq/20200109_U24/20191120_GM18534B/bam/",
	"pipeline_20190625_HGSVC24_rerun_GM18939" : "/g/korbel2/StrandSeq/20200109_U24/20191122_GM18939/bam/",
	"pipeline_20190625_HGSVC24_rerun_GM19650" : "/g/korbel2/StrandSeq/20200109_U24/20191120_GM19650A/bam/",
	"pipeline_20190625_HGSVC24_rerun_GM19983" : "/g/korbel2/StrandSeq/20200109_U24/20191024_GM19983/bam/",
	"pipeline_20190625_HGSVC24_rerun_GM20847" : "/g/korbel2/StrandSeq/20200109_U24/20191203_GM20847B/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG00096_HG00171" : "/g/korbel2/StrandSeq/20200109_U24/20191203_HG00096/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG00864" : "/g/korbel2/StrandSeq/20200109_U24/20191122_HG00864/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG01114" : "/g/korbel2/StrandSeq/20200109_U24/20190911_HG01114/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG01505" : "/g/korbel2/StrandSeq/20200109_U24/20191120_HG01505/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG01596" : "/g/korbel2/StrandSeq/20200109_U24/20191122_HG01596/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG02011" : "/g/korbel2/StrandSeq/20200109_U24/20191023_HG02011/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG02492" : "/g/korbel2/StrandSeq/20200109_U24/20200630_HG02492/strand-seq-pipeline/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG02587" : "/g/korbel2/StrandSeq/20200109_U24/20191011_HG02587/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG02818" : "/g/korbel2/StrandSeq//20200421_RO1_Samples/HG02818/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03009" : "/g/korbel2/StrandSeq/20200109_U24/20191105_HG03009/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03065" : "/g/korbel2/StrandSeq/20200109_U24/20191105_HG03065/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03125" : "/g/korbel2/StrandSeq/20200421_RO1_Samples/20200630_HG03125/strand-seq-pipeline/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03371" : "/g/korbel2/StrandSeq/20200109_U24/20191024_HG03371/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03486" : "/g/korbel2/StrandSeq/20200421_RO1_Samples/20200622_HG03486/strand-seq-pipeline/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03683" : "/g/korbel2/StrandSeq/20200109_U24/20191105_HG03683/bam/",
	"pipeline_20190625_HGSVC24_rerun_HG03732" : "/g/korbel2/StrandSeq/20200109_U24/20190917_HG03732/bam/"
}

final_dir = "/g/korbel2/weber/MosaiCatcher_files/HGSVC_correct/"


samples = {d : os.listdir(scratch_dir + d + "/bam/")[-1] for d in os.listdir(scratch_dir)}

pprint(samples)

if os.path.exists(final_dir):
    shutil.rmtree(final_dir)

[os.makedirs(final_dir + s + "/all/", exist_ok=True) for p,s in samples.items()]
[os.makedirs(final_dir + s + "/selected/", exist_ok=True) for p,s in samples.items()]

for p,s in tqdm(samples.items()):
    for f in os.listdir(map_dict[p]):
        if os.path.islink(final_dir + s + "/all/" + f) is False:
            os.symlink(map_dict[p] + f, final_dir + s + "/all/" + f)


for d in tqdm(os.listdir(scratch_dir)):
    for f in os.listdir(scratch_dir + d + "/bam/" + samples[d] + "/selected/"):
        # print(final_dir + samples[d] + "/selected/" + f)
        os.symlink(final_dir + samples[d] + "/all/" + f, final_dir + samples[d] + "/selected/" + f)


# list_files(final_dir)