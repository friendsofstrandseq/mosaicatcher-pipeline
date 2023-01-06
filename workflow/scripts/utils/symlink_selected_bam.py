import subprocess

if snakemake.config["use_light_data"] is False:
    subprocess.Popen(
        "ln -s  {input_bam} {output_bam}".format(input_bam=snakemake.input.bam, output_bam=snakemake.output.bam),
        shell=True,
        stdout=subprocess.PIPE,
    )
    subprocess.Popen(
        "ln -s  {input_bai} {output_bai}".format(input_bai=snakemake.input.bai, output_bai=snakemake.output.bai),
        shell=True,
        stdout=subprocess.PIPE,
    )
else:
    subprocess.Popen(
        "cp {input_bam} {output_bam}".format(input_bam=snakemake.input.bam, output_bam=snakemake.output.bam),
        shell=True,
        stdout=subprocess.PIPE,
    )
    subprocess.Popen(
        "cp {input_bai} {output_bai}".format(input_bai=snakemake.input.bai, output_bai=snakemake.output.bai),
        shell=True,
        stdout=subprocess.PIPE,
    )
