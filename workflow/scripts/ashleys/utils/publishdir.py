import subprocess, os, sys

# input_list = [
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/cell_selection/labels_raw.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/cell_selection/labels_raw.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/cell_selection/labels_raw.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/cell_selection/labels.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/cell_selection/labels.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/cell_selection/labels.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/counts/IMR90E6E7PD103s1p2x01.info_raw",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/counts/IMR90E6E7PD106s1p3x01.info_raw",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/counts/TAllPDX6340p6RELs1p1x.info_raw",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/counts/IMR90E6E7PD103s1p2x01.txt.raw.gz",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/counts/IMR90E6E7PD106s1p3x01.txt.raw.gz",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/counts/TAllPDX6340p6RELs1p1x.txt.raw.gz",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/plots/counts/CountComplete.raw.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/plots/counts/CountComplete.normalised.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/plots/counts/CountComplete.raw.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/plots/counts/CountComplete.normalised.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/plots/counts/CountComplete.raw.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/plots/counts/CountComplete.normalised.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/plots/plate/ashleys_plate_predictions.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/plots/plate/ashleys_plate_probabilities.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/plots/plate/ashleys_plate_predictions.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/plots/plate/ashleys_plate_probabilities.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/plots/plate/ashleys_plate_predictions.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/plots/plate/ashleys_plate_probabilities.pdf",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/cell_selection/labels_positive_control_corrected.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/cell_selection/labels_positive_control_corrected.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/cell_selection/labels_positive_control_corrected.tsv",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD103s1p2x01/config/bypass_cell.txt",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/IMR90E6E7PD106s1p3x01/config/bypass_cell.txt",
#     "/scratch/tweber/DATA/MC_DATA/STOCKS_DEV/2023-06-23-HGFLGAFX7/TAllPDX6340p6RELs1p1x/config/bypass_cell.txt",
# ]
# publishdir = sys.argv[1]
# run = sys.argv[2]
# data_location = sys.argv[3]
data_location = snakemake.config["data_location"]
publishdir = snakemake.config["publishdir"]
run = snakemake.wildcards.folder.split("/")[-1]
print(snakemake.wildcards.folder)
print(run)
os.makedirs(f"{publishdir}/{run}", exist_ok=True)

# for file in input_list:
for file in list(snakemake.input.list_publishdir):
    print(file)
    # subprocess.Popen("echo {file}".format(file=file), shell=True, stdout=subprocess.PIPE)
    # print(snakemake.config["publishdir"])
    sub_path = "/".join(file.replace(data_location, "").split("/")[:-1])
    print(sub_path)
    # sub_path = "/".join(file.replace(snakemake.config["data_location"], "").split("/")[:-2])
    sub_path = sub_path if "/" in sub_path else "/{}".format(sub_path)
    print(sub_path)
    # # folder_path = snakemake.config["publishdir"] + sub_path + "/"
    folder_path = f"{publishdir}/{run}{sub_path}"
    # # folder_path = snakemake.config["publishdir"] + "/".join(file.replace(snakemake.config["data_location"], ""))
    print(folder_path)
    os.makedirs(folder_path, exist_ok=True)

    # # subprocess.Popen("mkdir -p {folder_path}".format(folder_path=folder_path), shell=True, stdout=subprocess.PIPE)
    print("rsync --ignore-existing -avzh --progress {file} {folder_path}".format(file=file, folder_path=folder_path))
    subprocess.Popen(
        "rsync --ignore-existing -avzh --progress {file} {folder_path}".format(file=file, folder_path=folder_path),
        shell=True,
        stdout=subprocess.PIPE,
    )
