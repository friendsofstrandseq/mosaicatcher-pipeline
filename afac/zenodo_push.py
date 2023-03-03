import requests
import os, sys
import json
import pandas as pd

ACCESS_TOKEN = sys.argv[1]
title, desc = sys.argv[2], sys.argv[3]
filelist_path = sys.argv[4]
compressed = bool(sys.argv[5])

# GET
# r = requests.get("https://sandbox.zenodo.org/api/deposit/depositions", params={"access_token": ACCESS_TOKEN})
r = requests.get("https://zenodo.org/api/deposit/depositions", params={"access_token": ACCESS_TOKEN})
print("GET")

# POST DEPOSITION
headers = {"Content-Type": "application/json"}
params = {"access_token": ACCESS_TOKEN}
r = requests.post(
    # "https://sandbox.zenodo.org/api/deposit/depositions",
    "https://zenodo.org/api/deposit/depositions",
    params=params,
    json={},
    # Headers are not necessary here since "requests" automatically
    # adds "Content-Type: application/json", because we're using
    # the "json=" keyword argument
    # headers=headers,
    headers=headers,
)
print("POST DEPOSITION")

deposition_id = r.json()["id"]
bucket_url = r.json()["links"]["bucket"]
print("Deposition id: {}".format(deposition_id))
print("Bucket URL: {}".format(bucket_url))

# The target URL is a combination of the bucket link with the desired filename
# seperated by a slash.

filelist = pd.read_csv(filelist_path, names=["file"])["file"].values.tolist()
print("Files list:")
print(filelist)
print("Putting files into zenodo bucket")


def zenodo_put_into_bucket(fp):
    """
    Function to put files into the zenodo bucket
    """
    if "/" not in fp:
        fp = "./" + fp
    print(fp)
    filename = os.path.basename(fp)
    read_mode = "rb" if compressed is True else "r"
    with open(fp, read_mode) as f_zenodo:
        r = requests.put(
            "%s/%s" % (bucket_url, filename),
            data=f_zenodo,
            params=params,
        )


# FOR LOOP TO ITERATE OF FILE LIST
for file in filelist:
    zenodo_put_into_bucket(file)


# METADATA
data = {
    "metadata": {
        "title": title,
        "upload_type": "other",
        "description": desc,
        # "communities": "friendsofstrandseq",
        "creators": [{"name": "Weber, Thomas", "affiliation": "EMBL Heidelberg"}],
    }
}
r = requests.put(
    # "https://sandbox.zenodo.org/api/deposit/depositions/%s" % deposition_id,
    "https://zenodo.org/api/deposit/depositions/%s" % deposition_id,
    params={"access_token": ACCESS_TOKEN},
    data=json.dumps(data),
    headers=headers,
)
print("Adding metadata, title: {}, desc: {}".format(title, desc))


r = requests.post(
    "https://zenodo.org/api/deposit/depositions/%s/actions/publish" % deposition_id,
    params={"access_token": ACCESS_TOKEN}
    # "https://sandbox.zenodo.org/api/deposit/depositions/%s/actions/publish" % deposition_id, params={"access_token": ACCESS_TOKEN}
)
print("FINAL POST")
