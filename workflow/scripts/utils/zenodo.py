import requests
import os, sys

ACCESS_TOKEN = sys.argv[1]
path = sys.argv[2]
filename = os.path.basename(path)

headers = {"Content-Type": "application/json"}
params = {'access_token': ACCESS_TOKEN}

r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions',
                   params=params,
                   json={},
                   # Headers are not necessary here since "requests" automatically
                   # adds "Content-Type: application/json", because we're using
                   # the "json=" keyword argument
                   # headers=headers,
                   headers=headers)

# bucket_url = r.json()["links"]["bucket"]


# The target URL is a combination of the bucket link with the desired filename
# seperated by a slash.
# with open(path, "rb") as fp:
#     r = requests.put(
#         "%s/%s" % (bucket_url, filename),
#         data=fp,
#         params=params,
#     )

# Old API
# Get the deposition id from the previous response
deposition_id = r.json()['id']
data = {'name': 'myfirstfile.csv'}
files = {'file': open('/path/to/myfirstfile.csv', 'rb')}
r = requests.post('https://zenodo.org/api/deposit/depositions/%s/files' % deposition_id,
                  params={'access_token': ACCESS_TOKEN}, data=data,
                  files=files)
r.status_code

r.json()

r = requests.post('https://zenodo.org/api/deposit/depositions/%s/actions/publish' % deposition_id, params={'access_token': ACCESS_TOKEN} )