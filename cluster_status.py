#!/usr/bin/env python

# Code by Jelle Scholtalbers 
import subprocess
import sys
import re

jobid = sys.argv[1]

try:
    output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())
except subprocess.CalledProcessError:
    print("failed")
    sys.exit(0)

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
    try:
        output = str(subprocess.check_output("grep 'slurmstepd: error: Exceeded step memory limit at some point.' slurm-%s.out" % jobid, shell=True))
    except subprocess.CalledProcessError:
        # grep fails to find error (or fails to find log file): success
        print("success")
    else:
        print("failed")
    sys.exit(0)
elif any(s in output for s in running_status):
    print("running")
else:
    print("failed")

