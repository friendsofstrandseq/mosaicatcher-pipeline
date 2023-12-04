#!/bin/bash

# Check if a Dockerfile path is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path-to-Dockerfile>"
    exit 1
fi

DOCKERFILE=$1

# Check if the Dockerfile exists
if [ ! -f "$DOCKERFILE" ]; then
    echo "Dockerfile not found: $DOCKERFILE"
    exit 1
fi

# Extract the R environment variable
Renv=$(grep -P "\/rtools.*environment\.yaml" "$DOCKERFILE" | sed "s/\//\t/g" | cut -f 5)

# Check if Renv is extracted
if [ -z "$Renv" ]; then
    echo "R environment variable not found in the Dockerfile."
    exit 1
fi

# Append custom steps to the Dockerfile
{
    echo '\n'
    echo "# CUSTOM PART"
    echo "RUN wget https://zenodo.org/record/7697400/files/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz -P /workflow/data/ref_genomes/"
    echo "COPY /workflow/scripts/utils/install_R_package.R /conda-envs/"
    echo "RUN chmod -R 0777 /conda-envs/$Renv/lib/R/library && /conda-envs/$Renv/bin/Rscript /conda-envs/install_R_package.R /workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz"
} >>"$DOCKERFILE"

echo "Custom steps added to $DOCKERFILE"
