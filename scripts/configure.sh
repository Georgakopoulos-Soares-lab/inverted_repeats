#!/bin/bash


echo "Installing snakemake. Please wait..."
mamba install -c bioconda -c biopython snakemake==7.32.4

echo "Installing biopython..."
mamba install -c bioconda biopython

echo "Installing pipeline requirements..."
pip install -r requirements.txt

