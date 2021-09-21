#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..
RMD_DIR=$BASH_DIR/../rmd/

# 2. Create required conda environnements
conda env create -f $SM_DIR/../mw-lib/src/snakemake/envs/snakemake.yaml
#conda env create -f $SM_DIR/src/snakemake/envs/r_rmd_sarcoma.yaml

# 3. Find Snakemake targets in the Rmd files
TARGET=`grep '^out/cellranger <- "' $RMD_DIR/*.Rmd` 
echo $SMI

# 4. Run Snakemake on all targets in the Rmd files
eval "$(conda shell.bash hook)"
conda activate snakemake
cd $SM_DIR
snakemake -prk --rerun-incomplete --cores 16 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q lifescope -o log/qsub -e log/qsub -l nodes=1:ppn={threads},walltime=99:00:00" --use-conda --max-jobs-per-second 3 $TARGETS
