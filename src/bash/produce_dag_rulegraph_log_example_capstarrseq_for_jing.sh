#!/usr/bin/env bash

# 1. Set variables to make the script works from any directory
BASH_DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
SM_DIR=$BASH_DIR/../..

# 2. Create required conda environnements
eval "$(conda shell.bash hook)"
conda env create -f $SM_DIR/../mw-lib/src/snakemake/envs/snakemake.yaml
conda activate snakemake
cd $SM_DIR

# 3. Download Sequencing Summary containing the references for samples
wget -N https://dl.dropboxusercontent.com/s/g4hx58039uyou2k/Sequencing_summary.xlsx

# 4. Define Snakemake target
SMI="out/capstarrseq/merge_all_data_hProm/awk/extend_reads_314/bedtools/bamtobed/samtools/sort_-n/ln/alias/sst/all_samples/hg19/bam/Dao2017_K562_CapStarr_rep1_over_Dao2017_CapStarr_input.allData.tsv"

# 5. Produce Snakemake target (Only work in Sacapus)
snakemake -prk --rerun-incomplete --cores 7 --cluster "qsub -V -q tagc -o log/qsub -e log/qsub -l nodes=1:ppn={threads}" --use-conda $SMI

# 6. Produce Rulegraph, DAG and lists of all tasks
mkdir -p $SM_DIR/out/snakemake
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI --dag > $SM_DIR/out/snakemake/example_capstarrseq_jing_dag.dot
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI --rulegraph > $SM_DIR/out/snakemake/example_capstarrseq_jing_rulegraph.dot
snakemake -prk --rerun-incomplete --cores 2 --use-conda $SMI -fn > $SM_DIR/out/snakemake/example_capstarrseq_jing_tasks.txt

cd $SM_DIR/out/snakemake
dot -Tpdf -o example_capstarrseq_jing_dag.pdf example_capstarrseq_jing_dag.dot
dot -Tpdf -o example_capstarrseq_jing_rulegraph.pdf example_capstarrseq_jing_rulegraph.dot
