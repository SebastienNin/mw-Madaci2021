#!/bin/bash

# Activate a controlled environment for all users.
#. /gpfs/projects/spicuglia/mw/opt/miniconda/etc/profile.d/conda.sh
#conda activate /gpfs/tagc/home/gcharbonnier/mw/opt/miniconda/envs/dev

# Download the configuration table.
cd /gpfs/projects/spicuglia/mw-sst
wget -N https://dl.dropboxusercontent.com/s/g4hx58039uyou2k/Sequencing_summary.xlsx

# Run continuous integration tests. (obsolete now with github actions)
# Send an error by mail if one file can not be generated from inputs files and current workflow state.
#cd ../mw-sst-rulegraph
#snakemake --use-conda -f --rerun-incomplete -prk --verbose out/snakemake/stdout_--rulegraph_ci-sst.dot

# Run analysis
cd ../mw-sst
snakemake -prk --rerun-incomplete --cores 7 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q tagc -o log/qsub -e log/qsub -l nodes=1:ppn={threads}" --use-conda config_targets

# Run multiQC on available log files
# Currently this is just a test to check differences between directory of execution:
snakemake -prk --rerun-incomplete --cores 7 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q tagc -o log/qsub -e log/qsub -l nodes=1:ppn={threads}" --use-conda config_qc_targets

# Create output directories for the edge case where workflow is first executed with no target (maybe not required anymore, to test)
#mkdir -p out/ln/alias/sst/ out/config_targets/out/ln/alias/sst/

# Create a dereferenced tree with multiple hardlinks for the same files.
# This is required to get rsync sync only one time files that are duplicated in multiple trees.
find out/ln/alias/sst/ -type l -print0 | while read -r -d $'\0' file; do mkdir -p `dirname out/config_targets/$file` ; ln -Lf $file out/config_targets/$file ; done
#find out/ln/alias/sst/ -type l -print0 | while IFS= read -r -d $'\\0' file; do mkdir -p `dirname out/config_targets/$file` ; ln -Lf $file out/config_targets/$file ; done

# Rsync can now send changed files to the NAS
DESTINATION=/mnt/SynologySpicuglia/guillaume/illumina_sequencing_data
# Check first that the NAS is currently mounted (it is often unmounted after cluster reboot)
if [ -d $DESTINATION ]
then
	rsync --recursive --times -H -vPh --stats out/config_targets/out/ln/alias/sst/ $DESTINATION
	rsync --recursive --times -H -vPh --stats out/multiqc/dir_--interactive/ln/alias/sst/ $DESTINATION/multiQC/interactive
	rsync --recursive --times -H -vPh --stats out/multiqc/dir_--flat/ln/alias/sst/ $DESTINATION/multiQC/flat
	rsync --recursive --times -H -vPh --stats out/find/md5sum/ln/alias/sst/ $DESTINATION
else
	echo "NAS unreachable, ask SI to mount it again."
fi

# Huge processed files may now be removed from the output directory.
# Small ones should be conserved in order to have multiqc find the log files of previous runs,
# thus allowing to still have "done" samples in the QC report.
#find out -type f -size +100000k -delete

