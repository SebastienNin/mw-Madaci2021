#!/bin/sh
# Aim: a simple script to find all files in Salvatore Spicuglia's tree that 
# do not correspond to sample_name defined in Sequencing_summary.xlsx.
# This is useful to clean Sacapus and NAS directories of files corresponding
# to a modified sample_name in Sequencing_summary.xlsx.
# Currently, copy/paste sample_name for Sequencing_summary.xlsx to replace
# clean_illumina_sequencing_data_of_old_sample_names.sample_name content
# Modified from: 
# https://unix.stackexchange.com/questions/301717/how-to-find-all-symbolic-links-pointing-to-any-file-directory-inside-a-given-dir
set --
while IFS= read -r pattern
do
	set -- "$@" -o -name "$pattern"'*'
done < /gpfs/projects/spicuglia/mw-sst/src/bash/clean_illumina_sequencing_data_of_old_sample_names.sample_name

# bash magic to remove the first '-o'
if [ $# -ne 0 ]; then
	shift
	set -- "$@"
fi

echo "$@"
# This script is dangerous.
# Uncomment delete only if your are happy with selected files.
#find -L /gpfs/projects/spicuglia/mw-sst/out/ln/alias/sst -not \( "$@" \) -and \( -type l -o -type f \) #-delete

# additional files not listed (html, txt, json, tsv) are from multiQC as they do not contain the sampe_name
find -L /gpfs/projects/spicuglia/mw-sst/out/ln/alias/sst -not \( "$@" -o -name '*.html' -o -name '*.txt' -o -name '*json' -o -name '*.log' -o -name '*.tsv' \) -and \( -type l -o -type f \) #-delete
#find -L /mnt/SynologySpicuglia/guillaume/illumina_sequencing_data -not \( "$@" -o -name '*.html' -o -name '*.txt' -o -name '*json' -o -name '*.log' -o -name '*.tsv' \) -and \( -type l -o -type f \) -delete



