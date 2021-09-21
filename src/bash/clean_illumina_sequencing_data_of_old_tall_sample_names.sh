#!/bin/sh
# Aim: a simple script to find all files in Salvatore Spicuglia's tree that 
# do not correspond to sample_name defined in Sequencing_summary.xlsx.
# This is useful to clean Sacapus and NAS directories of files corresponding
# to a modified sample_name in Sequencing_summary.xlsx.
# Currently, copy/paste sample_name for Sequencing_summary.xlsx to replace
# clean_illumina_sequencing_data_of_old_sample_names.sample_name content
# Uncomment delete only if your are happy with selected files.
#find -L /gpfs/projects/spicuglia/mw/out/ln/alias/sst -not \( "$@" \) -and \( -type l -o -type f \) -delete

find -L /gpfs/projects/spicuglia/mw-sst/out/ln/alias/sst -name "T*_TALL_*" -and \( -type l -o -type f \) -delete
find -L /gpfs/projects/spicuglia/mw-sst/out/ln/alias/sst -name "T4_*" -and \( -type l -o -type f \) -delete
