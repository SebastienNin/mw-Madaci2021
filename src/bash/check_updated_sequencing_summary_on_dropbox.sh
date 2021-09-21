#!/bin/bash

remote_file="dl.dropboxusercontent.com/s/dae6otw5honour0/Sequencing_summary.xlsx"
local_file="mw/out/wget/https/dl.dropboxusercontent.com/s/dae6otw5honour0/Sequencing_summary.xlsx"

modified=$(curl --silent --head $remote_file | \
	             awk '/^Last-Modified/{print $0}' | \
				              sed 's/^Last-Modified: //')
remote_ctime=$(date --date="$modified" +%s)
local_ctime=$(stat -c %z "$local_file")
local_ctime=$(date --date="$local_ctime" +%s)
