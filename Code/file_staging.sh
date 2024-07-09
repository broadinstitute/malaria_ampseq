#!/bin/bash


###################################################################
##Copy files to the working directory and run the AmpSeq pipeline##
###################################################################

mkdir fq_dir
mkdir references
gsutil -m cp -r ~{sep = ' ' path_to_r1} fq_dir/
gsutil -m cp -r ~{sep = ' ' path_to_r2} fq_dir/
gstuil -m cp -r ~{sep = ' ' path_to_flist} references/
gstuil -m cp -r ~{sep = ' ' pr1} references/
gstuil -m cp -r ~{sep = ' ' pr2} references/
gstuil -m cp -r ~{sep = ' ' reference} references/

if [[ "~{reference2}" != '' ]]; then
	gsutil -m cp -r ~{sep = ' ' reference2} references/
fi

if [[ "~{path_to_snv}" != '' ]]; then
	gsutil -m cp -r ~{sep = ' ' path_to_snv} references/
fi
