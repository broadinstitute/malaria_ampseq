#!/usr/bin/env python

import os
import sys
import argparse
import json
import subprocess
import csv
import gzip
import pandas as pd

from Bio import SeqIO

import amplicon_decontamination as ad
import asv_to_cigar as ac

def main():
	"""
	Implementation of the amplicon decontamination pipeline for use and 
	distribution in TERRA v=1.0

	Usage: python Code/Amplicon_TerraPipeline.py --config config.json ...
	
	Returns:
	None. However, see the functions documentation.
	"""

	### LOAD ARGUMENTS

	#Parse command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help="Path to config.json file.", required =True)
	parser.add_argument('--terra', action="store_true", help="Specify whether the pipeline is being run in the Terra platform.")
	parser.add_argument('--meta', action="store_true", help="Specify if metadata file must be created. This flag runs the pipeline from the beginning.")
	parser.add_argument('--repo', action="store_true", help="Specify if the reports must be created.")
	parser.add_argument('--contamination', action="store_true", help="Specify if contamination detection is performed.")
	parser.add_argument('--adaptor_removal', action="store_true", help="Specify if adaptor removal needed.")
	parser.add_argument('--separate_reads', action="store_true", help="Specify if reads must be separated by size.")
	parser.add_argument('--primer_removal', action="store_true", help="Specify if primer removal needed.")
	parser.add_argument('--dada2', action="store_true", help="Specifiy if standard preprocess merge with DADA2 is performed.")
	parser.add_argument('--postproc_dada2', action="store_true", help="Specifiy if postProcess of DADA2 results is perfomed.")
	parser.add_argument('--asv_to_cigar', action="store_true", help="Specifiy if the ASV to CIGAR transformation is perfomed.")

	args = parser.parse_args()

	#Check minimum arguments and contracdicting flags
	if args.terra:
		print("Pipeline is running in Terra. Adjusted paths will be used.")
	else:
		print("Pipeline not running in Terra. Default paths will be used.")
				
	#Configuration aguments will be parsed from config.json
	with open(args.config, 'r') as config_file:
		config_inputs = json.load(config_file)
		path_to_fq = config_inputs['path_to_fq']
		path_to_flist = 'barcodes_matches.csv'
		pr1 = 'primers_fw.fasta'
		pr2 = 'primers_rv.fasta'
		path_to_snv = 'snv_filters.txt'
		if 'pattern_fw' in config_inputs.keys(): pattern_fw = config_inputs['pattern_fw']
		if 'pattern_rv' in config_inputs.keys(): pattern_rv = config_inputs['pattern_rv']
		if 'Class' in config_inputs.keys(): Class = config_inputs['Class']
		if 'maxEE' in config_inputs.keys(): maxEE = config_inputs['maxEE']
		if 'trimRight' in config_inputs.keys(): trimRight = config_inputs['trimRight']
		if 'minLen' in config_inputs.keys(): minLen = config_inputs['minLen']
		if 'truncQ' in config_inputs.keys(): truncQ = config_inputs['truncQ']
		if 'matchIDs' in config_inputs.keys(): matchIDs = config_inputs['matchIDs']
		if 'max_consist' in config_inputs.keys(): max_consist = config_inputs['max_consist']
		if 'omegaA' in config_inputs.keys(): omegaA = config_inputs['omegaA']
		if 'saveRdata' in config_inputs.keys(): saveRdata = config_inputs['saveRdata']
		if 'justConcatenate' in config_inputs.keys(): justConcatenate = config_inputs['justConcatenate']
		if 'maxMismatch' in config_inputs.keys(): maxMismatch = config_inputs['maxMismatch']
		if 'adjust_mode' in config_inputs.keys(): adjust_mode = config_inputs['adjust_mode']
		if 'no_ref' in config_inputs.keys(): no_ref = config_inputs['no_ref']
		if 'strain' in config_inputs.keys(): strain = config_inputs['strain']
		if 'strain2' in config_inputs.keys(): strain2 = config_inputs['strain2']
		if 'polyN' in config_inputs.keys(): polyN = int(config_inputs['polyN'])
		if 'min_reads' in config_inputs.keys(): min_reads = int(config_inputs['min_reads'])
		if 'min_samples' in config_inputs.keys(): min_samples = int(config_inputs['min_samples'])
		if 'max_snv_dist' in config_inputs.keys(): max_snv_dist = int(config_inputs['max_snv_dist'])
		if 'max_indel_dist' in config_inputs.keys(): max_indel_dist = int(config_inputs['max_indel_dist'])
		if 'include_failed' in config_inputs.keys(): include_failed = eval(config_inputs['include_failed'])
		if 'exclude_bimeras' in config_inputs.keys(): exclude_bimeras = eval(config_inputs['exclude_bimeras'])
		if 'verbose' in config_inputs.keys(): verbose = eval(config_inputs['verbose'])
		if 'adapter' in config_inputs.keys(): adapter = eval(config_inputs['adapter'])

	### PREPARE OUTPUT DIRECTORIES

	#Generate the path to the Results directory.
	global res_dir
	global rep_dir
	res_dir = os.path.abspath(os.path.join("Results"))
	rep_dir = os.path.abspath(os.path.join("Report"))

	#Restart Results and Report directory if the workflow runs form the start
	if args.meta:
		os.system("rm -rf " + res_dir)
		os.mkdir(res_dir)

	if args.repo:
		os.system("rm -rf " + rep_dir)
		os.mkdir(rep_dir)
	
	#Create metadata files
	if args.meta:
		print("Creating meta files")
		ad.flush_dir(res_dir, "Fq_metadata")
		ad.create_meta(path_to_fq, res_dir, "Fq_metadata", "rawfilelist.tsv", pattern_fw, pattern_rv)

		#List missing file from the complete list of files
		with open(path_to_flist, 'r') as input_file:
			next(input_file)
			samples = [line.split(',')[0] for line in input_file]

		with open('Results/Fq_metadata/rawfilelist.tsv', 'r') as raw_files:
			raw_file_samples = [line.split('\t')[0] for line in raw_files]

		missing_samples = [sample.strip() for sample in samples if sample not in raw_file_samples]

		with open('Results/missing_files.tsv', 'w') as output_file:
			output_file.write('\n'.join(missing_samples))

	### EXECUTE

	#Remove adaptors
	#Most sequences must be adaptor free; just in case, run this step to eliminate any lingering adaptors.
	if args.adaptor_removal:
		print("Removing adaptors")
		ad.flush_dir(res_dir, "AdaptorRem")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples:
			slist = sample.split()
			ad.adaptor_rem(slist[0], slist[1], slist[2], res_dir, "AdaptorRem")
	
		ad.create_meta(os.path.join(res_dir, "AdaptorRem"), res_dir, "AdaptorRem", "adaptorrem_meta.tsv",
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")

	#Merge forward and reverse reads with bbmerge. Only for reads that overlap and
	#Process merge report and generate RStudio plots
	if args.contamination:
		print("Concatenating reads for contamination detection")
		ad.flush_dir(res_dir, "Merge")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv") , 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ad.mergereads(slist[0], slist[1], slist[2], res_dir, "Merge")

		ad.flush_dir(rep_dir, "Merge")
		meta = open(os.path.join(res_dir, "Merge", "merge_meta.tsv"), 'r')
		samples = meta.readlines()
		
		print("Extracting fields from bbmerge reports")
		for sample in samples:
			slist = sample.split()
			ad.extract_bbmergefields(slist[0], slist[1], slist[3], path_to_flist, res_dir, rep_dir, "Merge", args.terra)

	#Determine if files must be cleaned from contamination. Note that running this protocol will regenerate, 
	#Fq_Metadata and AdaptorRem.
	if args.contamination and args.dada2:
		print("Making files with no contaminating reads")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		ad.flush_dir(res_dir, "Clean_Data_Repo")
		for sample in samples:
			slist = sample.split()
			tsv_name = slist[0]+"_final.tsv"
			path_to_match_report = os.path.join(rep_dir, "Merge", tsv_name)
			if os.path.isfile(path_to_match_report):	
				print(path_to_match_report)
				path_outfile_r1 = os.path.join(res_dir, "Clean_Data_Repo", slist[0]+'_L001_R1_001.fastq')
				path_outfile_r2 = os.path.join(res_dir, "Clean_Data_Repo", slist[0]+'_L001_R2_001.fastq')
				ad.filter_fastq_by_read_names(slist[1], slist[2], path_to_match_report, path_outfile_r1, path_outfile_r2)
			else:	
				print("File is empty")
				print(path_to_match_report)

		ad.flush_dir(res_dir, "Fq_metadata")
		ad.create_meta(os.path.join(res_dir, "Clean_Data_Repo"), res_dir, "Fq_metadata", "rawfilelist.tsv", pattern_fw, pattern_rv)

		print("Removing adaptors of clean read set")
		ad.flush_dir(res_dir, "AdaptorRem")
		meta = open(os.path.join(res_dir, "Fq_metadata", "rawfilelist.tsv"), 'r')
		samples = meta.readlines()
		
		for sample in samples:
			slist = sample.split()
			ad.adaptor_rem(slist[0], slist[1], slist[2], res_dir, "AdaptorRem")
	
		ad.create_meta(os.path.join(res_dir, "AdaptorRem"), res_dir, "AdaptorRem", "adaptorrem_meta.tsv",
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")

	if args.separate_reads:
		print("Entering demultiplexing algorithm")
		meta = open(os.path.join(res_dir, "AdaptorRem", "adaptorrem_meta.tsv"), 'r')
		samples = meta.readlines()
		ad.flush_dir(res_dir, "Demultiplex_by_Size")

		#Get the reads size
		print("Getting the read size")
		lengths_fw = []
		lengths_rv = []
		for sample in samples:
			longest_sequence_length_fw = ad.find_longest_sequence_length(sample.split()[1])
			longest_sequence_length_rv= ad.find_longest_sequence_length(sample.split()[2])
			lengths_fw.append(longest_sequence_length_fw)
			lengths_rv.append(longest_sequence_length_rv)

		percentile = lambda lengths: sorted(lengths)[int(len(lengths) * 0.95):]
		p_fw = percentile(lengths_fw)
		p_rv = percentile(lengths_rv)

		all_equal_fw = len(set(p_fw)) == 1
		all_equal_rv = len(set(p_rv)) == 1
		if all_equal_fw and all_equal_rv:
			print("Are all values in the top 95% percentile of read size equal?", all_equal_fw and all_equal_rv)
			print("The following values will be used as the standard size of the reads in this run:")
		else:
			print("Are all values in the top 95% percentile of read size equal?", all_equal_fw and all_equal_rv)
			print("Largest values found for the forward and reverse read will be used. These values are:")
		read_size_fw = sorted(lengths_fw)[-1]
		read_size_rv = sorted(lengths_rv)[-1]
		print("Forward read:", read_size_fw)
		print("Reverse read:", read_size_rv)
		print("If these sizes, do not match the expected size of your technology, consider rerurring the pipeline after manually providing the size of your reads")
		#Get and remove the adapter sequence
		print("Removing the adapter sequence from the primers")
		if adapter is None:
			adapter_fw = ad.find_common_subsequence(pr1)
			adapter_rv = ad.find_common_subsequence(pr2)
		else:
			adapter_fw = adapter
			adapter_rv = adapter
		print("Adapter forward:", adapter_fw)
		print("Adapter reverse:", adapter_rv)
		ad.remove_adapter(pr1, adapter_fw, 'primer_fw_no_adapter.fasta')
		ad.remove_adapter(pr2, adapter_rv, 'primer_rv_no_adapter.fasta')

		#Get the size of the reference ASVs
		print("Getting the size of the reference ASV")
		if os.path.exists("reference_panel_1.fasta"):
			ref_files = ["reference_panel_1.fasta"]
			if os.path.exists("reference_panel_2.fasta"):
				ref_files.append("reference_panel_2.fasta")
		else:
			sys.exit('At the least one reference is necessary to separate reads according to ASV target size.')
		asv_lengths = {}
		for ref_file in ref_files:
			for record in SeqIO.parse(ref_file, "fasta"):
				if record.id not in asv_lengths:
					asv_lengths[record.id] = len(record.seq)
				elif asv_lengths[record.id] < len(record.seq):
					#When two reference genomes are provided, if an ASV is already present, always reference the longest version.
					asv_lengths[record.id] = len(record.seq)

		with open(os.path.join(res_dir, "asv_lengths.tsv"), 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['ASV', 'Length'])
			for seq_id, length in asv_lengths.items():
				writer.writerow([seq_id, length])

		print("Demultiplexing reads by size of reads according to their target amplicon")
		#Make dictionary of wells
		sample_number = [f"S{i}" for i in range(1, 193)]
		illumina_well = [f"{char}{num}" for char in list("ABCDEFGH") for num in range(1, 13)]*2
		sample_dict = dict(zip(sample_number, illumina_well))

		for sample in samples:
			slist = sample.split()
			ad.demultiplex_per_size(slist[0], slist[1], slist[2], 'primer_fw_no_adapter.fasta', 'primer_rv_no_adapter.fasta', res_dir, "Demultiplex_by_Size", read_size_fw, read_size_rv, asv_lengths, args.contamination, sample_dict)
		
		#Create Metafile for reads with no overlap
		ad.create_meta(os.path.join(res_dir, "Demultiplex_by_Size"), res_dir, "Demultiplex_by_Size", "demux_nop_meta.tsv",
			pattern_fw="*_nop_L001_R1_001.fastq.gz", pattern_rv="*_nop_L001_R2_001.fastq.gz")
		ad.create_meta(os.path.join(res_dir, "Demultiplex_by_Size"), res_dir, "Demultiplex_by_Size", "demux_op_meta.tsv",
			pattern_fw="*_op_L001_R1_001.fastq.gz", pattern_rv="*_op_L001_R2_001.fastq.gz")

	#Remove primers
	#For a set where all reads have overlap
	if args.primer_removal:
		print("Removing primers")
		#Extract primer for the target without amplicons
		# Read input fasta file
		if args.contamination:
			fw = 'primer_fw_no_adapter.fasta'
			rv = 'primer_rv_no_adapter.fasta'
		else:
			fw = 'primers_fw.fasta'
			rv = 'primers_rv.fasta'
			
		records_fw = SeqIO.parse(fw, 'fasta')
		common_subsequences = ad.find_common_subsequences(records_fw)
		ad.write_common_subsequences_to_fasta(common_subsequences, 'amp_primer_fw.fasta')
		records_rv = SeqIO.parse(rv, 'fasta')
		common_subsequences = ad.find_common_subsequences(records_rv)
		ad.write_common_subsequences_to_fasta(common_subsequences, 'amp_primer_rv.fasta')

		ad.flush_dir(res_dir, "PrimerRem")

		#Trim primers off non-overlapping targets
		meta = open(os.path.join(res_dir, "Demultiplex_by_Size", "demux_nop_meta.tsv"), 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ad.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", "amp_primer_fw.fasta", "amp_primer_rv.fasta", "mixed_nop")

		#Metafile for trimmed non-op target reads
		ad.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv", 
			pattern_fw="*_mixed_nop_1.fq.gz", pattern_rv="*_mixed_nop_2.fq.gz")

		#Trim primers off overlapping targets
		meta = open(os.path.join(res_dir, "Demultiplex_by_Size", "demux_op_meta.tsv"), 'r')
		samples = meta.readlines()
		for sample in samples:
			slist = sample.split()
			ad.trim_primer(slist[0], slist[1], slist[2], res_dir, "PrimerRem", "amp_primer_fw.fasta", "amp_primer_rv.fasta", "mixed_op")

		#Metafile for trimmed overlapping target reads
		ad.create_meta(os.path.join(res_dir, "PrimerRem"), res_dir, "PrimerRem", "mixed_op_prim_meta.tsv",
			pattern_fw="*_mixed_op_1.fq.gz", pattern_rv="*_mixed_op_2.fq.gz")

#	#For a set that mixes reads with and without overlap
	if args.dada2:
		print("Running DADA2")
		#Run DADA2 on op targets
		ad.flush_dir(res_dir, "DADA2_OP", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_op_prim_meta.tsv")
		justConcatenate=0	
		ad.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch,saveRdata, res_dir, "DADA2_OP", args.terra)
		seqtab_op = os.path.join(res_dir, 'DADA2_OP', 'seqtab.tsv')
		bimera_op = os.path.join(res_dir, 'DADA2_OP', 'ASVBimeras.txt')

		#Run DADA2 on non-op targets
		ad.flush_dir(res_dir, "DADA2_NOP", "QProfile")
		path_to_meta = os.path.join(res_dir, "PrimerRem", "mixed_nop_prim_meta.tsv")
		justConcatenate=1	
		ad.run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch,saveRdata, res_dir, "DADA2_NOP", args.terra)
		seqtab_nop = os.path.join(res_dir, 'DADA2_NOP', 'seqtab.tsv')
		bimera_nop = os.path.join(res_dir, 'DADA2_NOP', 'ASVBimeras.txt')

		#ASV modification block for non-op targets and merge two ASV tables
		if os.path.exists("reference_panel_1.fasta"):
			if args.terra:
				path_to_program = os.path.join("/", "Code/adjustASV.R")
			else:
				path_to_program = os.path.join("Code/adjustASV.R")
			adjASV = ['Rscript', path_to_program, '-s', seqtab_nop, '-ref', "reference_panel_1.fasta",
			'-dist', adjust_mode,
			'-o', os.path.join(res_dir, 'DADA2_NOP', 'correctedASV.txt')]
			print(adjASV)
			procASV = subprocess.Popen(adjASV)
			procASV.wait()
			seqtab_corrected = os.path.join(res_dir, 'DADA2_NOP', 'seqtab_corrected.tsv')
			seqtab = ad.merge_seqtab(seqtab_op, seqtab_corrected)
			bimera = ad.merge_bimeras(bimera_op, bimera_nop)
		else:
			print('--reference file not found. skipping ASV correction..')
			seqtab = ad.merge_seqtab(seqtab_op, seqtab_nop)
			bimera = ad.merge_bimeras(bimera_op, bimera_nop)

		seqtab.to_csv(os.path.join(res_dir, 'seqtab.tsv'), sep = "\t", index=False)
		bimera.to_csv(os.path.join(res_dir, 'ASVBimeras.txt'), sep = "\t", index=False)

	if args.postproc_dada2:	
		print("Performing PostProc")	
		ad.flush_dir(res_dir, "PostProc_DADA2")
		
		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')

		if args.terra:
			path_to_program = os.path.join("/", "Code/postProc_dada2.R")
		else:
			path_to_program = os.path.join("Code/postProc_dada2.R")

		postProc = ['Rscript', path_to_program, 
				'-s', path_to_seqtab, 
				'-b', os.path.join(res_dir, 'ASVBimeras.txt'),
				'-snv', os.path.join(path_to_snv),
				'--indel_filter', '0.895',
				'-o', os.path.join(res_dir, 'PostProc_DADA2', 'ASVTable.txt'),
				'--fasta']

		if no_ref == 'True':
			postProc.extend(['-no_ref'])
		else:
			postProc.extend(['--reference', "reference_panel_1.fasta", '--strain', strain])
			if os.path.exists("reference_panel_2.fasta"):
				postProc.extend(['--reference2', "reference_panel_2.fasta", '--strain2', strain2])

		print(postProc)
		procASV = subprocess.Popen(postProc)
		procASV.wait()

	#ASV to CIGAR
	#Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.
	if args.asv_to_cigar:
		print("Converting ASVs to CIGARs")
		ad.flush_dir(res_dir, "ASV_to_CIGAR", "alingments")

		path_to_seqtab = os.path.join(res_dir, 'seqtab.tsv')
		path_to_fasta = os.path.join(res_dir, "PostProc_DADA2", "ASVSeqs.fasta") #Fasta file of ASV sequences from DADA2 pipeline"
		path_to_table = os.path.join(res_dir, "PostProc_DADA2", "ASVTable.txt") #ASV table from DADA2 pipeline
		path_to_out = os.path.join(res_dir, "CIGARVariants_Bfilter.out.tsv") #Output seqtab tsv file with amplicon/variant counts
		path_asv_to_cigar = os.path.join(res_dir, "ASV_to_CIGAR", "ASV_to_CIGAR.out.txt") #Output file for ASV -> CIGAR string table 
		path_to_amp_db = "reference_panel_1.fasta" #Amplicon sequence fasta file
		path_to_alignments = os.path.join(res_dir, "ASV_to_CIGAR", "alingments") #Directory to store ASV alignment files

		print(f"INFO: Loading {path_to_amp_db}")
		amplicons = ac.parse_amp_db(path_to_amp_db)
		if not amplicons:
			print(f"ERROR: No amplicons in {path_to_amp_db}")
			sys.exit(1)

		if os.path.exists("amp_mask.txt"):
			print(f"INFO: Loading amp_mask.txt")
			mask = ac.parse_dustmasker("amp_mask.txt")
		else:
			print(f"INFO: No mask data specified.")
			mask = {}

		print(f"INFO: Loading {path_to_fasta}")
		asvs = ac.get_asv_seqs(path_to_fasta)
		if not asvs:
			print(f"ERROR: No ASV sequences in {path_to_fasta}")
			sys.exit(1)

		print(f"INFO: Parsing {path_to_table} with total reads >= {min_reads}, samples >= {min_samples}, snv_dist <= {max_snv_dist}, indel_dist <= {max_indel_dist}")

		if include_failed:
			print("WARNING: Including ASVs that failed post-DADA2 filters! This is not recommended.")
		else:
			print("INFO: Excluding ASVs that failed post-DADA2 filters.")

		if exclude_bimeras:
			print("INFO: Excluding ASVs that DADA2 marked as bimeras.")

		bins = ac.parse_asv_table(path_to_table, min_reads=min_reads, min_samples=min_samples, max_snv_dist=max_snv_dist, max_indel_dist=max_indel_dist, include_failed=include_failed, exclude_bimeras=exclude_bimeras) #This function only matches to the first strain.
		if not bins:
			print(f"ERROR: No useable data in {path_to_table}")
			sys.exit(1)

		print(f"INFO: Writing amplicon fasta files to {path_to_alignments}")
		ac.write_amplicon_fastas(asvs, bins, amplicons, outdir=path_to_alignments)

		print("INFO: Running MUSCLE aligner on amplicon fasta files. Please wait...")
		ac.run_muscle(bins, outdir=path_to_alignments)

		print("INFO: Parsing alignments to CIGAR strings")
		cigars = ac.parse_alignments(bins, mask=mask, min_homopolymer_length=polyN, outdir=path_to_alignments, verbose=False)
		if not cigars:
			print("ERROR: could not determine CIGAR strings")
			sys.exit(1)

		if path_asv_to_cigar:
			ac.write_cigar_strings(cigars, path_asv_to_cigar)
			print(f"INFO: Wrote ASV->CIGAR table to {path_asv_to_cigar}")

		print(f"INFO: Converting DADA2 seqtab file {path_to_seqtab} to {path_to_out}")
		if ac.convert_seqtab(path_to_seqtab, cigars, path_to_out):
			print("INFO: Completed successfully!")
		
if __name__ == "__main__":
	main()
