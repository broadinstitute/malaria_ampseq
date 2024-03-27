"""
Amplicon decontamination scripts
"""

import pandas as pd
import os, fnmatch
import subprocess
import sys
import shutil
import glob
import gzip

from Bio import SeqIO

#GENERAL USE SCRIPTS SECTION
def flush_dir(parent_dir, *dirnames):
	"""
	Remove all files and subdirectories from a directory, and create a new empty
	directory with the same name. Multiple subdirectiories may be provided.

	Args:
	parent_dir (str): The path of the parent directory.
	dirname (str): The name of the directory to be flushed.

	Returns:
	None
	"""

	dirpath = os.path.join(parent_dir, *dirnames)
	shutil.rmtree(dirpath, ignore_errors=True)
	os.makedirs(dirpath)
	return ()

def create_meta(path_to_fq, parent_dir, dirname, filename, pattern_fw, pattern_rv):
	"""
	Creates a metadata file with the list of files to process and their paths.

	Args:
	path_to_fq: string, path to the directory containing the input fastq files
	parent_dir: string, path to the parent directory where the output files will be stored
	dirname: string, name of the subdirectory where the output files will be stored
	filename: string, name of the output metadata file
	pattern_fw: string, pattern to match forward reads
	pattern_rv: string, pattern to match reverse reads

	Returns:
	None

	Example usage:

	create_meta('/path/to/input/fastq/', '/path/to/output/', 'subdir', 'metadata.tsv', '*_R1.fastq', '*_R2.fastq')

	This will search for all files in /path/to/input/fastq/ that end with '_R1.fastq' and '_R2.fastq',
	create a metadata file named 'metadata.tsv', and store it in a subdirectory named 'subdir' within 
	/path/to/output/.
	"""

	filelist = os.listdir(path_to_fq)

	dirpath = os.path.join(parent_dir, dirname) 
	outfile = os.path.join(dirpath, filename)

	meta_df = []
	for file_fw in glob.glob(os.path.join(path_to_fq, pattern_fw)):
		sampleid = os.path.basename(file_fw).split(pattern_fw[1:])[0]
		file_rv = os.path.join(path_to_fq, sampleid + pattern_rv[1:])
		if os.path.isfile(file_rv):
			meta_df.append((sampleid, file_fw, file_rv))

	with open(outfile, 'w') as f:
		for row in meta_df:
			f.write('\t'.join(row) + '\n')

	print(f"Meta file generated at location {outfile}")
	return()

def gzip_file(input_filename, output_filename):
	"""
	Compresses a file using gzip compression.

	Args:
	- input_filename (str): Path to the input file.
	- output_filename (str): Path to the output compressed file.

	Returns:
	- None

	Raises:
	- FileNotFoundError, IOError, OSError, gzip.BadGzipFile

	Compresses the input file using gzip compression and saves the result to the output file.
	The input file is removed after compression.

	Example:
	gzip_file("input.txt", "output.txt.gz")
	"""
	with open(input_filename, 'rb') as f_in:
		with gzip.open(output_filename, 'wb') as f_out:
			f_out.writelines(f_in)
	os.remove(input_filename)

#ADAPTOR AND PRIMER REMOVAL SECTION						
def adaptor_rem(sampleid, fileF, fileR, res_dir, subdir, qvalue = 5, length = 20): 
	"""
	Runs Trim Galore to remove adaptors and trim low-quality reads from paired-end fastq files.

	Args:
	sampleid (str): The base name for the output files.
	fileF (str): The path to the forward-read fastq file.
	fileR (str): The path to the reverse-read fastq file.
	res_dir (str): The path to the directory where the output files will be saved.
	subdir (str): The name of the subdirectory within the results directory where output files should be written
	qvalue (int, optional): The minimum quality score for trimming. Defaults to 5.
	length (int, optional): The minimum length of the reads to keep after trimming. Defaults to 20.

	Returns:
	None
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):
		output_dir = os.path.join(res_dir, subdir)
		cmd = ['trim_galore', '--paired', '--gzip', '--quality', f'{qvalue}', '--length', 
		f'{length}', '--output_dir', f'{output_dir}', '--basename', f'{sampleid}', f'{fileF}', 
		f'{fileR}']
		proc = subprocess.Popen(cmd)
		proc.wait()
	else:
		sys.exit('Adaptor Removal halted : one or both of the fastq files not found! Exiting...')
	return()

def trim_primer(sampleid, fileF, fileR, res_dir, subdir, pr1, pr2, prefix, keep_untrimmed=False):
	"""
	Trim primers from paired-end fastq files using cutadapt.

	Args:
	sampleid (str): Sample identifier.
	fileF (str): Path to input forward fastq file.
	fileR (str): Path to input reverse fastq file.
	res_dir (str): Path to output directory.
	subdir (str): The name of the subdirectory within the results directory where output files should be written
	pr1 (str): Path to primer sequence file for forward read.
	pr2 (str): Path to primer sequence file for reverse read.
	prefix (str): Prefix to use for output filenames.
	keep_untrimmed (bool, optional): If True, keep untrimmed reads in separate files. Default is False.

	Returns:
	None
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):

		cmd = ['cutadapt', '-g', f'file:{pr1}', '-G', f'file:{pr2}',
			'-o', os.path.join(res_dir, subdir, f'{sampleid}_{prefix}_1.fq.gz'),
			'-p', os.path.join(res_dir, subdir, f'{sampleid}_{prefix}_2.fq.gz'),
			'--pair-adapters', '--action=trim']

		if keep_untrimmed:
			cmd.extend(['--untrimmed-output', os.path.join(res_dir, subdir, f'{sampleid}_temp_1.fq.gz'),
	    			'--untrimmed-paired-output', os.path.join(res_dir, subdir, f'{sampleid}_temp_2.fq.gz')])
		else:
			cmd.append('--discard-untrimmed')

		cmd.extend([fileF, fileR])
		print(cmd)
		proc = subprocess.Popen(cmd)
		proc.wait()

	else:
		sys.exit('Pre-process halted : one or both of the fastq files not found! Exiting...')
	return()

#RUN DADA2 SECTION
def run_dada2(path_to_meta, path_to_fq, path_to_flist, Class, maxEE, trimRight, minLen, truncQ, matchIDs, max_consist, omegaA, justConcatenate, maxMismatch, saveRdata, res_dir, subdir, terra):
	"""
	Runs the DADA2 pipeline on the input files using the specified parameters.

	Args:
	path_to_meta (str): the path to the metadata file containing sample information.
	path_to_fq (str): the path to the raw fastq.gz files.
	path_to_flist (str): the path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample
	res_dir (str): the path to the directory where results will be saved.
	subdir (str): the name of the subdirectory where the output files will be saved.
	Class (str): the name of the column in the metadata file that contains the sample class information.
	maxEE (float): the maximum expected error rate.
	trimRight (int): the number of bases to trim from the 3' end of the reads.
	minLen (int): the minimum length of reads to retain after trimming.
	truncQ (int): the quality threshold for truncating reads.
	matchIDs (bool): boolean to request DADA2 to match ids on fastqs to make sure reads on forward and reverse end are in same order.
	max_consist (int): the maximum number of mismatches allowed in the overlap region for merging paired-end reads.
	omegaA (float): the alpha parameter for the consensus quality score.
	justConcatenate (int): whether to just concatenate the forward and reverse reads without merging them.
	maxMismatch (int): the maximum number of mismatches allowed during merging.
	saveRdata (str): whether to save the intermediate R data files.
	terra (bool): boolean to indicate if modified paths to terra must be used.

	Returns:
	None
	"""

	if os.path.isfile(path_to_meta):

		if subdir in ['DADA2', 'DADA2_OP', 'DADA2_NOP']:
			program = 'runDADA2.R'
		else:
			program = 'runDADA2contamination.R' 

		bimera = '--bimera'
		if terra:
			path_to_program = os.path.join("/", "Code/", program)
			platform = '--terra'
		else:	
			path_to_program = os.path.join("Code/", program)

		cmd = ['Rscript', path_to_program,
		'-p', f'{path_to_meta}',
		'-r', f'{path_to_fq}',
		'-b', f'{path_to_flist}',
		'-d', os.path.join(res_dir, subdir),
		'-o', os.path.join(res_dir, subdir, 'seqtab.tsv'),
		'-c', f'{Class}',
		'-ee', f'{maxEE}',
		'-tR', f'{trimRight}',
		'-mL', f'{minLen}',
		'-tQ', f'{truncQ}',
		'-id', f'{matchIDs}',
		'-mC', f'{max_consist}',
		'-wA', f'{omegaA}',
		'-jC', f'{justConcatenate}',
		'-mM', f'{maxMismatch}',
		'-s', f'{saveRdata}',
		f'{bimera}'] 

		if terra:
			cmd.append(f'{platform}')

		print(cmd)
		proc = subprocess.Popen(cmd)
		proc.wait()
	else:
		sys.exit('DADA2 halted : No path to meta file provided! Exiting...')
	return()

def merge_seqtab(path_op, path_nop):
	"""
	Merges overlapping and non-overlapping dada2 tables into a single table.

	Parameters:
	path_op (str): File path to the overlapping dada2 table in CSV format.
	path_nop (str): File path to the non-overlapping dada2 table in CSV format.

	Returns:
	seqtab (pd.DataFrame): Merged dada2 table containing both overlapping and non-overlapping data.
	"""

	if os.path.isfile(path_op) and os.path.isfile(path_nop):
		seqtab_op = pd.read_csv(path_op, sep = "\t", index_col=0)
		seqtab_nop = pd.read_csv(path_nop, sep = "\t", index_col=0)
		seqtab = seqtab_op.merge(seqtab_nop, left_index=True, right_index=True, how='outer')
		seqtab.reset_index(inplace=True)
	else:
		sys.exit('Overlapping and/or non-overlapping dada2 tables not found! Exiting...')

	return(seqtab)

def merge_seqtab_cont(path_op, path_nop):
	"""
	Merges overlapping and non-overlapping dada2 tables into a single table for decontamination pipeline.

	Parameters:
	path_op (str): File path to the overlapping dada2 table in CSV format.
	path_nop (str): File path to the non-overlapping dada2 table in CSV format.

	Returns:
	seqtab (pd.DataFrame): Merged dada2 table containing both overlapping and non-overlapping data.
	"""

	if os.path.isfile(path_op) and os.path.isfile(path_nop):
		seqtab_op = pd.read_csv(path_op, sep = "\t", index_col=0)
		seqtab_nop = pd.read_csv(path_nop, sep = "\t", index_col=0)
		seqtab = seqtab_op.merge(seqtab_nop, left_index=True, right_index=True, how='outer')
	else:
		sys.exit('Overlapping and/or non-overlapping dada2 tables not found! Exiting...')

	return(seqtab)

def merge_bimeras(path_op, path_nop):
	"""
	Merges overlapping and non-overlapping bimera tables from dada2 tables.

	Parameters:
	path_op (str): File path to the overlapping dada2 bimera table in CSV format.
	path_nop (str): File path to the non-overlapping dada2 bimera table in CSV format.

	Returns:
	bimeratab (pd.DataFrame): Merged dada2 bimera table containing both overlapping and non-overlapping bimera data.
	"""

	if os.path.isfile(path_op) and os.path.isfile(path_nop):
		bimera_op = pd.read_csv(path_op, sep='\t')
		bimera_nop = pd.read_csv(path_nop, sep='\t')
		bimeras = pd.concat([bimera_op, bimera_nop], ignore_index=True)
	else:
		sys.exit('Overlapping and/or non-overlapping dada2 bimera tables not found! Exiting...')

	return(bimeras)

#RUN THE AMPLICON DECONTAMINATION SECTION OF THE PIPELINE
def mergereads(sampleid, fileF, fileR, res_dir, subdir):
	"""
	This function uses bbmerge.sh to merge paired-end reads from two fastq files
	(fileF and fileR) into a single fastq file. It also generates two other fastq
	files of unmerged reads. The output files are saved in the specified res_dir
	and subdir directory paths. The function also creates a metadata file 
	(merge_meta.tsv) containing the sample ID, output file name, and standard 
	output and error logs. If either the forward or reverse fastq files are not 
	found, the function exits with an error message. This functions is optimized 
	for reads that are at the least 200bp long and amplicons 100 bp long or longer.
	Merging shorter reads will require chaning this parameters in the config file.
	
	Args:
	sampleid: a string representing the sample identifier.
	fileF: a string representing the file path of the forward reads.
	fileR: a string representing the file path of the reverse reads.
	res_dir: a string representing the directory path of the results.
	subdir: a string representing the subdirectory path of the results.

	Returns: None

	Example usage:

	mergereads("sample1", "/path/to/forward.fastq", "/path/to/reverse.fastq", "/path/to/results", "subdirectory")
	"""

	if os.path.isfile(fileF) and os.path.isfile(fileR):	
		file_nameout = os.path.join(res_dir, subdir, f"{sampleid}_stdout.txt")
		file_nameerr = os.path.join(res_dir, subdir, f"{sampleid}_stderr.txt")
		output_file_path = os.path.join(res_dir, subdir, f"{sampleid}_merged.fastq")
		output_unmerged_f_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedf.fastq")
		output_unmerged_r_path = os.path.join(res_dir, subdir, f"{sampleid}_unmergedr.fastq")
		meta_file_path = os.path.join(res_dir, subdir, f"merge_meta.tsv")

		original_stdout = sys.stdout
		original_stderr = sys.stderr

		sys.stdout = open(file_nameout, "w")
		sys.stderr = open(file_nameerr, "w")

		with open(meta_file_path, "a") as meta_file:
			meta_file.write(f"{sampleid}\t{output_file_path}\t{file_nameout}\t{file_nameerr}\n")
		
		cmd = ['fuse.sh', 
			f'in1={fileF}',
			f'in2={fileR}',
			f'out={output_file_path}',
			f'fusepairs=t',
			f'pad=10']

		print(cmd)
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()

		sys.stdout = original_stdout
		sys.stderr = original_stderr
	else:
		sys.exit('BBmerge halted : one or both of the fastq files not found! Exiting..')
	return()

def extract_bbmergefields(sampleid, mergefile, bbreportfile, path_to_flist, res_dir, rep_dir, subdir, terra):
	"""
	Extracts relevant data from a bbmerge report file and saves it to a tab-separated file.

	Args:
	sampleid: the ID of the sample being processed
	mergefile: the path to the file with the merged reads
	bbreportfile: the path to the bbmerge report file
	path_to_flist: the path to a csv file with the sample_id,Forward,Reverse, where Forward and Reverse are columns with the barcodes for the sample
	res_dir: the path to the main results directory
	rep_dir: the path to the reports directory within the results directory
	subdir: the name of the subdirectory within the results directory where output files should be written
	terra: boolean to indicate if modified paths to terra must be used.

	Returns:
	None

	Example Usage:

	extract_bbmergefields("Sample1", "/path/to/bbmerge.fastq", "/path/to/bbmerge_report.txt", "path/to/barcodes_match.csv", "/path/to/results", "/path/to/reports", "bbmerge", terra)
	"""

	if os.path.isfile(bbreportfile) and os.path.isfile(mergefile):				
		bbmergedata = {}
		bbmergedata['sampleid'] = sampleid

		if terra:
			path_to_program = os.path.join("/", "Code/", "runBBMergecontamination.R")
			platform = '--terra'
		else:	
			path_to_program = os.path.join("Code/", "runBBMergecontamination.R")

		cmd = ['Rscript', path_to_program,
		'-p', f'{mergefile}',
		'-d', os.path.join(rep_dir, subdir),
		'-b', path_to_flist]

		if terra:
			cmd.append(f'{platform}')

		print(cmd)
		proc = subprocess.Popen(cmd)
		proc.wait()
	else:
		sys.exit('Extract bbmerge report halted : bbmerge report file not found! Exiting..')
	return()

def filter_fastq_by_read_names(input_fastq_1, input_fastq_2, tsv_file, output_fastq_1, output_fastq_2):
	"""
	Filter reads from paired-end FASTQ files based on read names provided in a TSV file and save them to new output files.
	"""
	def load_read_names_from_tsv(tsv_file):
		"""
		Load read names from the first column of a TSV file.
		"""
		read_names = set()
		with open(tsv_file, 'r') as tsv:
			for line in tsv:
				read_name = line.strip().split('\t')[0]
				read_name = read_name.strip().split(' ')[0]
				match_status = line.strip().split('\t')[7]
				if match_status == 'Match':
					read_names.add(read_name[1:])
		return read_names

	def filter_fastq(input_fastq, output_fastq, read_names):
		"""
		Filter reads from input FASTQ file based on read names and save them to output FASTQ file.
		"""
		with open(output_fastq, 'w') as output_file:
			with gzip.open(input_fastq, "rt") as handle:
				for record in SeqIO.parse(handle, 'fastq'):
					if record.id.strip() in read_names:
						SeqIO.write(record, output_file, 'fastq')

		output_file_gzip = output_fastq+".gz"
		gzip_file(output_fastq, output_file_gzip)

	# Load read names from TSV file
	read_names = load_read_names_from_tsv(tsv_file)

	# Filter and save reads from paired-end FASTQ files
	filter_fastq(input_fastq_1, output_fastq_1, read_names)
	filter_fastq(input_fastq_2, output_fastq_2, read_names)

def hamming_distance(str1, str2):
	"""
	Computes the Hamming distance between two strings of equal length.

	Parameters:
	str1 (str): First input string.
	str2 (str): Second input string.

	Returns:
	int: Hamming distance between the two input strings.

	Raises:
	ValueError: If the input strings have different lengths.
	"""
	if len(str1) != len(str2):
		raise ValueError("Input strings must have the same length")

	return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))

#SEPARATE READS FOR CONTATENATION OR MERGING
#Find the adapter
def find_common_subsequence(fasta_file):
	"""
	Finds the common subsequence among sequences in a FASTA file.

	Parameters:
	fasta_file (str): Path to the input FASTA file.

	Returns:
	str: Common subsequence shared by all sequences.

	Raises:
	ValueError: If all sequences in the FASTA file are identical.
	"""
	#The adapter must be the common subsequence to all sequences in the primer file
	sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
	shortest_sequence = min(sequences, key=len)
	starting_length = len(shortest_sequence)
	# Iterate over other sequences to find common subsequence
	for sequence in sequences:
		common_subsequence = []
		for char_shortest, char_sequence in zip(shortest_sequence, sequence):
			if char_shortest == char_sequence:
				common_subsequence.append(char_shortest)
			else:
				break
		shortest_sequence = "".join(common_subsequence)

	if len(shortest_sequence) == starting_length:
		sys.exit('All sequences in the primer file are the same and the adapter cannot be determined.\nCheck your primer file or manually provide the adapter sequence.')
	else:
		return shortest_sequence

def find_longest_sequence_length(fastq_file):
	max_length = 0
	counter = 0 #There is no need to loop through the whole file. Once the max_length has remained stable over 99 iteration, the size of the read is fairly certain. Break the loop.
	with gzip.open(fastq_file, 'rt') as file:
		for line_num, line in enumerate(file, start=1):
			if line_num % 4 == 2:
				sequence_length = len(line.strip())
				max_length = max(sequence_length, max_length)
				counter = counter + 1 if sequence_length >= max_length else 0
				if counter == 99:
					break
	return max_length

def remove_adapter(fasta_file, sequence_to_remove, output_file):
	"""
	Removes a specific sequence adapter from each sequence in a FASTA file and writes the modified sequences to an output file.

	Parameters:
	fasta_file (str): Path to the input FASTA file.
	sequence_to_remove (str): Adapter sequence to remove from each sequence.
	output_file (str): Path to the output file where modified sequences will be written.
	"""
	with open(fasta_file, 'r') as input_file, open(output_file, 'w') as output_file:
		current_sequence = ''
		current_header = ''

		for line in input_file:
			if line.startswith('>'):
				if current_sequence:
					modified_sequence = current_sequence.replace(sequence_to_remove, '')
					output_file.write(f"{current_header}\n{modified_sequence}\n")

				current_header = line.strip()
				current_sequence = ''
			else:
				current_sequence += line.strip()

		if current_sequence:
			modified_sequence = current_sequence.replace(sequence_to_remove, '')
			output_file.write(f"{current_header}\n{modified_sequence}\n")

def demultiplex_per_size(sampleid, fileF, fileR, pr1, pr2, res_dir, subdir, read_size_fw, read_size_rv, asv_lengths, ci, sample_dict, mismatches = 2):
	"""
	Demultiplexes paired-end FASTQ files based on primer sequences and read sizes.

	Parameters:
	sampleid (str): Sample identifier.
	fileF (str): Path to the forward (R1) FASTQ file.
	fileR (str): Path to the reverse (R2) FASTQ file.
	pr1 (str): Path to the forward primer FASTA file.
	pr2 (str): Path to the reverse primer FASTA file.
	res_dir (str): Directory to save demultiplexed files.
	subdir (str): Subdirectory within res_dir to save sample files.
	read_size_fw (int): Expected read size for forward reads.
	read_size_rv (int): Expected read size for reverse reads.
	asv_lengths (list): List of expected lengths for ASV sequences.
	ci (bool): Subprocess for list with inline barcodes.
	sample_dict (dict): dictionary of samples per well.
	mismatches (int): Number of mismatches allowed in primer matching.

	Returns:
	None
	"""
	print("Demultiplexing: " + sampleid)

	#Retrieve Primer
	
	#TO DO: INSTEAD OF PADDING THE NAME. MAKE A NESTED DICTIONARY FOR EACH AMPLICON TARGET. OR A LIST. THIS WILL SAVE THE PROBLEM OF DEALING WITH INLINED OR NOT INLINED BARCODES.
	primer_dict_fw = {}
	primer_dict_rv = {}
	
	illumina_cols = list(range(1, 13))
	illumina_rows = list("ABCDEFGH")
	fw_asvs = []
	rv_asvs = []
	sample_s = sampleid.split("_")[1]

	print("Obtaining the primer")
	if ci:
		with open(pr1, 'r') as forward_fasta:
			for record in SeqIO.parse(forward_fasta, 'fasta'):
				fw_asvs.append(record.id)
		targets_fw = len(set(fw_asvs))
    
		with open(pr2, 'r') as reverse_fasta:
			for record in SeqIO.parse(reverse_fasta, 'fasta'):
				rv_asvs.append(record.id)
		targets_rv = len(set(rv_asvs))

		padded_names_fw = [f'{asv}_{pad}' for asv, pad in zip(fw_asvs, illumina_cols*targets_fw)]
		padded_names_rv = [f'{asv}_{pad}' for asv, pad in zip(rv_asvs, illumina_rows*targets_rv)]
		asvs = set(fw_asvs)

		with open(pr1, 'r') as forward_fasta:
			for counter, forward_record in enumerate(SeqIO.parse(forward_fasta, 'fasta')):	
				primer_dict_fw[padded_names_fw[counter]] = forward_record.seq

		with open(pr2, 'r') as reverse_fasta:
			for counter, reverse_record in enumerate(SeqIO.parse(reverse_fasta, 'fasta')):
				primer_dict_rv[padded_names_rv[counter]] = reverse_record.seq
	else:	
		dict_fw = {}
		dict_rv = {}
		with open(pr1, 'r') as forward_fasta:
			for forward_record in SeqIO.parse(forward_fasta, 'fasta'):	
				dict_fw[forward_record.id] = forward_record.seq

		with open(pr2, 'r') as reverse_fasta:
			for reverse_record in SeqIO.parse(reverse_fasta, 'fasta'):
				dict_rv[reverse_record.id] = reverse_record.seq

		primer_dict_fw = {f"{asv}_{index}": recordseq for asv, recordseq in dict_fw.items() for index in illumina_cols}
		primer_dict_rv = {f"{asv}_{index}": recordseq for asv, recordseq in dict_rv.items() for index in illumina_rows}
		asvs = set(dict_fw.keys())

	output_fastq_fw_nop = os.path.join(res_dir, subdir, f"{sampleid}_nop_L001_R1_001.fastq")
	output_fastq_rv_nop = os.path.join(res_dir, subdir, f"{sampleid}_nop_L001_R2_001.fastq")
	output_fastq_fw_op = os.path.join(res_dir, subdir, f"{sampleid}_op_L001_R1_001.fastq")
	output_fastq_rv_op = os.path.join(res_dir, subdir, f"{sampleid}_op_L001_R2_001.fastq")

	with gzip.open(os.path.join("Fastq", fileF), 'rt') as forward_fastq, gzip.open(os.path.join("Fastq", fileR), 'rt') as reverse_fastq:
		for forward_record, reverse_record in zip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq, 'fastq')):
			for asv in asvs:
				well = sample_dict[sample_s]
				row_s = well[0]
				col_s = well[1:]

				key_fw = f"{asv}_{col_s}"
				key_rv = f"{asv}_{row_s}"

				primer_fw = primer_dict_fw[key_fw]
				primer_rv = primer_dict_rv[key_rv]		
				len_fw = len(primer_fw)
				len_rv = len(primer_rv)

				forward_read_len_no_primer = len(forward_record.seq) - len_fw
				reverse_read_len_no_primer = len(reverse_record.seq) - len_rv

				#Trimgalore clips 3' ends of bad quality, which sometimes makes the record.seq shorter than the primer.
				#In such cases, the hamming distance will fail and the read must be eliminated.
				if forward_read_len_no_primer > 0 and reverse_read_len_no_primer > 0:
					hamming_distance_fw = hamming_distance(primer_fw, str(forward_record.seq[:len_fw]))
					hamming_distance_rv = hamming_distance(primer_rv, str(reverse_record.seq[:len_rv]))
					if hamming_distance_fw <= mismatches and hamming_distance_rv <= mismatches:
						usable_read_length = forward_read_len_no_primer + reverse_read_len_no_primer
						if asv_lengths[asv] < usable_read_length - 20:
							#Store reads with at least 10bp overlap
							with open(output_fastq_fw_op, 'a') as output_file_fw, open(output_fastq_rv_op, 'a') as output_file_rv:
								SeqIO.write(forward_record, output_file_fw, 'fastq')
								SeqIO.write(reverse_record, output_file_rv, 'fastq')
							break
						elif asv_lengths[asv] > usable_read_length:
							#Store reads with absolute no overlap
							with open(output_fastq_fw_nop, 'a') as output_file_fw, open(output_fastq_rv_nop, 'a') as output_file_rv:
								SeqIO.write(forward_record, output_file_fw, 'fastq')
								SeqIO.write(reverse_record, output_file_rv, 'fastq')
							break
						else:
							#Clip and store reads with intermediate overlap (between 1 and 10bps)
							seq_fw_tmp = forward_record.seq[:-10]
							seq_rv_tmp = reverse_record.seq[:-10]
							forward_record_modified = forward_record.__class__(seq_fw_tmp, id=forward_record.id, description=forward_record.description)
							reverse_record_modified = reverse_record.__class__(seq_rv_tmp, id=reverse_record.id, description=reverse_record.description)
							forward_record_modified.letter_annotations["phred_quality"] = forward_record.letter_annotations["phred_quality"][:-10]
							reverse_record_modified.letter_annotations["phred_quality"] = reverse_record.letter_annotations["phred_quality"][:-10]
							# Write the modified sequences to output files
							with open(output_fastq_fw_nop, 'a') as output_file_fw, open(output_fastq_rv_nop, 'a') as output_file_rv:
								SeqIO.write(forward_record_modified, output_file_fw, 'fastq')
								SeqIO.write(reverse_record_modified, output_file_rv, 'fastq')
							break

	    # Check and gzip the output files
	for output_file_fw, output_file_rv in [(output_fastq_fw_nop, output_fastq_rv_nop), (output_fastq_fw_op, output_fastq_rv_op)]:
		if os.path.exists(output_file_fw) and os.path.exists(output_file_rv):
			gzip_file(output_file_fw, f"{output_file_fw}.gz")
			gzip_file(output_file_rv, f"{output_file_rv}.gz")

def find_common_subsequences(records):
	"""
	Finds common subsequences among sequences with the same identifier in a list of records.

	Parameters:
	records (list): List of SeqRecord objects.

	Returns:
	dict: Dictionary mapping sequence identifiers to their common subsequences.
	"""
	sequences_by_name = {}
	for record in records:
		name = record.id
		sequence = str(record.seq)[::-1]
		if name in sequences_by_name:
			sequences_by_name[name].append(sequence)
		else:
			sequences_by_name[name] = [sequence]

	common_subsequences = {}
	for name, sequences in sequences_by_name.items():
		common_subsequences[name] = sequences[0]
		for seq in sequences[1:]:
			common_subsequences[name] = ''.join(x for x, y in zip(common_subsequences[name], seq) if x == y)
		tmp = common_subsequences[name][::-1]
		common_subsequences[name] = tmp
	return common_subsequences

def write_common_subsequences_to_fasta(common_subsequences, output_file):
	"""
	Writes common subsequences to a FASTA file.

	Parameters:
	common_subsequences (dict): Dictionary mapping sequence identifiers to their common subsequences.
	output_file (str): Path to the output FASTA file.

	Returns:
	None
	"""
	with open(output_file, 'w') as out_fasta:
		for name, sequence in common_subsequences.items():
			out_fasta.write(f'>{name}\n{sequence}\n')
