"""
Align ASVs to target amplicon reference and report variants as CIGAR strings
"""
import argparse
import sys
import os
import subprocess
import re

from Bio import SeqIO, AlignIO

# parse amplicon dust mask info
def parse_dustmasker(mask_info):
	"""
	Parse DUST accloc format mask info
	Description:
	This function parses DUST accloc format mask information from a file and organizes it into a dictionary. It extracts gene names, start positions, and end positions from the input file and constructs a dictionary where each gene is associated with a set of positions that are masked.

	Parameters:
	mask_info (str): The file path to the DUST accloc format mask information file.
Returns:
	mask (dict): A dictionary containing the parsed mask information. The keys are gene names, and the values are sets of positions that are masked for each gene.
	"""
	if not mask_info:
		return
	mask = {}
	with open(mask_info) as f:
		for line in f:
			line = line.strip().split("\t")
			gene = line[0].split(":")[0][1:]
			if gene not in mask:
				mask[gene] = set()
			start = int(line[1])+1 # mask info is 0-based, but we want 1-based
			end = int(line[2])+2 # +1 for 1-based and +1 to include last pos in range
			mask[gene].update(list(range(start, end))) # add all pos in between start and end
	if not mask:
		print("ERROR: No mask data loaded! Is the file the correct format?")
		sys.exit(1)
	return mask

# parse amplicon database
def parse_amp_db(fasta_file):
	"""
	Load sequences from fasta file of amplicons
	Description:
	This function loads sequences from a FASTA file containing amplicons and creates a dictionary where the sequence identifiers (IDs) are associated with their respective sequence objects.

	Parameters:
	fasta_file (str, optional): The path to the FASTA file containing amplicon sequences.
	Returns:
	amplicons (dict): A dictionary where the keys are sequence IDs and the values are the corresponding sequence objects.
	"""

	amplicons = {}
	for seq in SeqIO.parse(fasta_file, "fasta"):
		amplicons[seq.id] = seq
	return amplicons

# parse asv to amplicon table
def parse_asv_table(file, min_reads=0, min_samples=0, max_snv_dist=-1, max_indel_dist=-1, include_failed=False, exclude_bimeras=False):
	"""
	Parse DADA2 ASV table format
	Description:
	This function parses a DADA2 ASV table in a specific format. It filters the data based on various criteria such as minimum number of reads, minimum number of samples, maximum SNV (Single Nucleotide Variant) distance, maximum indel (insertion/deletion) distance, inclusion/exclusion of failed reads, and exclusion of bimeras. The function organizes the ASVs (Amplicon Sequence Variants) into bins based on the target gene or amplicon they belong to.

	Parameters:

	file (str): The path to the ASV table file.
	min_reads (int, optional): The minimum number of reads required for an ASV to be considered. Default is 0.
	min_samples (int, optional): The minimum number of samples an ASV must be present in to be considered. Default is 0.
	max_snv_dist (int, optional): The maximum allowed SNV distance for an ASV. ASVs with a greater distance will be excluded. Default is -1, indicating no threshold.
	max_indel_dist (int, optional): The maximum allowed indel distance for an ASV. ASVs with a greater distance will be excluded. Default is -1, indicating no threshold.
	include_failed (bool, optional): Whether to include ASVs that have failed post-DADA2 filters. Default is False.
	exclude_bimeras (bool, optional): Whether to exclude ASVs that have been identified as bimeras by DADA2. Default is False.
	Returns:

	bins (dict): A dictionary where the keys are amplicons (target genes) and the values are lists of ASVs assigned to each amplicon.
	"""

	bins = {}
	with open(file) as f:
		f.readline()
		for line in f:
			line = line.strip().split("\t")
			# total reads
			print(line)
			if int(line[2]) < min_reads: 
				continue # skip if too few total reads
			# total samples
			if int(line[3]) < min_samples: 
				continue # skip if in too few samples
			# minimum SNV distance
			if max_snv_dist >= 0 and int(line[6]) > max_snv_dist:
				continue # skip if snv distance > threshold
			# minimum indel distance
			if max_indel_dist >= 0 and int(line[7]) > max_indel_dist:
				continue # skip if indel distance > threshold
			# check for failing the snv_filter and indel_filter
			if not include_failed and (line[-3] == "FAIL" or line[-2] == "FAIL"):
				continue # failed post-DADA2 filters
			# check for dada2 calling asv a bimera
			if exclude_bimeras and line[-1] == "TRUE":
				continue # skip if dada2 called bimera
			ASV = line[0] # (e.g. ASV123)
			amplicon = line[5] # target gene/amplicon
			if amplicon not in bins:
				bins[amplicon] = []
			bins[amplicon].append(ASV)
	return bins # bins is dict of amplicon -> list of ASVs assigned to amplicon

# parse ASV fasta file
def get_asv_seqs(file):
	"""
	Load ASV sequences from fasta file
	Description:
	This function loads ASV sequences from a FASTA file and returns them as a dictionary where the sequence identifiers (IDs) are associated with their respective sequence objects.

	Parameters:

	file (str): The path to the FASTA file containing ASV sequences.
	Returns:

	asv_seqs (dict): A dictionary where the keys are sequence IDs and the values are the corresponding sequence objects.
	"""
	return {seq.id: seq for seq in SeqIO.parse(file, "fasta")}

# write amplicon fasta files
def write_amplicon_fastas(asvs, bins, amplicons, outdir="ASVs"):
	"""
	Write one fasta file per amplicon, containing reference sequence and assigned ASVs
	Description:
	This function writes one FASTA file per amplicon, including the reference sequence and assigned ASVs (Amplicon Sequence Variants). Each amplicon file is created in a specified output directory.

	Parameters:

	asvs (dict): A dictionary containing the ASVs, where the keys are ASV IDs and the values are the corresponding sequence objects.
	bins (dict): A dictionary where the keys are amplicons (target genes) and the values are lists of ASVs assigned to each amplicon.
	amplicons (dict): A dictionary containing the amplicon reference sequences, where the keys are amplicon names and the values are the corresponding sequence objects.
	outdir (str, optional): The path to the output directory where the amplicon FASTA files will be written. Default is "ASVs".
	Returns:
	None

	"""

	for amplicon in bins:
		if amplicon not in amplicons:
			print(f"WARNING: {amplicon} target not found in amplicon sequence database")
			continue
		with open(os.path.join(outdir, f"{amplicon}.fasta"), "w") as w:
			SeqIO.write(amplicons[amplicon], w, "fasta")
			SeqIO.write([asvs[asv] for asv in bins[amplicon]], w, "fasta")

# run muscle for each amplicon
def run_muscle(bins, outdir="ASVs"):
	"""
	Iterate through amplicons, aligning each one with MUSCLE
	Description:
	This function iterates through the amplicons in the bins dictionary and performs multiple sequence alignment (MSA) for each amplicon using the MUSCLE tool. The aligned sequences are saved as MSA files in a specified output directory.

	Parameters:

	bins (dict): A dictionary where the keys are amplicons (target genes) and the values are lists of ASVs assigned to each amplicon.
	outdir (str, optional): The path to the output directory where the MSA files will be saved. Default is "ASVs".
	Returns:
	None

	"""
	for amplicon in bins:
		fasta = os.path.join(outdir, f"{amplicon}.fasta")
		if not os.path.isfile(fasta):
			print(f"ERROR: Could not find {fasta}")
			continue
		msa = os.path.join(outdir, f"{amplicon}.msa")
		subprocess.run(["muscle", "-in", fasta, "-out", msa], capture_output=True)

# get coords of homopolymer runs
def _get_homopolymer_runs(seq, min_length=5):
	"""
	Detect and report homopolymer runs of minimum length
	Description:
	This function detects and reports homopolymer runs of a minimum length in a given sequence. A homopolymer run refers to consecutive repeated nucleotides or characters in the sequence.

	Parameters:

	seq (str): The input sequence in which homopolymer runs are to be detected.
	min_length (int, optional): The minimum length of the homopolymer run to be reported. Default is 5.
	Returns:

	runs (set): A set of indices representing the positions of the detected homopolymer runs in the input sequence.

	"""
	runs = set()
	prev = None
	run = 1
	start = None
	last_non_gap = None
	for i in range(len(seq)):
		if seq[i] == "-":
			continue
		if seq[i] == prev:
			if not start:
				if i > 1 and seq[i-2] == '-':
					# gap at start of run
					j = i - 2
					while j >= 0:
						if seq[j] != "-":
							start = j+1 # start is the start of the gap
							break
						j -= 1
					else:
						start = 0
				else:
					start = last_non_gap
			run += 1
		else:
			if run >= min_length:
				runs.update(list(range(start, i)))
			run = 1
			start = None
		prev = seq[i]
		last_non_gap = i
	
	return runs

# parse muscle alignment
def parse_alignment(alignment, mask={}, min_homopolymer_length=5, amplicon=None, verbose=False):
	"""
	Parse amplicon alignment file, converting ASV to CIGAR string
	
	Description:
	This function parses an amplicon alignment file and converts ASV (Amplicon Sequence Variant) alignments to CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings. It processes the alignment, identifies indels, and handles masking of positions. The function returns a dictionary mapping ASV IDs to their respective CIGAR strings.

	Parameters:

	alignment (str): The path to the amplicon alignment file.
	mask (dict, optional): A dictionary containing masking information for specific amplicon sequences. Default is an empty dictionary ({}).
	min_homopolymer_length (int, optional): The minimum length of a homopolymer run to be reported. Default is 5.
	amplicon (str, optional): The ID of the amplicon reference sequence. Default is None.
	Returns:

	asv_to_cigar (dict): A dictionary mapping ASV IDs to their corresponding CIGAR strings.
	"""
	aln = AlignIO.read(alignment, "fasta")
	# sort such that amplicon reference is first in alignment
	aln.sort(key = lambda record: (record.id != amplicon, record.id))
	anchor = aln[0]
	if anchor.id != amplicon:
		print(f"ERROR: No anchor gene for {alignment}")
		# don't parse if amplicon reference not in alignment (this shouldn't happen)
		return

	if min_homopolymer_length > 1:
		# detect homopolymer runs in reference sequence
		homopolymer_runs = _get_homopolymer_runs(aln[0], min_length=min_homopolymer_length)

	if len(anchor.seq.lstrip("-")) != aln.get_alignment_length():
		print(f"WARNING: {os.path.basename(alignment)} extends beyond 5' end of reference sequence!")
	elif len(anchor.seq.rstrip("-")) != aln.get_alignment_length():
		print(f"WARNING: {os.path.basename(alignment)} extends beyond 3' end of reference sequence!")

	masked = mask.get(aln[0].id, None)

	asv_to_cigar = {}
	for seq in aln[1:]:
		pos = 1 # start at position 1 in anchor sequence
		cigar = ""  # cigar string to output, start empty
		indel = False # indicate alignment column in an indel
		masking = False # indicate alignment column is being masked
		for i in range(aln.get_alignment_length()):
			# if anchor pos masked, or next base in anchor is masked and anchor position is a gap
			if masked and (pos in masked or (pos+1 in masked and anchor[i] == '-')):
				if verbose and seq.id == aln[1].id:
					if not masking:
						print(f"INFO: Skipping masked positions starting at {pos} in {os.path.basename(alignment)}")
						if anchor[i] == '-':
							print(f"INFO: Gap in alignment at start of masked region!")
						masking = True
					elif pos not in masked:
						print(f"INFO: Ending masked positions at {pos-1} in {os.path.basename(alignment)}")
						masking = False
			elif min_homopolymer_length > 1 and i in homopolymer_runs:
				if verbose and seq.id == aln[1].id:
					if i and i-1 not in homopolymer_runs:
						print(f"INFO: Skipping homopolymer run (poly-{anchor[i]}) beginning at position {pos} in {os.path.basename(alignment)}")
					elif i+1 not in homopolymer_runs:
						print(f"INFO: End of homopolymer run (poly-{anchor[i]}) at position {pos} in {os.path.basename(alignment)}")
			elif seq[i] != anchor[i]:
				if anchor[i] == "-":
					if not indel:
						indel = True
						cigar += f"{pos}I="
						if i:
							for j in range(1,len(anchor)-i):
								if anchor[i-j] != "-":
									cigar += anchor[i-j]
									break
					cigar += seq[i]
				elif seq[i] == "-":
					if not indel:
						indel = True
						cigar += f"{pos}D="
					cigar += f"{anchor[i]}"
				else:
					cigar += f"{pos}{seq[i]}"
					indel = False
			else:
				indel = False
			if anchor[i] != '-':
				pos += 1

		if not cigar:
			cigar = "."
		asv_to_cigar[seq.id] = cigar
	return asv_to_cigar

# get variants per amplicon per position
def parse_alignments(bins, mask={}, min_homopolymer_length=5, outdir="ASVs", verbose=False):
	"""
	Parse multi-sequence alignment fasta file from MUSCLE
	Description:
	This function parses a multi-sequence alignment FASTA file generated by MUSCLE and extracts CIGAR strings for each amplicon. It processes the alignments for each amplicon using the parse_alignment() function and returns a dictionary containing the amplicon IDs as keys and their corresponding CIGAR strings as values.

	Parameters:

	bins (dict): A dictionary containing the amplicons as keys and the associated ASVs (Amplicon Sequence Variants) as values. The ASVs can be stored as a list, set, or any iterable.
	mask (dict, optional): A dictionary containing masking information for specific amplicon sequences. Default is an empty dictionary ({}).
	min_homopolymer_length (int, optional): The minimum length of a homopolymer run to be reported. Default is 5.
	outdir (str, optional): The directory where the MUSCLE alignment files are located. Default is "ASVs".
	Returns:

	cigars (dict): A dictionary mapping amplicon IDs to their respective CIGAR strings.

	"""
	cigars = {}
	for amplicon in sorted(bins):
		msa = os.path.join(outdir, f"{amplicon}.msa")
		if not os.path.isfile(msa):
			print(f"ERROR: Could not find {msa}")
			continue
		# store CIGAR strings per amplicon in dict
		cigars[amplicon] = parse_alignment(msa, mask=mask, min_homopolymer_length=min_homopolymer_length, amplicon=amplicon, verbose=False)
	
	return cigars

# write table of asv -> amplicon/cigar
def write_cigar_strings(cigars, out):
	"""
	Write conversion table from ASV to CIGAR string
	Description:
	This function writes a conversion table from ASV (Amplicon Sequence Variant) to CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings. It takes a dictionary of CIGAR strings for each amplicon and writes them to a specified output file in a tab-separated format. Each line of the output file contains the ASV ID, the corresponding amplicon target, and the CIGAR string.

	Parameters:

	cigars (dict): A dictionary containing the CIGAR strings for each amplicon. The keys are amplicon IDs, and the values are dictionaries mapping ASV IDs to their respective CIGAR strings.
	out (str): The path of the output file to write the conversion table.
	Returns:

	None
	"""
	number = re.compile(r"\d+")
	with open(out, 'w') as w:
		# write tab file with ASV, amplicon target, and CIGAR string
		w.write("ASV\tAmplicon\tCIGAR\n")
		for amplicon in sorted(cigars):
			# sort on ASV number
			for ASV in sorted(cigars[amplicon], key = lambda x: int(number.search(x).group())):
				w.write(f"{ASV}\t{amplicon}\t{cigars[amplicon][ASV]}\n")


def convert_seqtab(file, cigars, out):
	"""
	Parse seqtab file, converting ASVs to CIGAR strings
	Description:
	This function parses a seqtab file, which contains ASV counts per sample, and converts the ASVs to their corresponding CIGAR strings using a provided dictionary of CIGAR strings for each amplicon. It then writes the converted data to an output file in tab-separated format, with each line representing a sample and its associated CIGAR string counts.

	Parameters:

	file (str): The path to the seqtab file to be parsed and converted.
	cigars (dict): A dictionary containing the CIGAR strings for each amplicon. The keys are amplicon IDs, and the values are dictionaries mapping ASV IDs to their respective CIGAR strings.
	out (str): The path of the output file to write the converted data.
	Returns:

	True if the conversion and writing process completes successfully.
	"""

	# get dict of ASVs -> amplicon/CIGAR
	asv_to_cigar = {}
	variants = set()
	for amplicon in sorted(cigars):
		for ASV in sorted(cigars[amplicon]):
			variant = f"{amplicon},{cigars[amplicon][ASV]}"
			asv_to_cigar[ASV] = variant
			variants.add(variant)
	
	if not variants:
		print("ERROR: No haplotypes to convert!")
		return
	
	total_reads = {}
	# parse seqtab file
	with open(file) as f:
		seqtab = {}
		f.readline()
		for line in f:
			line = line.strip().split("\t")
			sample = line[0]
			seqtab[sample] = {}
			# iterate through each ASV (i.e. ASV1, ASV2, ... ASVn)
			print(line)
			print(line[1:])
			for i, count in enumerate(line[1:]):
				asv = f"ASV{i+1}" # don't use actual sequence
				variant = asv_to_cigar.get(asv)
				if not variant:
					continue # ASV was filtered out
				if count == '': 
					count = '0.0' #Correct NA entries left in the seqtab
				# sum ASVs per sample that are the same variant
				count = float(count)
				if variant not in total_reads:
					total_reads[variant] = 0
				total_reads[variant] += count
				if variant not in seqtab[sample]:
					seqtab[sample][variant] = 0
				seqtab[sample][variant] += count
		
		if not seqtab:
			print("ERROR: No seqtab data to write!")
			return

		# write output file (sort variants first)
		variants = sorted(list(variants), key=lambda variant: total_reads.get(variant, 0), reverse=True)
		with open(out, "w") as w:
			# write header
			w.write("sample\t" + "\t".join(variants) + "\n")
			# write one sample per line
			for sample in sorted(seqtab):
				w.write(f"{sample}\t" + "\t".join([f"{seqtab[sample].get(variant, 0)}" for variant in variants]) + "\n")
		
			return True

