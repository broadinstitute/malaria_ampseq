#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def main():
    """
    Gets the marker table from a FASTA or a BED file of the amplicons; this markers table is required for
    the MHap-Analysis pipeline.

    Usage: python Code/markersTable_from_bed.py -i <BED/FASTA file> -o <output file>

    Returns:
    None.
    """

    ## Load arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help="Input path to amplicon BED/FASTA file.", required=True)
    parser.add_argument('-o', '--output', help="Path for markers table output", required=True)
    args = parser.parse_args()

    if args.input.lower().endswith(".fasta"):
        # Assumes that the FASTA ID is of the format "{chromosome}:{start}-{end}"
        amplicons = []
        chroms = []
        starts = []
        ends = []
        lengths = []
        positions = []
        with open(args.input, 'r') as f:
             for line in f:
                  line = line.strip()
                  if line.startswith('>'):
                       header = line[1:]
                       h_split = header.split(':')
                       pos_split = h_split[1].split('-')
                       chrom = h_split[0]
                       start = int(pos_split[0])
                       end = int(pos_split[1])
                       pos = (start + end) // 2
                       length = end - start + 1

                       amplicons.append(header)
                       chroms.append(chrom)
                       starts.append(start)
                       ends.append(end)
                       positions.append(pos)
                       lengths.append(length)
                                     
        markers_table = pd.DataFrame({
             'amplicon': amplicons,
             'chromosome': chroms,
             'start': starts,
             'end': ends,
             'pos': positions,
             'length': lengths     
        })

    elif args.input.lower().endswith(".bed"):
        amplicon_bed = pd.read_csv(filepath_or_buffer=args.input, sep=None, engine="python", header=None, names=["chromosome", "start", "end"])
        amplicon_bed['amplicon'] = amplicon_bed['chromosome'] + ":" + amplicon_bed["start"].astype(str) + "-" + amplicon_bed["end"].astype(str)
        amplicon_bed['pos'] = (amplicon_bed['start'] + amplicon_bed['end'] ) // 2
        amplicon_bed['length'] = amplicon_bed['end'] - amplicon_bed['start'] + 1

        markers_table = amplicon_bed[['amplicon', 'chromosome', 'start', 'end', 'pos', 'length']]
    else:
         print("ERROR: Please provide either a FASTA or BED file of your amplicon regions.")
         sys.exit(1)
    
    markers_table.to_csv(args.output, index=False)

if __name__ == "__main__":
	main()