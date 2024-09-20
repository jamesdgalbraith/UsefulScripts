#!/usr/bin/env python

from re import split
import argparse
from time import time
start_time = time()
parser = argparse.ArgumentParser(description='Calculate mean coverage over windows using output of bedtools genomecov')
parser.add_argument('-i', '--depth_file', type=str, required=True,
                    help='Input bam')
parser.add_argument('-o', '--out', type=str, default='100kb_cov_win_means.txt',
                    help='Output tsv)')
parser.add_argument('-c', '--window_size', type=int, default=100000,
                    help='Size of window (default 100kbp)')
args = parser.parse_args()

window_size = args.window_size

with open(args.out, "w") as o:
    # Write header
    o.write("\t".join(["id", "pos", "cov.win.mean\n"]))
    with open(args.depth_file, "r") as in_file:
        for i, line in enumerate(in_file):
            delim = split(r"\t", line)
            # Set initial contig name and coverage count
            if i == 0:
                contig = delim[0]
                coverage_count = 0
            
            # Check if new contig, if so do maths, write and reset line counter
            if contig != delim[0]:
                # Calculate contig position, calculate mean coverage, write to file
                contig_position = int((current_contig_position//window_size)*window_size + window_size/2 )
                mean_coverage = round(coverage_count/(current_contig_position%window_size), 2)
                o.write("\t".join([contig, str(contig_position), str(mean_coverage)+"\n"])) 
                # Reset contig name, counter and position
                contig = delim[0] # SET NEW CONTIG NAME
                coverage_count = int(delim[2]) # RESET COVERAGE
                current_contig_position = int(delim[1]) # RESET COVERAGE POSITION
            # Check if end of window
            elif int(delim[1]) % window_size == 0:
                # Calculate contig position, calculate mean coverage, write to file
                contig_position = int((current_contig_position//window_size)*window_size + window_size/2)
                mean_coverage = round(coverage_count/(current_contig_position%window_size), 2)
                o.write("\t".join([contig, str(contig_position), str(mean_coverage)+"\n"])) 
                # Reset counter and position
                coverage_count = int(delim[2]) # RESET COVERAGE
                current_contig_position = int(delim[1])
            else:
                # Add coverage and record current position
                coverage_count += int(delim[2]) # ADD COVERAGE TO COUNT
                current_contig_position = int(delim[1])

    # Calculate final contig position, calculate mean coverage, write to file
    contig_position = int((current_contig_position//window_size)*window_size + window_size/2)
    mean_coverage = round(coverage_count/(current_contig_position%window_size), 2)
    o.write("\t".join([contig, str(contig_position), str(mean_coverage)+"\n"])) 

print("Calculating mean covergae over windows took "+str(time()-start_time)+ " seconds")