# adapted from https://gist.github.com/lexnederbragt/3689ee2301493c34c8ab#file-trf2gff-py

# Library  needed for dealing with system
import sys

# Open file
with open(sys.argv[1]) as infile:
  # All below loops over file
  for line in infile:
    # Split input into columns based on spaces
    ele = line.strip().split(" ")
    # Grab sequence name from lines beginning with "@"
    if line.startswith('@'):
      seq_name = ele[0][1:]
    # Print sequence name and select other columns for other lines
    elif ele[0].isdigit():
      [lstart, lstop, lperiod, rstart, rstop, rperiod, loop, pcmatches, pcindels, score, pcAT, pcGC, pcATpair, pcGCpair, pcGTpair, ctr2, av_ctr2, lseq, rseq] = ele
      bed_line = [seq_name, lstart, rstop]
      print('\t'.join(bed_line))
