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
    else:
      [start, stop, period, copies, 
             consensus_size, perc_match, perc_indels, 
             align_score, perc_A, perc_C, perc_G, perc_T, 
             entropy, cons_seq, repeat_seq, left_flank, right_flank] = ele
      bed_line = [seq_name, start, stop]
      print '\t'.join(bed_line)

