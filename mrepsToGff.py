def argsParser():

    import argparse

    parser = argparse.ArgumentParser(description='Converts mreps output to gff and fasta')

    parser.add_argument('-i', '--mreps', type=str, required=True,
                        help='Input mreps data')
    parser.add_argument('-g', '--gff', type=str, required=True,
                        help='Path to GFF output')
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='Path to FASTA output')
    
    args = parser.parse_args()
    return(args.mreps, args.gff, args.fasta)

def mrepsConverter(mreps, gff, fasta):
    
    from re import sub

    with open(mreps, 'r') as in_mreps, open(gff, 'w') as out_gff, open(fasta, 'w') as out_fasta:
        scaffold=""
        for line in in_mreps:
            if "Processing sequence" in line:
                scaffold=sub("'", "", line.split()[2])
            if "->" in line and "from" not in line:
                line_split=line.split()
                start = line_split[0]
                end = line_split[2]
                width = line_split[4]
                period = sub("<", "", sub(">", "", line_split[5]))
                n_expected = sub("\\[", "", sub("\\]", "", line_split[6]))
                error=line_split[7]
                sequence=line_split[8]
                out_gff.write("\t".join([scaffold, "mreps", "Satellite", start, end, ".", "+", ".", "Period="+period+";Number="+n_expected+";Error="+error+"\n"]))
                out_fasta.write(">"+scaffold+":"+start+"-"+end+"_(Period="+period+")\n")
                out_fasta.write(line_split[8]+"\n")

def main():
    # Get arguments
    mreps, gff, fasta = argsParser()

    # Convert mreps to gff and fasta
    mrepsConverter(mreps, gff, fasta)

if __name__ == '__main__':
    __version__ = '0.1'
    try:
        main()
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user\n")
        exit(-1)