def argsParser():
    import argparse
    parser = argparse.ArgumentParser(description='Counts hetrozygous and homozygous sites from VCF\ngenerated with mpileup (from single BAM)')
    parser.add_argument('-v', '--in_vcf', type=str, required=True,
                        help='Input vcf')
    parser.add_argument('-d', '--min_depth', type=int, required=True,
                        help='Minimum depth for calling)')
    parser.add_argument('-D', '--max_depth', type=int, required=True,
                        help='Maximum depth for calling)')
    parser.add_argument('-o', '--out_tsv', type=str, required=True,
                        help='Output tsv')
    parser.add_argument('-w', '--window_size', type=int, default=100000,
                        help='Window size')
    args = parser.parse_args()
    return(args.in_vcf, args.min_depth, args.max_depth, args.out_tsv, args.window_size)

def site_counter(in_vcf, out_txt, min_depth, max_depth, binning_width):
    from math import floor
    from re import sub
    with open(out_txt, 'w') as txt_out:
        txt_out.write("\t".join(["scaffold", "position", "homozygous_sites", "heterozygous_sites", "not_called\n"]))
        with open(in_vcf, 'r') as vcf_in:
            scaffold = ""
            bin = 0
            i = 0
            for line in vcf_in:
                if line[0] != "#":
                    delim = line.split("\t")
                    metadata = dict(x.split("=") for x in delim[7].split(";"))
                    # Check depth and if in required scaffolds if good write to file
                    if int(metadata["DP"]) >= min_depth and int(metadata["DP"]) <= max_depth:
                        if i == 0:
                            scaffold = delim[0]
                            bin = floor(int(delim[1])/binning_width)
                            hom, het, not_called = 0, 0, 0
                        # Update on new scaffold
                        if delim[0] != scaffold:
                            # Write old data
                            txt_out.write("\t".join([scaffold, str(bin*binning_width), str(hom), str(het), str(not_called)+'\n']))
                            # Reset values
                            scaffold = delim[0]
                            bin = floor(int(delim[1])/binning_width)
                            hom, het, not_called = 0, 0, 0
                        # Update on new bin
                        elif bin != floor(int(delim[1])/binning_width):
                            # Write old data
                            txt_out.write("\t".join([scaffold, str(bin*binning_width), str(hom), str(het), str(not_called)+'\n']))
                            # Reset values
                            bin = floor(int(delim[1])/binning_width)
                            hom, het, not_called = 0, 0, 0

                        # Check depth and if in required scaffolds if good write to file
                        if int(metadata["DP"]) >= min_depth and int(metadata["DP"]) <= max_depth:
                            genotype = sub("\n", "", delim[9].split(":")[0]).split("/")
                            if genotype[0] != "." and genotype[0] != genotype[1]:
                                het+=1
                            elif genotype[0] != "." and genotype[0] == genotype[1]:
                                hom+=1
                            else:
                                not_called+=1
                        i+=1
        # Write final line
        txt_out.write("\t".join([scaffold, str(bin*binning_width), str(hom), str(het), str(not_called)+'\n']))

def main():
    # Get arguments
    in_vcf, min_depth, max_depth, out_txt, binning_width, = argsParser()
    # go along VCF
    site_counter(in_vcf, out_txt, min_depth, max_depth, binning_width)

if __name__ == '__main__':
    __version__ = '0.3'
    try:
        main()
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user\n")
        exit(-1)