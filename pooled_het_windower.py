def argsParser():
    import argparse
    parser = argparse.ArgumentParser(description='Calculates pooled heterozygozity given VCF generated with mpileup from a single bam')
    parser.add_argument('-v', '--in_vcf', type=str, required=True,
                        help='Input vcf')
    parser.add_argument('-d', '--min_depth', type=int, required=True,
                        help='Minimum depth for calling)')
    parser.add_argument('-D', '--max_depth', type=int, required=True,
                        help='Maximum depth for calling)')
    parser.add_argument('-o', '--out_tsv', type=str, required=True,
                        help='Output tsv')
    parser.add_argument('-w', '--window_size', type=int, default=50000,
                        help='Window size')
    parser.add_argument('-c', '--minimum_coverage', type=int, default=25,
                        help='Minimum percentage of window size to be covered for Hp to be calculated (default 25)')
    args = parser.parse_args()
    return(args.in_vcf, args.min_depth, args.max_depth, args.out_tsv, args.window_size)

def sum_maj_min(in_vcf, out_txt, min_depth, max_depth, window_size):
    from math import floor
    from re import sub
    with open(out_txt, 'w') as txt_out:
        txt_out.write("\t".join(["scaffold", "position", "major_allele", "minor_allele", "total", "n\n"]))    
        with open(in_vcf, 'r') as vcf_in:
            scaffold = ""
            bin = 0
            i = 0
            for line in vcf_in:
                if line[0] != "#":
                    delim = line.split("\t")
                    metadata = dict(x.split("=") for x in delim[7].split(";"))
                    genotype = sub("\n", "", delim[9].split(":")[0]).split("/")
                    # Check depth and if genotype called
                    # If in new bin or scaffold write to file and restart
                    if int(metadata["DP"]) >= min_depth and int(metadata["DP"]) <= max_depth and genotype[0] != "." :
                        if i == 0:
                            scaffold = delim[0]
                            bin = floor(int(delim[1])/window_size)
                            sum_maj, sum_min, n, sum_tot = 0, 0, 0, 0
                        # Update on new scaffold
                        if delim[0] != scaffold:
                            # Write old data
                            txt_out.write("\t".join([scaffold, str(bin*window_size), str(sum_maj), str(sum_min), str(sum_tot), str(n)+'\n']))
                            # Reset values
                            scaffold = delim[0]
                            bin = floor(int(delim[1])/window_size)
                            sum_maj, sum_min, n, sum_tot = 0, 0, 0, 0
                        # Update on new bin
                        elif bin != floor(int(delim[1])/window_size):
                            # Write old data
                            txt_out.write("\t".join([scaffold, str(bin*window_size), str(sum_maj), str(sum_min), str(sum_tot), str(n)+'\n']))
                            # Reset values
                            bin = floor(int(delim[1])/window_size)
                            sum_maj, sum_min, n, sum_tot = 0, 0, 0, 0
                            
                        # Check depth and if in required scaffolds if good write to file
                        depth_list = metadata["DP4"].split(',')
                        ref_n = int(depth_list[0])+int(depth_list[1])
                        alt_n = int(depth_list[2])+int(depth_list[3])
                        
                        # Determine major/minor, add to sum
                        if ref_n > alt_n:
                            sum_maj = sum_maj + ref_n
                            sum_min = sum_min + alt_n
                        else:
                            sum_maj = sum_maj + alt_n
                            sum_min = sum_min + ref_n
                        sum_tot = sum_tot + ref_n + alt_n
                        n+=1
                        i+=1
        # Write final line
        txt_out.write("\t".join([scaffold, str(bin*window_size), str(sum_maj), str(sum_min), str(sum_tot), str(n)+'\n']))

def pooledCalc(txt_out, window_size, minimum_coverage):
    import pandas as pd
    from statistics import mean, stdev
    
    # read in data
    major_minor_count = pd.read_csv(txt_out, sep = "\t")
    # filter out windows with lower than desired coverage (default 25%)
    major_minor_count = major_minor_count[major_minor_count["n"] >= (minimum_coverage/100)*window_size].copy().reset_index()
    # Calculate pooled heterozygosity
    major_minor_count['Hp'] = (2 * major_minor_count['major_allele'] * major_minor_count['minor_allele']) / pow((major_minor_count['major_allele'] + major_minor_count['minor_allele']), 2)
    # Calculate Z-score pooled heterozygosity
    major_minor_count['ZHp'] = (major_minor_count['Hp'] - mean(major_minor_count['Hp']))/stdev(major_minor_count['Hp'])
    # write to file
    major_minor_count.to_csv(txt_out, sep = "\t", header = True, index = False)

def main():
    # Get arguments
    in_vcf, min_depth, max_depth, out_txt, window_size, minimum_coverage = argsParser()
    # go along VCF
    sum_maj_min(in_vcf, out_txt, min_depth, max_depth, window_size)
    # Calculate Hp and ZHp from initial output
    pooledCalc(out_txt, window_size, minimum_coverage)

if __name__ == '__main__':
    __version__ = '0.2'
    try:
        main()
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user\n")
        exit(-1)