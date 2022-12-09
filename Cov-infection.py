# imports
import argparse, os, utils, mixtas, sys

# parse arguments (-B)
parser = argparse.ArgumentParser(description="Script to infer a potential co-infection")

# INPUT parameters
parser.add_argument("bamfile", help="bam file to extract htz positions")
parser.add_argument("covfile", help="cov file to extract coverage")
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--file_sep", help="File separator", default="\t")

# OPTIONAL
parser.add_argument("--min_DP", help="minimum frequency (depth) to accept a SNP", default=10,
                    type=int)
parser.add_argument("--min_HOM", help="minimum proportion for homocygosis",
                    default=0.85, type=float)
parser.add_argument("--ambiguity", help="min value to segregate", default=0.55, type=float)
parser.add_argument("--max_mean_htz", help="maximum htz proportion",
                    default=0.75, type=float)
parser.add_argument("--max_std_htz", help="maximum deviation of htz proportion",
                    default=0.08, type=float)
parser.add_argument("--max_extra_std_htz", help="maximum extra deviation of htz proportion",
                    default=0.015, type=float)
parser.add_argument("--SNPs_out", help="percentage of SNPs out mean htz +- (std htz + extra_std)",
                    default=0.3, type=float)

parser.add_argument("--pangolin", help="pangolin annotation", 
                    action="store_true")
parser.add_argument("--snipit", help="snipit analysis", 
                    action="store_true")
parser.add_argument("--episode", help="Extra recent episodes (samples) as fasta sequences", 
                    default=[], required=False, action='append')

# main
def main():

    # Script directory
    abs_path = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(abs_path)

    # check arguments
    args = parser.parse_args()

    if utils.check_argmunets(args, script_dir):
        exit(1)

    else:

        # Check out dir
        utils.check_create_dir(args.out_dir)

        # bam file (SNPs)
        ap_bamfile = os.path.abspath(args.bamfile)
        name_bam = os.path.basename(ap_bamfile).rstrip(".rg.markdup.trimmed.sorted.bam")
        
        # get tsv from bam file
        tsv_file = mixtas.bam2tsv(args.bamfile, args, script_dir, name_bam)
        name_tsv = os.path.basename(tsv_file).rstrip(".tsv")
        dir_name_tsv = os.path.join(args.out_dir, name_tsv)

        # Parse mutations.
        mut_dir = os.path.join(script_dir, "mutations")
        mutations = utils.parse_mut(mut_dir, args)

        # parse coverage file
        cov_d = utils.parse_covfile(args.covfile)

    # parse variant calling file
    df = mixtas.parse_vcf(args, tsv_file, name_tsv, mutations)

    # get alingment
    mixtas.get_alingment(args, script_dir, name_tsv, df, mutations, cov_d)

    # include previous or post episodes
    if args.episode:
        mixtas.compare_episode(args, name_tsv, script_dir)

    # Get stats
    utils.quality_control(df, args, mutations, name_tsv, dir_name_tsv)
    
    # If all OK
    exit(0)

if __name__ == "__main__":
    main()
