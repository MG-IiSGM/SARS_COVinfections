import argparse
import os
import utils
import mixtas
import sys

# avoid __pycache__ folder
sys.dont_write_bytecode = True

# parse arguments (-B)
parser = argparse.ArgumentParser()

parser.add_argument("tsv_file", help="tsv file to extract htz positions")
parser.add_argument("--mut_dir", help="directory where all mutation file are located",
                     required=False, default="../mutations")
parser.add_argument("--ref_genome", help="file with reference genome", required=False,
                    default="../COVID_ref.fasta")
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--file_sep", help="File separator", default="\t")

parser.add_argument("--get_SNP", help="Get SNVs from a certain lineage", 
                    default=[], required=False, action='append')
parser.add_argument("--episode", help="Extra recent episodes (samples) as fasta sequences", 
                    default=[], required=False, action='append')
parser.add_argument("--min_DP", help="minimum frequency (depth) to accept a SNV", default=15,
                    type=int)
parser.add_argument("--min_HOM", help="minimum proportion for homocigosis",
                    default=0.85, type=float)
parser.add_argument("--ambiguity", help="min value to segregate", default=0.56, type=float)

parser.add_argument("--pangolin", help="pangolin annotation", 
                    action="store_true")
parser.add_argument("--snipit", help="snipit analysis", 
                    action="store_true")

# main
def main():

    # check arguments
    args = parser.parse_args()

    if utils.check_argmunets(args):
        exit(1)

    else:

        # Check out dir
        utils.check_create_dir(args.out_dir)

        # tsv file (SNPs)
        name_tsv = os.path.basename(args.tsv_file).rstrip(".tsv")

        # Parse mutations.
        mutations = utils.parse_mut(args)

    # main function
    df = mixtas.get_lineage(args, name_tsv, mutations)

    # get HTZ/HOM positions
    HTZ_SNVs, HOM_SNVs = mixtas.get_HTZ(df, args, name_tsv, mutations)

    # get alingment
    mixtas.get_alingment(args, name_tsv, HTZ_SNVs, HOM_SNVs, df, mutations)

    # include previous or post episodes
    if args.episode:
        mixtas.compare_episode(args, name_tsv)

if __name__ == "__main__":
    main()
