import os, argparse, sys
import pandas as pd
import utils
import dataframe_image as dfi

parser = argparse.ArgumentParser(description="Script to infer a potential co-infection")

parser.add_argument("fasta_file", help="fasta file")
parser.add_argument("--out_dir", help="Output directory", default=".")

parser.add_argument("--include", help="Get SNPs from a certain lineage", 
                    default=[], required=False, action='append')

def main():

    # Script directory
    abs_path = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(abs_path)

    # check arguments
    args = parser.parse_args()

    # Parse mutations.
    mut_dir = os.path.join(script_dir, "mutations")
    mutations = utils.parse_mut(mut_dir, args)

    # Parse reference sequence
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")
    l_ref_sequence, header, ref_sequence = utils.parse_fasta(ref_genome)

    # parse input fasta
    l_sample, header_sample, sample_sequence = utils.parse_fasta(args.fasta_file)
    fasta_name = header_sample[1:]

    # out_aln_dir
    dir_name_tsv = os.path.join(args.out_dir, fasta_name)
    utils.check_create_dir(dir_name_tsv)

    # Convert alingment to dataframe
    coordinates = list(range(1, len(l_ref_sequence) + 1))
    df_aln = pd.DataFrame([l_ref_sequence, l_sample], columns = coordinates,
        index=["Reference", fasta_name])
    df_aln_t = df_aln.T
    df_aln_SNV = df_aln_t[df_aln_t["Reference"] != df_aln_t[fasta_name]].T

    # include lineages to df
    df_concat = utils.include_lineages(args, df_aln_SNV, mutations)

    # Store df
    df_concat.T.to_csv("%s_aln.csv" %(dir_name_tsv + "/" + fasta_name), sep=",")

    # Color Dataframe
    dfi.export(df_concat.style.apply(utils.color_df, axis = 0),
                    "%s_aln.png" %(dir_name_tsv + "/" + fasta_name),
                    max_cols=-1)

if __name__ == "__main__":
    main()
