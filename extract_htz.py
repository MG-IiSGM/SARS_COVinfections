# imports
import os
import pandas as pd
import numpy as np
import utils
import argparse
import warnings
warnings.filterwarnings("ignore")

# parse arguments (-B)
parser = argparse.ArgumentParser()

parser.add_argument("--tsv_file", help="tsv file to extract htz positions", 
                    default=[], required=True, action='append')
parser.add_argument("--cov_file", help="tsv file to extract htz positions", 
                    default=[], required=True, action='append')
parser.add_argument("--mut_dir", help="directory where all mutation file are located", required=False,
                    default="/media/NASII/Datos/ANALYSIS/MISC/covid_analysis/coinfecciones/PROYECTO_MIXTAS/mutations")
parser.add_argument("--ref_genome", help="file with reference genome", required=False,
                    default="/media/NASII/Datos/ANALYSIS/MISC/covid_analysis/coinfecciones/PROYECTO_MIXTAS/SCRIPTS/COVID_ref.fasta")

parser.add_argument("--get_SNVs", help="Get SNVs from a certain lineage", 
                    default=[], required=False, action='append')
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--out_name", help="alignment out name", default="alingment")
parser.add_argument("--file_sep", help="File separator", default="\t")
parser.add_argument("--min_DP_freq", help="minimum frequency (depth) to accept a SNV", default=10)
parser.add_argument("--min_homo_prop", help="minimum proportion for htz", default=0.92)

# check arguments
args = parser.parse_args()
if utils.check_argmunets(args):
    exit(1)
else:
    # Parse mutations.
    mutations = utils.parse_mut(args)
    # Check out dir
    utils.check_create_dir(args.out_dir)
    # alignment out_file
    out_aln = open(args.out_dir + "/" + args.out_name + "_aln.fasta", "w")
    # parse coverage to dictionary
    cov_d = utils.parse_coverage_file(args)

to_write = ""
for file in args.tsv_file:

    # Parse reference sequence
    l_ref_sequence, header, ref_sequence = utils.parse_ref_genome(args)

    # Filter SNVs with a ALT_DP < min_freq
    # Add genome SNV genome location
    df = utils.filter_tsv(args, file, drop_features=True)

    # Label SNVs to lineage
    df["LINEAGE"] = utils.indetify_variants(df, mutations)
    df.fillna(0, inplace=True)

    # Store
    name_file = os.path.basename(file).split(".")[0]
    df.to_csv("%s_lineage.csv" %(args.out_dir + "/" + name_file), index=False, sep=",")
    df.explode("LINEAGE").to_csv("%s_lineage_sep.csv" %(args.out_dir + "/" + name_file), index=False, sep=",")

    # General stats
    utils.get_general_stats(df, name_file, args.out_dir, args.min_homo_prop)
    # Generate alingment
    indel_range = []
    total_dp = 0
    for _, row in df.iterrows():

        coordinate = row["POS"] - 1

        if coordinate in indel_range and row["TOTAL_DP"] < total_dp:
            continue

        elif row["ALT_FREQ"] >= args.min_homo_prop:
            if "-" in row["ALT"] or "+" in row["ALT"]:
                if row["ALT"][1] == "A":
                    l_ref_sequence[coordinate] = "T"
                elif row["ALT"][1] == "T":
                    l_ref_sequence[coordinate] = "A"
                elif row["ALT"][1] == "C":
                    l_ref_sequence[coordinate] = "G"
                elif row["ALT"][1] == "G":
                    l_ref_sequence[coordinate] = "C"
                indel_range = list(range(coordinate, coordinate + len(row["ALT"])))
                total_dp = row["TOTAL_DP"]
            else:
                l_ref_sequence[coordinate] = row["ALT"]

        elif row["REF_FREQ"] >= args.min_homo_prop:
            l_ref_sequence[coordinate] = row["REF"]
        elif row["ALT_FREQ"] >= 0.50:
            l_ref_sequence[coordinate] = "L"
            indel_range = list(range(coordinate, coordinate + len(row["ALT"])))
            total_dp = row["TOTAL_DP"]
        elif row["REF_FREQ"] >= 0.50:
            l_ref_sequence[coordinate] = "R"
            indel_range = list(range(coordinate, coordinate + len(row["ALT"])))
            total_dp = row["TOTAL_DP"]

    # Write alignment
    to_write += ">" + name_file + "\n" + "".join(l_ref_sequence) + "\n"

reference = header + "\n" + ref_sequence + "\n" 
out_aln.write(reference)
out_aln.write(to_write)
out_aln.close()

# IDS and sequences
aln_file = args.out_dir + "/" + args.out_name + "_aln.fasta"
ids, sequences = utils.parse_aln(aln_file)
reference = sequences[0]
id_ref = ids[0]
rest_seqs = sequences[1:]
rest_ids = ids[1:]

# Extract SNPs
snp_dict, snp_counter, aln_ref= utils.find_snps(reference, rest_seqs, rest_ids)

# Set df
# Introduce coverage
df = utils.aln2matrix(ids, aln_ref, sequences, cov_d, args)


# Plot SNPs
df_concat = utils.plot_SNPs(args, df, mutations)

# Get alingment stats
first_sample_name = args.out_dir + "/" + args.tsv_file[0].replace(".tsv", "") + "_lineage.csv"
second_sample_name = args.out_dir + "/" + args.tsv_file[1].replace(".tsv", "") + "_lineage.csv"
utils.get_alingment_stats(df, args, df_concat, first_sample_name, second_sample_name)




# # Get general stats
# utils.get_general_stats(df_nf, name_file, out_dir)

# if args.get_SNVs != "":

#     df_nf = utils.filter_tsv(args, drop_features=True, filter_NoHTZ=False)

#     # Indetify if a SNV corresponding to a certain lineage
#     df_nf["LINEAGE"] = utils.indetify_variants(df_nf, mutations)
#     df_nf.fillna(0, inplace=True)
#     df_nf = df_nf.explode("LINEAGE")
#     df_nf = df_nf[df_nf.LINEAGE == args.get_SNVs]
#     df_nf.replace(0, "", inplace=True)
#     df_nf.to_csv("%s_%s.csv" %(out_dir + "/" + name_file, args.get_SNVs), index=False, sep=",")

# # Filter HTZ SNVs
# if args.select_HTZ_SNVs:

#     df = utils.filter_tsv(args, drop_features=True, filter_NoHTZ=True)
#     df.to_csv("%s_filter.csv" %(out_dir + "/" + name_file), index=False, sep=",")

#     # Indetify if a SNV corresponding to a certain lineage
#     df["LINEAGE"] = utils.indetify_variants(df, mutations)

#     # Sort columns
#     df = df[["POS", "REF", "ALT", "REF_DP", "REF_FREQ",
#             "ALT_DP", "ALT_FREQ", "GEN", "LINEAGE",
#             "REF_CODON", "REF_AA", "ALT_CODON", "ALT_AA"]]
#     df.to_csv("%s_all.csv" %(out_dir + "/" + name_file), index=False, sep=",")

#     # Discard not identified INDELS
#     df.LINEAGE = df.LINEAGE.apply(lambda y: np.nan if len(y)==0 else y)
#     df.fillna(0, inplace=True)
#     df = df[(df.LINEAGE != 0) | (df.REF_CODON != 0)]
#     df.replace(0, "", inplace=True)
#     df.to_csv("%s_INDEL.csv" %(out_dir + "/" + name_file), index=False, sep=",")
