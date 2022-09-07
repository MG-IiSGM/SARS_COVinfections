# imports
import dataframe_image as dfi
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
parser.add_argument("--get_HTZ", help="Get HTZ and HOM positions", 
                    action="store_true")
parser.add_argument("--discard_INDEL", help="Discard INDEL positions", 
                    action="store_true")
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--out_name", help="alignment out name", default="alingment")
parser.add_argument("--file_sep", help="File separator", default="\t")
parser.add_argument("--min_DP", help="minimum frequency (depth) to accept a SNV", default=15)
parser.add_argument("--min_prop", help="minimum proportion for htz", default=0.15)
parser.add_argument("--max_prop", help="minimum proportion for htz", default=0.85)

# check arguments
args = parser.parse_args()

if utils.check_argmunets(args):
    exit(1)
else:

    # Parse mutations.
    mutations = utils.parse_mut(args)

    # Check out dir
    utils.check_create_dir(args.out_dir)

    # parse coverage to dictionary
    cov_d = utils.parse_coverage_file(args)

for tsv in args.tsv_file:

    # create out/tsv dir
    name_tsv = os.path.basename(tsv).replace(".tsv", "")
    dir_name_tsv = args.out_dir + "/" + name_tsv
    dir_name_tsv_explode = dir_name_tsv + "/" + "explode"
    utils.check_create_dir(dir_name_tsv)
    utils.check_create_dir(dir_name_tsv_explode)

    # Read tsv
    df = pd.read_csv(tsv, sep=args.file_sep)

    # Drop duplicates
    df.drop_duplicates(subset=["POS"], keep="first", inplace=True)

    if args.discard_INDEL:
        df = df[df.ALT.str.len() == 1]

    # Select SNVs with a minimum DP
    df = df[df.TOTAL_DP > args.min_DP]

    # Set gen of each SNV
    df["GEN"] = [utils.pos_to_gen(pos) for pos in df.POS]

    # Round to 3 ALT_FREQ
    df.ALT_FREQ = df.ALT_FREQ.round(3)

    # Set REF_FREQ
    df["REF_FREQ"] = df.parallel_apply(lambda row: 1 - row.ALT_FREQ, axis = 1)
    df.REF_FREQ = df.REF_FREQ.round(3)

    # Drop non relevant features
    df.drop(columns=["REGION", "REF_RV", "REF_QUAL", "ALT_RV",
                     "ALT_QUAL", "PVAL", "PASS", "GFF_FEATURE"], 
                     inplace=True)

    df = df[["POS", "REF", "ALT", "TOTAL_DP", "REF_DP", "REF_FREQ",
            "ALT_DP", "ALT_FREQ", "REF_CODON", "REF_AA", "ALT_CODON",
            "ALT_AA", "GEN"]]

    # Label SNVs to lineage
    df["LINEAGE"] = utils.indetify_variants(df, mutations)

    # separate lineage label 
    df_explode = df.explode("LINEAGE")

    # Change [] ""
    df.LINEAGE = df.LINEAGE.apply(lambda y: np.nan if len(y)==0 else y)
    df.fillna("", inplace=True)

    # Store df
    df.to_csv("%s_lineage.csv" %(dir_name_tsv + "/" + name_tsv),
             index=False, sep=",")
    df_explode.to_csv("%s_exp_lineage.csv" %(dir_name_tsv_explode + "/" + name_tsv),
             index=False, sep=",")

    if args.get_HTZ:

        # HTZ SNVs
        HTZ_SNVs = df[(df.ALT_FREQ < args.max_prop) & (df.ALT_FREQ > args.min_prop)]
        HTZ_SNVs_explode = HTZ_SNVs.explode("LINEAGE")

        # Store df
        HTZ_SNVs.to_csv("%s_htz.csv" %(dir_name_tsv + "/" + name_tsv),
                index=False, sep=",")
        HTZ_SNVs_explode.to_csv("%s_exp_htz.csv" %(dir_name_tsv_explode + "/" + name_tsv),
             index=False, sep=",")

        # total no. htz SNVs
        total_htz = HTZ_SNVs.shape[0]
        
        # get mean and std HTZ proportion
        mean_ALT_HTZ_prop = round(HTZ_SNVs["ALT_FREQ"].mean(), 3)
        std_ALT_HTZ_prop = round(HTZ_SNVs["ALT_FREQ"].std(), 3)

        # HOM SNVs
        HOM_SNVs = df[(df.ALT_FREQ > args.max_prop) | (df.ALT_FREQ < args.min_prop)]
        HOM_SNVs_explode = HOM_SNVs.explode("LINEAGE")

        # Store df
        HOM_SNVs.to_csv("%s_hom.csv" %(dir_name_tsv + "/" + name_tsv),
                index=False, sep=",")
        HOM_SNVs_explode.to_csv("%s_exp_hom.csv" %(dir_name_tsv_explode + "/" + name_tsv),
             index=False, sep=",")
        
        # Store stats HTZ
        stats = open("%s_stats.csv" %(dir_name_tsv_explode + "/" + name_tsv), "w")
        to_write = "Lineage,Total_HTZ,mean,std\n"
        to_write += "ALL" + "," + str(total_htz) + "," + str(mean_ALT_HTZ_prop) + "," + str(std_ALT_HTZ_prop) + "\n"

        for variant in mutations:
            df_variant = HTZ_SNVs_explode[HTZ_SNVs_explode["LINEAGE"] == variant]
            if df_variant.shape[0] != 1:
                to_write += variant + "," + str(df_variant.shape[0]) + "," + str(round(df_variant["ALT_FREQ"].mean(), 3)) + "," + str(round(df_variant["ALT_FREQ"].std(), 3)) + "\n"
            else:
                to_write += variant + ",,,\n"

        stats.write(to_write)
        stats.close()


        ######### ALINGMENT ########
        # Parse reference sequence
        l_ref_sequence, header, ref_sequence = utils.parse_ref_genome(args)

        # List with Sample1 (<ALT_FREQ), Sample1+2 and Sample2 (>ALT_FREQ)
        genomes = [l_ref_sequence.copy(), l_ref_sequence.copy(), l_ref_sequence.copy()]

        for i in range(len(genomes)):
            # genome
            genome = genomes[i]

            # Introduce HOM positions in ref_genome
            for _, row in HOM_SNVs.iterrows():

                # position
                coordinate = row["POS"] - 1

                # If SNV
                if row["REF_FREQ"] < row["ALT_FREQ"]:
                    genome[coordinate] = row["ALT"]
            
            # Introduce HTZ positions in ref_genome
            for _, row in HTZ_SNVs.iterrows():
                
                # positionhom
                coordinate = row["POS"] - 1

                # If Sample1
                if i == 0:
                    if row["REF_FREQ"] > row["ALT_FREQ"]:
                        genome[coordinate] = row["ALT"]
                    else:
                        genome[coordinate] = row["REF"]
                
                # If Sample1+2
                elif i == 1:
                    if row["REF_FREQ"] > row["ALT_FREQ"]:
                        genome[coordinate] = row["REF"] + "/" + row["ALT"]
                    else:
                        genome[coordinate] = row["ALT"] + "/" + row["REF"]
                
                # If Sample2
                else:
                    if row["REF_FREQ"] > row["ALT_FREQ"]:
                        genome[coordinate] = row["REF"]
                    else:
                        genome[coordinate] = row["ALT"]
        
        # Convert alingment to dataframe
        coordinates = list(range(1, len(genomes[0]) + 1))
        df_aln = pd.DataFrame([l_ref_sequence] + genomes, columns = coordinates,
         index=["Reference", "Sample1", "Sample1+2", "Sample2"])
        df_aln_t = df_aln.T
        df_aln_SNV = df_aln_t[(df_aln_t["Reference"] != df_aln_t["Sample1"]) | 
                        (df_aln_t["Reference"] != df_aln_t["Sample2"])].T
        
        # Set lineage
        # Dictonary to store SNPs-location
        l_SNPs = {}

        for variant in mutations:
            # BA.2

            # List to store mutations asociated to the lineage
            l_SNPs[variant] = []

            # Positions associated with the lineage
            mut_dict = mutations[variant]

            for i in range(len(df_aln_SNV.columns)):

                pos = list(df_aln_SNV.columns)[i]

                if pos in mut_dict:
                    l_SNPs[variant].append(variant)
                else:
                    l_SNPs[variant].append("")

        # Create pandas Dataframe with positions
        l_SNPs_df = pd.DataFrame()
        for variant in mutations:
            l_SNPs_df[variant] = l_SNPs[variant]
        l_SNPs_df.index = df_aln_SNV.columns

        # Select only lineages specified in get_SNVs
        if len(args.get_SNVs):
            for variant in mutations:
                if variant not in args.get_SNVs:
                    l_SNPs_df.drop(columns=[variant], inplace=True)

        # Concatenate with SNV alignment
        df_concat = pd.concat([df_aln_SNV, l_SNPs_df.T], sort=False)

        # Store df
        df_concat.to_csv("%s_aln.csv" %(dir_name_tsv + "/" + name_tsv), sep=",")


        # Color Dataframe
        dfi.export(df_concat.style.apply(utils.color_df, axis = 0),
                     "%s_aln.png" %(dir_name_tsv + "/" + name_tsv),
                     max_cols=-1, )

exit(1)


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
