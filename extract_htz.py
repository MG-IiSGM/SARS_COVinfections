# imports
import dataframe_image as dfi
from scipy import stats
import subprocess
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
                    required=True)
parser.add_argument("--mut_dir", help="directory where all mutation file are located",
                     required=False, default="../mutations")
parser.add_argument("--ref_genome", help="file with reference genome", required=False,
                    default="../COVID_ref.fasta")

parser.add_argument("--get_SNVs", help="Get SNVs from a certain lineage", 
                    default=[], required=False, action='append')
parser.add_argument("--episode", help="Extra recent episodes (samples) as fasta sequences", 
                    default=[], required=False, action='append')
parser.add_argument("--get_HTZ", help="Get HTZ and HOM positions", 
                    action="store_true")
parser.add_argument("--discard_INDEL", help="Discard INDEL positions", 
                    action="store_true")
parser.add_argument("--pangolin", help="pangolin annotation", 
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

# tsv file
tsv = args.tsv_file

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

    # HOM SNVs
    HOM_SNVs = df[(df.ALT_FREQ > args.max_prop) | (df.ALT_FREQ < args.min_prop)]
    HOM_SNVs_explode = HOM_SNVs.explode("LINEAGE")

    # Store df
    HOM_SNVs.to_csv("%s_hom.csv" %(dir_name_tsv + "/" + name_tsv),
            index=False, sep=",")
    HOM_SNVs_explode.to_csv("%s_exp_hom.csv" %(dir_name_tsv_explode + "/" + name_tsv),
            index=False, sep=",")
    
    ############################
    # Store stats HTZ
    # total no. htz SNVs
    total_htz = HTZ_SNVs.shape[0]
    
    # get mean and std HTZ proportion
    upper_HTZ_prop_l = []
    lower_HTZ_prop_l = []
    for _, row in HTZ_SNVs.iterrows():
        if row["ALT_FREQ"] > row["REF_FREQ"]:
            upper_HTZ_prop_l.append(row["ALT_FREQ"])
            lower_HTZ_prop_l.append(row["REF_FREQ"])
        else:
            upper_HTZ_prop_l.append(row["REF_FREQ"])
            lower_HTZ_prop_l.append(row["ALT_FREQ"])

    t, p_value = stats.ttest_ind(upper_HTZ_prop_l, lower_HTZ_prop_l, equal_var=False)
    mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 3)
    std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 3)
    median_ALT_HTZ_prop = round(np.median(upper_HTZ_prop_l), 3)
    var_ALT_HTZ_prop = round(np.var(upper_HTZ_prop_l), 3)
    min_ALT_HTZ_prop = np.min(upper_HTZ_prop_l)
    max_ALT_HTZ_prop = np.max(upper_HTZ_prop_l)
    
    stats = open("%s_stats.csv" %(dir_name_tsv_explode + "/" + name_tsv), "w")
    to_write = "Lineage,Total_HTZ,mean,std,median,var,min,max,p_value\n"
    to_write += name_tsv + "," + str(total_htz) + "," + str(mean_ALT_HTZ_prop) + \
                "," + str(std_ALT_HTZ_prop) + "," + str(median_ALT_HTZ_prop) + \
                    "," + str(var_ALT_HTZ_prop) + "," + str(min_ALT_HTZ_prop) + "," + str(max_ALT_HTZ_prop) + \
                    "," + str(p_value) + "\n"

    for variant in mutations:
        df_variant = HTZ_SNVs_explode[HTZ_SNVs_explode["LINEAGE"] == variant]
        if df_variant.shape[0] != 1:
            to_write += variant + "," + str(df_variant.shape[0]) + "," + str(round(df_variant["ALT_FREQ"].mean(), 3)) + "," + str(round(df_variant["ALT_FREQ"].std(), 3)) + ",,,,,\n"
        else:
            to_write += variant + ",,,,,,,,\n"
    
    # Not lineage
    df_non_variant = HTZ_SNVs_explode[HTZ_SNVs_explode["LINEAGE"] == ""]
    NL_upper_HTZ_prop_l = []
    NL_low_HTZ_prop_l = []
    for _, row in df_non_variant.iterrows():
        if row["ALT_FREQ"] > row["REF_FREQ"]:
            NL_upper_HTZ_prop_l.append(row["ALT_FREQ"])
            NL_low_HTZ_prop_l.append(row["REF_FREQ"])
        else:
            NL_upper_HTZ_prop_l.append(row["REF_FREQ"])
            NL_low_HTZ_prop_l.append(row["ALT_FREQ"])
    
    mean_ALT_HTZ_prop_no_variant = round(np.mean(NL_upper_HTZ_prop_l), 3)
    std_ALT_HTZ_prop_no_variant = round(np.std(NL_upper_HTZ_prop_l), 3)
    to_write += "No_variant" + "," + str(df_non_variant.shape[0]) + "," + str(mean_ALT_HTZ_prop_no_variant) + "," + str(std_ALT_HTZ_prop_no_variant) + "," + ",,,," + "\n"

    stats.write(to_write)
    stats.close()


    ######### ALINGMENT ########
    # out_aln_dir
    out_aln_dir = dir_name_tsv + "/ALN"
    utils.check_create_dir(out_aln_dir)

    # Parse reference sequence
    l_ref_sequence, header, ref_sequence = utils.parse_fasta(args.ref_genome)

    # List with Sample1 (<ALT_FREQ), Sample1+2 and Sample2 (>ALT_FREQ)
    genomes = [l_ref_sequence.copy(), l_ref_sequence.copy(), l_ref_sequence.copy()]
    sequences = [l_ref_sequence.copy(), l_ref_sequence.copy(), l_ref_sequence.copy()]

    for i in range(len(genomes)):
        # genome
        genome = genomes[i]
        # sequence
        sequence = sequences[i]

        # Introduce HOM positions in ref_genome
        for _, row in HOM_SNVs.iterrows():

            # position
            coordinate = row["POS"] - 1

            # If SNV
            if row["REF_FREQ"] < row["ALT_FREQ"]:
                genome[coordinate] = row["ALT"]
                sequence[coordinate] = row["ALT"]

        
        # Introduce HTZ positions in ref_genome
        for _, row in HTZ_SNVs.iterrows():
            
            # positionhom
            coordinate = row["POS"] - 1

            # If Sample1
            if i == 0:
                if row["REF_FREQ"] > row["ALT_FREQ"]:
                    genome[coordinate] = row["ALT"] + " (" + str(row["ALT_FREQ"]) + ")"
                    sequence[coordinate] = row["ALT"]
                else:
                    genome[coordinate] = row["REF"] + " (" + str(row["REF_FREQ"]) + ")"
                    sequence[coordinate] = row["REF"]
            
            # If Sample1+2
            elif i == 1:
                if row["REF_FREQ"] > row["ALT_FREQ"]:
                    genome[coordinate] = row["REF"] + "/" + row["ALT"]
                else:
                    genome[coordinate] = row["ALT"] + "/" + row["REF"]
            
            # If Sample2
            else:
                if row["REF_FREQ"] > row["ALT_FREQ"]:
                    genome[coordinate] = row["REF"] + " (" + str(row["REF_FREQ"]) + ")"
                    sequence[coordinate] = row["REF"]
                else:
                    genome[coordinate] = row["ALT"] + " (" + str(row["ALT_FREQ"]) + ")"
                    sequence[coordinate] = row["ALT"]
    
    # Store sequences to fasta files
    # out sequences dir
    out_seq_dir = dir_name_tsv + "/Sequences"
    utils.check_create_dir(out_seq_dir)

    # Sample1
    sample1 = open(out_seq_dir + "/Sample1.fasta", "w")
    to_write = ">Sample1\n" + "".join(sequences[0]) + "\n"
    sample1.write(to_write)
    sample1.close()

    if args.pangolin:
        subprocess.run(["pangolin", out_seq_dir + "/Sample1.fasta", "--outdir", out_seq_dir,
                            "--outfile", "Sample1_pangolin.csv",
                            "--max-ambig", "0.6"])

    # Sample2
    sample2 = open(out_seq_dir + "/Sample2.fasta", "w")
    to_write = ">Sample2\n" + "".join(sequences[2]) + "\n"
    sample2.write(to_write)
    sample2.close()

    if args.pangolin:
        subprocess.run(["pangolin", out_seq_dir + "/Sample2.fasta", "--outdir", out_seq_dir,
                            "--outfile", "Sample2_pangolin.csv",
                            "--max-ambig", "0.6"])
    
    ######################
    # Convert alingment to dataframe
    coordinates = list(range(1, len(genomes[0]) + 1))
    df_aln = pd.DataFrame([l_ref_sequence] + genomes, columns = coordinates,
        index=["Reference", "Sample1", "Sample1+2", "Sample2"])
    df_aln_t = df_aln.T
    df_aln_SNV = df_aln_t[(df_aln_t["Reference"] != df_aln_t["Sample1"]) | 
                    (df_aln_t["Reference"] != df_aln_t["Sample2"])].T
    
    # Set position depth
    position_dp = []
    position_l = []

    for column in df_aln_SNV.columns:
        position_l.append(column)
        position_dp.append(str(df[df["POS"] == column].TOTAL_DP.values[0]))
    
    # Concat to df_aln_SNV
    DP_df = pd.DataFrame()
    DP_df["Total_DP"] = position_dp
    DP_df.index = position_l
    df_aln_SNV = pd.concat([df_aln_SNV, DP_df.T], sort=False)

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
    df_concat.to_csv("%s_aln.csv" %(out_aln_dir + "/" + name_tsv), sep=",")

    # Color Dataframe
    dfi.export(df_concat.style.apply(utils.color_df, axis = 0),
                    "%s_aln.png" %(out_aln_dir + "/" + name_tsv),
                    max_cols=-1, )
    
    # Select only HTZ positions in Alingment
    df_concat_t = df_concat.T
    df_concat_HTZ = df_concat_t[df_concat_t["Sample1"] != df_concat_t["Sample2"]]
    df_concat_HTZ = df_concat_HTZ.T

    # Store df
    df_concat_HTZ.to_csv("%s_HTZ_aln.csv" %(out_aln_dir + "/" + name_tsv), sep=",")

    # Color Dataframe
    dfi.export(df_concat_HTZ.style.apply(utils.color_df, axis = 0),
                    "%s_HTZ_aln.png" %(out_aln_dir + "/" + name_tsv),
                    max_cols=-1, )
    ###############################

if args.episode:

    # out episode dir
    out_epi_dir = dir_name_tsv + "/Episode"
    utils.check_create_dir(out_epi_dir)

    # Construct alignment
    alin_episode = open(out_epi_dir + "/episode_aln.fasta", "w")

    # Parse reference sequence
    l_ref_sequence, header_ref, ref_sequence = utils.parse_fasta(args.ref_genome)

    # Sample1
    S1_seq_l, header_S1, S1_sequence = utils.parse_fasta(out_seq_dir + "/Sample1.fasta")

    # Sample2
    S2_seq_l, header_S2, S2_sequence = utils.parse_fasta(out_seq_dir + "/Sample2.fasta")

    # d to store other episodes
    d_episodes = {}
    for sample in args.episode:
        l_seq, header_seq, seq = utils.parse_fasta(sample)
        d_episodes[header_seq] = seq
    
    # Write header + sequence
    to_write = header_ref + "\n" + ref_sequence + "\n" + header_S1 + "\n" + \
                S1_sequence + "\n" + header_S2 + "\n" + S2_sequence + "\n"
    
    for epi in d_episodes:
        to_write += epi + "\n" + d_episodes[epi] + "\n"
    
    alin_episode.write(to_write)
    alin_episode.close()

    # Aling with mafft
    os.system("mafft " + out_epi_dir + "/episode_aln.fasta " + "> " + out_epi_dir + "/episode_aln.aln")

    # snipit
    subprocess.run(["snipit", out_epi_dir + "/episode_aln.aln", "-f", "pdf",
                            "--flip-vertical", "-o",
                            out_epi_dir + "/episode_aln"])
    
    # Compare samples
    l_samples = {}

    f = open(out_epi_dir + "/episode_aln.aln", "r")
    header = ""
    ref_sequence = ""
    for line in f:
        if line.startswith(">"):
            if header != "":
                l_samples[header] = list(ref_sequence)
            header = line.strip().split(" ")[0][1:]
            ref_sequence = ""
        else:
            ref_sequence += line.strip().upper()
    l_samples[header] = list(ref_sequence)
            
    f.close()

    df_episodes = pd.DataFrame(list(l_samples.values()), index = l_samples.keys()).T
    
    compare_row = []
    for column in df_episodes.columns[1:]:
        l = []
        for c in df_episodes.columns[1:]:
            l.append(sum((df_episodes[column] != df_episodes[c]) & 
            ((df_episodes[column] != "N") & (df_episodes[column] != "-")) &
             (df_episodes[c] != "N") & (df_episodes[c] != "-")))
        compare_row.append(l)
    
    compare_df = pd.DataFrame(compare_row, columns = list(df_episodes.columns[1:]), 
                index = list(df_episodes.columns[1:]))
    
    compare_df.to_csv(out_epi_dir + "/episode_compare.csv", sep=",")

exit(1)


to_write = ""
for file in args.tsv_file:

    # Parse reference sequence
    l_ref_sequence, header, ref_sequence = utils.parse_fasta(args)

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
