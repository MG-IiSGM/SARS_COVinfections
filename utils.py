import os
import pandas as pd
import numpy as np
from scipy import stats
from pandarallel import pandarallel
import matplotlib.pyplot as plt
import multiprocessing
nproc = multiprocessing.cpu_count()
pandarallel.initialize(nb_workers=nproc, verbose=0)

# check if directory exists
def check_create_dir(path):

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

# check initial arguments
def check_argmunets(args):

    # Check mutation directory
    if not os.path.isdir(args.mut_dir):
        print("%s: No such file or directory" %args.mut_dir)
        return 1
    
    # Check reference genome
    elif not os.path.isfile(args.ref_genome):
        print("%s: No such file or directory" %args.ref_genome)
        return 1
    
    # Check tsv files
    elif not os.path.isfile(args.tsv_file):
        print("%s: No such file or directory" %args.tsv_file)
        return 1
    
    # If all OK
    else:
        return 0

# read fasta function
def parse_fasta(file):

    # reference fasta
    fasta_ref = open(file, "r")
    ref_sequence = ""
    header = ""
    for line in fasta_ref:
        if not line.startswith(">"):
            ref_sequence += line.strip()
        else:
            header = line.strip().split(" ")[0]

    # reference genome
    l_ref_sequence = list(ref_sequence)
    fasta_ref.close()

    return l_ref_sequence, header, ref_sequence

# convert positions to gen
def pos_to_gen(pos):

    if pos > 265 and pos < 21556:
        return "ORF1ab" # "ORF1ab":[266, 21555]
    elif pos > 21562 and pos < 25385:
        return "S" # "S":[21563, 25384]
    elif pos > 25392 and pos < 26221:
        return "ORF3a" # "ORF3a":[25393, 26220]
    elif pos > 26244 and pos < 26473:
        return "E" #  "E":[26245, 26472]
    elif pos > 26522 and pos < 27192:
        return "M" # "M":[26523, 27191]
    elif pos > 27201 and pos < 27388:
        return "ORF6" # "ORF6":[27202, 27387]
    elif pos > 27393 and pos < 27760:
        return "ORF7a" # "ORF7a":[27394, 27759]
    elif pos > 27755 and pos < 27888:
        return "ORF7b" # "ORF7b":[27756, 27887]
    elif pos > 27893 and pos < 28260:
        return "ORF8" # "ORF8":[27894, 28259]
    elif pos > 28273 and pos < 29534:
        return "N" # "N":[28274, 29533]
    elif pos > 29557 and pos < 29675:
        return "ORF10" # "ORF10":[29558, 29674]

def plot_proportions(HTZ_SNVs, name_stats_file):

    # extract high and low percentages
    high = []
    low = []
    for r, a in zip(HTZ_SNVs.REF_FREQ, HTZ_SNVs.ALT_FREQ):
        if r > a:
            high.append(r)
            low.append(a)
        else:
            high.append(a)
            low.append(r)
    
    # plot
    plt.figure(figsize=(15, 12))
    plt.title("HTZ positions")
    l_pos = list(HTZ_SNVs.POS)
    for i in range(len(l_pos)):
        pos = l_pos[i]
        plt.bar(str(pos), high[i] + low[i], color="green")
        plt.bar(str(pos), high[i], color = "lightblue")
    plt.xlabel("Positions")
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1])
    plt.xticks(rotation = 90)
    
    # 0.5 horizontal line
    plt.axhline(y=0.5, color='r', linestyle='-')
    # mean low horizontal line
    high_mean = round(np.mean(high), 2)
    high_std = round(np.std(high), 2)
    plt.axhline(y=high_mean + high_std, color='black', linestyle='--')
    plt.axhline(y=high_mean, color='black', linestyle='-')
    plt.axhline(y=high_mean - high_std, color='black', linestyle='--')
    plt.savefig("%s.png" %(name_stats_file))


def get_HTZ_stats(HTZ_SNVs, HTZ_SNVs_explode, mutations, name_tsv, dir_name_tsv):

    # out_dir_stats
    dir_name_tsv_stats = os.path.join(dir_name_tsv, "Stats")
    check_create_dir(dir_name_tsv_stats)

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

    # statistics
    t, p_value = stats.ttest_ind(upper_HTZ_prop_l, lower_HTZ_prop_l, equal_var=False)
    mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 3)
    std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 3)
    median_ALT_HTZ_prop = round(np.median(upper_HTZ_prop_l), 3)
    var_ALT_HTZ_prop = round(np.var(upper_HTZ_prop_l), 3)
    try:
        min_ALT_HTZ_prop = np.min(upper_HTZ_prop_l)
    except:
        min_ALT_HTZ_prop = 0
    try:
        max_ALT_HTZ_prop = np.max(upper_HTZ_prop_l)
    except:
        max_ALT_HTZ_prop = 0
    
    name_stats_file = os.path.join(dir_name_tsv_stats, name_tsv)
    stats_file = open("%s_stats.csv" %(name_stats_file), "w")
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

    stats_file.write(to_write)
    stats_file.close()

    # plot HTZ pos
    plot_proportions(HTZ_SNVs, name_stats_file)

def include_lineages(args, df_aln_SNV, mutations):
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

    # Select only lineages specified in get_SNPs
    if len(args.get_SNP):
        for variant in mutations:
            if variant not in args.get_SNP:
                l_SNPs_df.drop(columns=[variant], inplace=True)

    # Concatenate with SNP alignment
    df_concat = pd.concat([df_aln_SNV, l_SNPs_df.T], sort=False)

    return df_concat

def mini_compare(out_epi_dir):

    # Compare samples
    l_samples = {}

    f = open(os.path.join(out_epi_dir, "episode_aln.aln"), "r")
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
            l.append(sum(
                (df_episodes[column] != df_episodes[c]) & 
            ((df_episodes[column] != "N") & (df_episodes[column] != "-")) &
             (df_episodes[c] != "N") & (df_episodes[c] != "-") &
             (df_episodes[column] != "X") & (df_episodes[c] != "X")
                ))
        compare_row.append(l)
    
    compare_df = pd.DataFrame(compare_row, columns = list(df_episodes.columns[1:]), 
                index = list(df_episodes.columns[1:]))
    
    compare_df.to_csv(os.path.join(out_epi_dir, "episode_compare.csv"), sep=",")

def color_df(row):

    # Colours
    colours = ["steelblue","indianred","darkseagreen","skyblue", "yellow", "white"]
    highlight = 'background-color: '
    l_colors = [highlight + "#cbaca4" + ";"]

    # row
    ref_cell = row[0]
    for i in range(len(row)):
        cell = row[i]
        
        # Not ref
        if i:
            if "/" in cell:
                l_colors += [highlight + colours[4] + ";"]

            # Only colour SNVs 
            elif i > 4 or (len(cell) and cell[0] not in ["A", "C", "G", "T"]):
                l_colors += [highlight + colours[-1] + ";"]

            elif len(cell) and cell[0] != ref_cell:
                if cell[0] == "A":
                    l_colors += [highlight + colours[0] + ";"]
                elif cell[0] == "C":
                    l_colors += [highlight + colours[1] + ";"]
                elif cell[0] == "G":
                    l_colors += [highlight + colours[2] + ";"]
                elif cell[0] == "T":
                    l_colors += [highlight + colours[3] + ";"]
                else:
                    l_colors += [highlight + colours[0] + ";"]
            else:
                l_colors += [highlight + "#cbaca4" + ";"]

    return l_colors

def parse_mut(args):

    # Dictionary to store mutations
    d_mut = {}

    # List mutation files
    mut_files = [file for file in os.listdir(args.mut_dir) if file.endswith(".csv")]

    # Parse files to get mutations
    for file in mut_files:
        d_var = {}
        variant = file.replace(".csv", "")

        f = open(os.path.join(args.mut_dir, file), "r")
        for l in f:
            line = l.strip().replace('"', '').split(",")

            # aminoacid mutation (A34L)
            aa_mut = line[0]
            # base mutation (A3456G)
            base_mut = line[1]
            # Position in genome (3456)
            pos = int(line[1][1:-1])
            # Store -> {3456: [A34L, A3456G]}
            d_var[pos] = [aa_mut, base_mut]
        f.close()

        # Store mutations related to variant
        # {BA.1: {3456: [A34L, A3456G]}}
        d_mut[variant] = d_var
    
    return d_mut

def indetify_variants(df, mutations):

    # List to store lineages
    lineage = []

    for variant in mutations:
        #variant: BA.2

        # variant mutations
        mut_dict = mutations[variant]

        # Label SNVs
        for i in range(len(df.POS)):
            pos = list(df.POS)[i]

            if not pos in list(mut_dict.keys()) and len(lineage) < df.shape[0]:
                lineage.append([])

            elif pos in list(mut_dict.keys()) and len(lineage) < df.shape[0]:
                lineage.append([variant])

            elif pos in list(mut_dict.keys()):
                lineage[i] += [variant]

    return lineage
