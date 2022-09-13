import os, multiprocessing
import pandas as pd
import numpy as np
from pandarallel import pandarallel
import matplotlib.pyplot as plt

nproc = multiprocessing.cpu_count()
pandarallel.initialize(nb_workers=nproc, verbose=0)

# check if directory exists
def check_create_dir(path):

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

# check initial arguments
def check_argmunets(args, script_dir):

    # Check mutation directory
    mut_dir = os.path.join(script_dir, "mutations")
    if not os.path.isdir(mut_dir):
        print("%s: No such file or directory" %mut_dir)
        return 1
    
    # Check reference genome
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")
    if not os.path.isfile(ref_genome):
        print("%s: No such file or directory" %ref_genome)
        return 1
    
    # Check tsv files
    if not os.path.isfile(args.tsv_file):
        print("%s: No such file or directory" %args.tsv_file)
        return 1
    
    # If all OK
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

def plot_proportions(HTZ_SNVs, name_stats_file, variant = False):

    # extract high and low percentages
    high = []
    low = []
    variant_name = os.path.basename(name_stats_file)

    if not variant:
        
        high = HTZ_SNVs[["REF_FREQ", "ALT_FREQ"]].max(axis=1).to_list()
        low = HTZ_SNVs[["REF_FREQ", "ALT_FREQ"]].min(axis=1).to_list()
    
    else:

        high = HTZ_SNVs["ALT_FREQ"].to_list()
        low = HTZ_SNVs["REF_FREQ"].to_list()
    
    # plot
    if len(high) < 5:
        width = 8
    elif len(high) < 15:
        width = 20
    else:
        width = 30

    plt.figure(figsize=(width,15))
    plt.title("HTZ positions %s" %(variant_name))
    l_pos = list(HTZ_SNVs.POS)
    for i in range(len(l_pos)):
        pos = l_pos[i]
            
        plt.bar(str(pos), high[i] + low[i], color="darkgreen")
        if i == 0:
            plt.bar(str(pos), high[i], color = "lightblue", label=variant_name)
        else:
            plt.bar(str(pos), high[i], color = "lightblue")
    plt.xlabel("Positions")
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1])
    plt.xticks(rotation = 90)
    plt.legend()
    
    # 0.5 horizontal line
    plt.axhline(y=0.5, color='r', linestyle='-')
    # mean low horizontal line
    high_mean = round(np.mean(high), 2)
    high_std = round(np.std(high), 2)
    plt.axhline(y=high_mean + high_std, color='black', linestyle='--')
    plt.axhline(y=high_mean, color='black', linestyle='-')
    plt.axhline(y=high_mean - high_std, color='black', linestyle='--')
    plt.savefig("%s.png" %(name_stats_file))

def infer_infection(name_stats_file, HOM_SNVs, mutations, name_tsv, dir_name_tsv):
    
    # out_dir_stats
    dir_name_tsv_stats = os.path.join(dir_name_tsv, "Stats")
    check_create_dir(dir_name_tsv_stats)

    # out file
    name_file = os.path.join(dir_name_tsv_stats, name_tsv)
    out_file = open("%s_inference.txt" %(name_file), "w")

    # Read stats csv
    stats_df = pd.read_csv(name_stats_file, sep=",", index_col="Lineage").T

    potential_lineages = []

    # CHECK if sample ontains > 90% lineage markers (total)
    l_marker_perc = []
    for lineage in stats_df.columns[1:]:

        if lineage in mutations:
            n_SNPs = len(mutations[lineage])

            markers_percentage = round(stats_df[lineage].values[0] / n_SNPs, 3)
            if markers_percentage > 0.9:
                l_marker_perc.append(markers_percentage)
                potential_lineages.append(lineage)
            else:
                to_write = "\n%s harbours %s" %(lineage, str(markers_percentage)) + " of markers\n"
                out_file.write(to_write)

    # CHECK if % Non_variant > 0.6 or < 0.4
    if not len(potential_lineages):

        if stats_df["Non_variant"].values[3] > 0.6 or \
                stats_df["Non_variant"].values[3] < 0.4:
            out_file.write("\nPotential co-infection between same lineage\n")
        else:
            out_file.write("\nNo co-infection\n")
        out_file.close()
        return 
    else:
        out_file.write("\nPotential lineages (SNP %)\n")
        for i in range(len(potential_lineages)):
            lineage = potential_lineages[i]
            out_file.write(lineage)
            out_file.write(" ")
            out_file.write(str(l_marker_perc[i]))
            out_file.write("\n")
    
    potential_coinfection = []
    for i in range(len(potential_lineages)):
        lineage1 = potential_lineages[i]

        for e in range(len(potential_lineages)):
            lineage2 = potential_lineages[e]
            complementary = False
            homo = False

            if e > i:
                
                # Check if % complementarity
                sum_percentages = stats_df[lineage1].values[3] + \
                    stats_df[lineage2].values[3]
                if sum_percentages > 0.93 and sum_percentages < 1.03:
                    complementary = True
                else:
                    to_write = "\n%s and %s have %s" %(lineage1, lineage2, str(round(sum_percentages, 2))) + \
                    " complementarity (lineage1 %" + " + lineage2 %)\n"
                    out_file.write(to_write)
                
                # Check if 90 % in Homozigosys:
                # List of non empty lineages in HOM_SNVs dataframe
                share_lineages = [elem for elem in HOM_SNVs["LINEAGE"].to_list() if len(elem)]
                n_pos = 0
                n_total_pos = 0
                for index in range(len(share_lineages)):
                    l_lineages = share_lineages[index]

                    if lineage1 in l_lineages:
                        n_total_pos += 1
                    if lineage1 in l_lineages and lineage2 in l_lineages:
                        n_pos += 1

                hom_percentage = n_pos / n_total_pos
                if hom_percentage > 0.9:
                    homo = True
                else:
                    to_write = "\n%s and %s have %s" %(lineage1, lineage2, str(round(hom_percentage, 2))) + \
                    " positions shared in Homo\n"
                    out_file.write(to_write)

                if complementary and homo:
                    potential_coinfection.append([lineage1, lineage2])
    
    # CHECK if % Non_variant > 0.6 or < 0.4
    if not len(potential_coinfection):

        if stats_df["Non_variant"].values[3] > 0.6 or \
                stats_df["Non_variant"].values[3] < 0.4:
            out_file.write("\nPotential co-infection between same lineage\n")
        else:
            out_file.write("\nNo co-infection\n")
        out_file.close()
        return 
    else:
        sum_total_SNP = []
        out_file.write("\nPotential co-infections pairs (Total SNPs both)\n")
        for pair in potential_coinfection:
            l1 = pair[0]
            l2 = pair[1]
            sum_total_SNP.append(stats_df[l1].values[0] + stats_df[l2].values[0])
            to_write = " - ".join(pair) + " " + str(stats_df[l1].values[0] + stats_df[l2].values[0])
            out_file.write(to_write)
            out_file.write("\n")

        
        out_file.write("\nThe best potential co-infection\n")
        best_indexes = np.argwhere(sum_total_SNP == np.amax(sum_total_SNP)).flatten().tolist()
        for index in best_indexes:
            out_file.write(" - ".join(potential_coinfection[index]))
            out_file.write("\n")
    
    out_file.close()    


def get_HTZ_stats(df, HTZ_SNVs, HOM_SNVs, mutations, name_tsv, dir_name_tsv):

    # out_dir_stats
    dir_name_tsv_stats = os.path.join(dir_name_tsv, "Stats")
    check_create_dir(dir_name_tsv_stats)

    # stats file
    name_stats_file = os.path.join(dir_name_tsv_stats, name_tsv)
    name_stats_file_ = name_stats_file[:]
    stats_file = open("%s_stats.csv" %(name_stats_file), "w")
    # header
    to_write = "Lineage,total_SNP,total_HOM,total_HTZ,mean,std,median,var,min,max\n"

    # explode
    HTZ_SNVs_explode = HTZ_SNVs.explode("LINEAGE")
    HOM_SNVs_explode = HOM_SNVs.explode("LINEAGE")

    # total no. SNPs
    total_SNPs = df.shape[0]

    # total no. htz SNPs
    total_htz = HTZ_SNVs.shape[0]
    
    # get stats HTZ proportion
    if total_htz:

        upper_HTZ_prop_l = HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()

        # statistics
        mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 3)
        std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 3)
        median_ALT_HTZ_prop = round(np.median(upper_HTZ_prop_l), 3)
        var_ALT_HTZ_prop = round(np.var(upper_HTZ_prop_l), 3)
        min_ALT_HTZ_prop = np.min(upper_HTZ_prop_l)
        max_ALT_HTZ_prop = np.max(upper_HTZ_prop_l)

        to_write += name_tsv + "," + str(total_SNPs) + "," + str(HOM_SNVs.shape[0]) + \
                "," + str(total_htz) + "," + \
                str(mean_ALT_HTZ_prop) + \
                "," + str(std_ALT_HTZ_prop) + "," + str(median_ALT_HTZ_prop) + \
                    "," + str(var_ALT_HTZ_prop) + "," + str(min_ALT_HTZ_prop) + \
                        "," + str(max_ALT_HTZ_prop) + "\n"

        # plot HTZ pos
        plot_proportions(HTZ_SNVs, name_stats_file)

    else:
        to_write += name_tsv + str(total_SNPs) + "," + str(0) + "," + str(0) + ",,,,,,\n"

    # variant stats
    for variant in mutations:

        df_variant_hom = HOM_SNVs_explode[HOM_SNVs_explode["LINEAGE"] == variant]
        df_variant_htz = HTZ_SNVs_explode[HTZ_SNVs_explode["LINEAGE"] == variant]

        # total no. of SNPs
        total_variant_SNP = df_variant_hom.shape[0] + df_variant_htz.shape[0]

        if df_variant_htz.shape[0]:
            l_variant_htz = df_variant_htz["ALT_FREQ"]

            mean_ALT_HTZ_prop = round(np.mean(l_variant_htz), 3)
            std_ALT_HTZ_prop = round(np.std(l_variant_htz), 3)
            median_ALT_HTZ_prop = round(np.median(l_variant_htz), 3)
            var_ALT_HTZ_prop = round(np.var(l_variant_htz), 3)
            min_ALT_HTZ_prop = np.min(l_variant_htz)
            max_ALT_HTZ_prop = np.max(l_variant_htz)

            to_write += variant + "," + str(total_variant_SNP) + "," + str(df_variant_hom.shape[0]) + \
                "," + str(df_variant_htz.shape[0]) + "," + \
                str(mean_ALT_HTZ_prop) + \
                "," + str(std_ALT_HTZ_prop) + "," + str(median_ALT_HTZ_prop) + \
                    "," + str(var_ALT_HTZ_prop) + "," + str(min_ALT_HTZ_prop) + \
                        "," + str(max_ALT_HTZ_prop) + "\n"
            
            name_stats_file = os.path.join(dir_name_tsv_stats, variant)
            # plot HTZ pos
            plot_proportions(df_variant_htz, name_stats_file, variant=True)

        else:
            to_write += variant + "," + str(total_variant_SNP) + "," + str(0) + "," + str(0) + ",,,,,,\n"
    
    # Not lineage
    df_non_variant_hom = HOM_SNVs_explode[HOM_SNVs_explode["LINEAGE"] == ""]
    df_non_variant_htz = HTZ_SNVs_explode[HTZ_SNVs_explode["LINEAGE"] == ""]

    # total no. of SNPs
    total_non_variant_SNP = df_non_variant_hom.shape[0] + df_non_variant_htz.shape[0]

    if df_non_variant_htz.shape[0]:

        NL_upper_HTZ_prop_l = df_non_variant_htz[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()
        
        # statistics
        mean_ALT_HTZ_prop = round(np.mean(NL_upper_HTZ_prop_l), 3)
        std_ALT_HTZ_prop = round(np.std(NL_upper_HTZ_prop_l), 3)
        median_ALT_HTZ_prop = round(np.median(NL_upper_HTZ_prop_l), 3)
        var_ALT_HTZ_prop = round(np.var(NL_upper_HTZ_prop_l), 3)
        min_ALT_HTZ_prop = np.min(NL_upper_HTZ_prop_l)
        max_ALT_HTZ_prop = np.max(NL_upper_HTZ_prop_l)

        to_write += "Non_variant" + "," + str(total_non_variant_SNP) + "," + str(df_non_variant_hom.shape[0]) + \
            "," + str(df_non_variant_htz.shape[0]) + "," + \
                str(mean_ALT_HTZ_prop) + \
                "," + str(std_ALT_HTZ_prop) + "," + str(median_ALT_HTZ_prop) + \
                    "," + str(var_ALT_HTZ_prop) + "," + str(min_ALT_HTZ_prop) + \
                        "," + str(max_ALT_HTZ_prop) + "\n"

        name_stats_file = os.path.join(dir_name_tsv_stats, "Non_variant")
        # plot HTZ pos
        plot_proportions(df_non_variant_htz, name_stats_file)
    
    else:
        to_write += "Non_variant" + str(total_non_variant_SNP) + "," + str(0) + "," + str(0) + ",,,,,,\n"

    stats_file.write(to_write)
    stats_file.close()

    return "%s_stats.csv" %(name_stats_file_)

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

    # Concatenate with SNP alignment
    df_concat = pd.concat([df_aln_SNV, l_SNPs_df.T], sort=False)

    return df_concat

def mini_compare(aln_name):

    # Parse alingment
    l_samples = {}

    # aln directory
    aln_dir = os.path.dirname(aln_name)

    # Parse alingment
    f = open(aln_name, "r")
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
    
    # COMAPARE ROWS
    samples = list(df_episodes.columns[1:])
    row = len(samples)
    matrix = np.zeros((row, row - 1))

    for i in range(len(samples)):
        column = samples[i]

        for e in range(len(samples)):
            c = samples[e]

            # Fill just one part of the matrix
            if e > i:
                matrix[e, i] = sum(
                (df_episodes[column] != df_episodes[c]) & 
            ((df_episodes[column] != "N") & (df_episodes[column] != "-")) &
             (df_episodes[c] != "N") & (df_episodes[c] != "-") &
             (df_episodes[column] != "X") & (df_episodes[c] != "X")
                )

                if e < row - 1:
                    matrix[i, e] = "nan"

    compare_df = pd.DataFrame(matrix, columns=samples[:-1], index=samples)    
    compare_df.to_csv(os.path.join(aln_dir, "episode_compare.csv"), sep=",")

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

def parse_mut(mut_dir, args):

    # Dictionary to store mutations
    d_mut = {}

    # List mutation files
    mut_files = [file for file in os.listdir(mut_dir) if file.endswith(".csv")]

    # Parse files to get mutations
    for file in mut_files:
        d_var = {}
        variant = ".".join(file.split(".")[:-1])

        # Select only include variants
        if len(args.include) and variant not in args.include:
            continue

        f = open(os.path.join(mut_dir, file), "r")
        for l in f:
            line = l.strip().replace('"', '').split(",")

            # aminoacid mutation (A34L)
            aa_mut = line[0]
            # base mutation (A3456G)
            base_mut = line[1]

            # if indel skip
            if base_mut.startswith("-"):
                continue

            # Position in genome (3456)
            pos = int(line[1][1:-1])
            # Store -> {3456: [A34L, A3456G]}
            d_var[pos] = [aa_mut, base_mut]
        f.close()

        # Store mutations related to variant
        # {BA.1: {3456: [A34L, A3456G]}}
        d_mut[variant] = d_var
    
    return d_mut

def indetify_variants(df, mutations, args):

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
    
    # If not lineages
    if len(lineage) == 0:
        lineage = [[]] * len(df.POS)

    return lineage
