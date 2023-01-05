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

# check fasta
def check_arg(args, script_dir):

    # Check reference genome
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")
    if not os.path.isfile(ref_genome):
        print("%s: No such file or directory" %ref_genome)
        return 1
    
    # Check fasta files
    if not os.path.isfile(args.fasta_file):
        print("%s: No such file or directory" %args.fasta_file)
        return 1
    
    # If all OK
    return 0

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
    
    # Check gff file
    gff_file = os.path.join(script_dir, "NC_045512.2.gff3")
    if not os.path.isfile(gff_file):
        print("%s: No such file or directory" %gff_file)
        return 1
    
    # Check bam files
    if not os.path.isfile(args.bamfile):
        print("%s: No such file or directory" %args.bamfile)
        return 1
    
    # Check cov files
    if not os.path.isfile(args.covfile):
        print("%s: No such file or directory" %args.covfile)
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

def parse_covfile(file):

    # Dictionary to store coverage
    cov_d = {}
    file = open(file, "r")
    for line in file:
        l = line.strip().split("\t")
        cov_d[int(l[1])] = int(l[2])
    file.close()

    return cov_d

def discard_SNP_in_DEL(df):

    # list with positions to remove
    pos_to_remove = []
    indel_pos = {}

    # previous postion features
    prev_pos = 0
    indel_len = 0
    indel_DP = 0

    for index, row in df.iterrows():
        pos = row["POS"]

        # check if DEL
        if len(row["ALT"]) > 1 and row["ALT"][0] == "-":
            prev_pos = pos
            indel_len = len(row["ALT"])
            indel_DP = row["ALT_DP"]
            for p in range(prev_pos, prev_pos + indel_len):
                if p in indel_pos:
                    if indel_pos[p] < indel_DP:
                        indel_pos[p] = indel_DP
                else:
                    indel_pos[p] = indel_DP
        
        # If not Del
        elif len(row["ALT"]) == 1:
            # if pos in DEL remove
            if pos in indel_pos:
                if indel_pos[pos] > row["ALT_DP"]:
                    pos_to_remove.append(index)
    
    # remove pos in INDEL
    df.drop(index = pos_to_remove, inplace=True)

    return df

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


def quality_control(df, args, mutations, name_tsv, dir_name_tsv):

    # out_dir_stats
    dir_name_tsv_stats = os.path.join(dir_name_tsv, "Stats")
    check_create_dir(dir_name_tsv_stats)

    # stats file
    name_stats_file = os.path.join(dir_name_tsv_stats, name_tsv)
    stats_file = open("%s_stats.csv" %(name_stats_file), "w")

    ##### HTZ PROPORTION
    # std = 0.08
    upper_std = 0.015
    # max_prop = 0.75
    # points = 0

    # fields
    fields = ["Muestra", "total_HTZ", "mean_htz_proportion",
                "std_htz_proportion", "N_SNPs_between_std", "%_SNPs_between_std", 
                "N_SNPs_higher_mean_std", "mean_dist_higher_mean_std", "std_dist_higher_mean_std",
                "N_SNPs_lower_mean_std", "mean_dist_lower_mean_std", "std_dist_lower_mean_std", 
                "Confidence_segregation"]
    row = [name_tsv]

    # get htz snps
    HTZ_SNVs = df[(df.ALT_FREQ <= args.min_HOM) & (df.ALT_FREQ >= (1 - args.min_HOM))]
    # number htz SNPs
    n_HTZ_SNPs = HTZ_SNVs.shape[0]

    if n_HTZ_SNPs:

        # select upper proportion htz
        upper_HTZ_prop_l = HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].max(axis=1).to_list()
        # Statistics
        mean_ALT_HTZ_prop = round(np.mean(upper_HTZ_prop_l), 2)
        std_ALT_HTZ_prop = round(np.std(upper_HTZ_prop_l), 2)
        N_SNPs_between_std = len([ p for p in upper_HTZ_prop_l 
                                                    if p <= mean_ALT_HTZ_prop + (std_ALT_HTZ_prop + upper_std) and
                                                    p >= mean_ALT_HTZ_prop - (std_ALT_HTZ_prop + upper_std)])
        SNPs_in_mean_limits = round(N_SNPs_between_std / len(upper_HTZ_prop_l), 2)
        
        # Info SNPs out interval mean +- std
        ## HIGHER
        N_SNPs_higher_mean_std = [p - (mean_ALT_HTZ_prop + std_ALT_HTZ_prop) for p in upper_HTZ_prop_l if p > mean_ALT_HTZ_prop + (std_ALT_HTZ_prop + upper_std)]
        mean_dist_higher_mean_std = round(np.mean(N_SNPs_higher_mean_std), 2)
        std_dist_higher_mean_std = round(np.std(N_SNPs_higher_mean_std), 2)

        ## LOWER
        N_SNPs_lower_mean_std = [(mean_ALT_HTZ_prop - std_ALT_HTZ_prop) - p for p in upper_HTZ_prop_l if p < mean_ALT_HTZ_prop - (std_ALT_HTZ_prop + upper_std)]
        mean_dist_lower_mean_std = round(np.mean(N_SNPs_lower_mean_std), 2)
        std_dist_lower_mean_std = round(np.std(N_SNPs_lower_mean_std), 2)

        row += [str(n_HTZ_SNPs), str(mean_ALT_HTZ_prop), str(std_ALT_HTZ_prop), str(N_SNPs_between_std), str(SNPs_in_mean_limits),
                str(len(N_SNPs_higher_mean_std)), str(mean_dist_higher_mean_std), str(std_dist_higher_mean_std),
                str(len(N_SNPs_lower_mean_std)), str(mean_dist_lower_mean_std), str(std_dist_lower_mean_std)]

        # If possible segregation
        if mean_ALT_HTZ_prop > 0.6:
            row += ["1"]
        else:
            row += ["0"]
        # if mean_ALT_HTZ_prop < 0.75 and std_ALT_HTZ_prop <= args.max_std_htz and SNPs_in_mean_limits >= (1 - args.SNPs_out):
        #     points += 1
    
    else:
        row += ["0"] * (len(row) - 1)

    ##### PANGOLIN
    if args.pangolin:

        # pango_fields
        fields += ["Min_pangolin", "Min_conflict", "Max_pangolin", "Max_conflict"]

        # parse pangolin files
        out_seq_dir = os.path.join(dir_name_tsv, "Sequences")
        min_file = os.path.join(out_seq_dir, name_tsv + "_1_pangolin.csv")
        min_df = pd.read_csv(min_file, sep=",")
        row += [str(min_df["lineage"].values[0]), str(min_df["conflict"].values[0])]

        max_file = os.path.join(out_seq_dir, name_tsv + "_2_pangolin.csv")
        max_df = pd.read_csv(max_file, sep=",")
        row += [str(max_df["lineage"].values[0]), str(max_df["conflict"].values[0])]

        # if min_df["conflict"].values[0] == 0 and max_df["conflict"].values[0] == 0:
        #     points += 1
    
    ##### HTZ DISTRIBUTION
    # number htz not related to lineage
    not_lineage_HTZ_SNVs = HTZ_SNVs[(HTZ_SNVs["LINEAGE"] == "") & 
                            ((HTZ_SNVs["REF_FREQ"] > args.ambiguity) |
                                ((HTZ_SNVs["ALT_FREQ"] > args.ambiguity)))]

    max_index = not_lineage_HTZ_SNVs[["ALT_FREQ", "REF_FREQ"]].idxmax(axis=1).to_list()
    
    fields += ["N_SNPs_not_lineage", "%_SNPs_not_lineage", "N_SNPs_min", "N_SNPs_max"]

    if len(max_index):

        row += [str(len(max_index)), str(round(len(max_index)/n_HTZ_SNPs, 2)),
                str(round(max_index.count("REF_FREQ") / len(max_index), 2)),
                str(round(max_index.count("ALT_FREQ") / len(max_index), 2))]
        
        # if len(max_index) > 5 and \
        #     round(max_index.count("REF_FREQ") / len(max_index), 2) > max_prop:
        #     points -= 3 

        # elif round(max_index.count("REF_FREQ") / len(max_index), 2) <= max_prop and \
        #         round(max_index.count("ALT_FREQ") / len(max_index), 2) <= max_prop:
        #     points += 1
    else:
        row += ["0", "0", "0", "0"]
    
    # fields += ["Points"]
    # row += [str(points)]

    #### WRITE FILE
    to_write = ",".join(fields) + "\n" + ",".join(row) + "\n"
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
    matrix = np.zeros((row, row))

    for i in range(len(samples)):
        column = samples[i]

        for e in range(len(samples)):
            c = samples[e]

            # Fill the matrix

            matrix[e, i] = sum(
                (df_episodes[column] != df_episodes[c]) & 
            ((df_episodes[column] != "N") & (df_episodes[column] != "-")) &
             (df_episodes[c] != "N") & (df_episodes[c] != "-")
                )

    compare_df = pd.DataFrame(matrix, columns=samples, index=samples)    
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

def indetify_variants(df, mutations):

    # List to store lineages
    lineage = []

    # Label SNVs
    for i in range(len(df.POS)):
        pos = list(df.POS)[i]

        for variant in mutations:
            #variant: BA.2

            # variant mutations
            mut_dict = mutations[variant]

            if not pos in list(mut_dict.keys()):
                if len(lineage) < df.shape[0]:
                    lineage.append([])
            
            else:
                pos_dict = mut_dict[pos]
                ref = pos_dict[1][0]
                alt = pos_dict[1][-1]

                if ref == list(df.REF)[i] and alt == list(df.ALT)[i]:
                    if len(lineage) < df.shape[0]:
                        lineage.append([variant])

                    else:
                        lineage[i] += [variant]
    
    # If not lineages
    if len(lineage) == 0:
        lineage = [[]] * len(df.POS)

    return lineage

def fasta2compare(script_dir, args, out_seq_dir, name_fasta):

    # ref genome
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")

    # Get tsv
    df = get_SNPs_fasta(ref_genome, out_seq_dir, name_fasta)
    df.to_csv(os.path.join(out_seq_dir, name_fasta) + ".tsv", sep="\t", index=False)

    # get cov
    get_cov_fasta(out_seq_dir, name_fasta)

def get_cov_fasta(out_seq_dir, name_fasta):

    # input fasta
    fasta_file = os.path.join(out_seq_dir, name_fasta) +  ".fasta"
    l_sample, sample_header, sample_sequence = parse_fasta(fasta_file)

    cov_file = open(os.path.join(out_seq_dir, name_fasta) + ".cov", "w")
    for e in range(len(l_sample)):
        pos = l_sample[e]
        
        if pos == "N" or pos == "X":
            to_write = "NC_045512.2\t%s\t0\n" %(str(e + 1))
        else:
            to_write = "NC_045512.2\t%s\t600\n" %(str(e + 1))
        cov_file.write(to_write)
    cov_file.close()

def get_SNPs_fasta(ref_genome, out_seq_dir, name_fasta):

    # parse ref genome
    l_ref_sequence, ref_header, ref_sequence = parse_fasta(ref_genome)

    # input fasta
    fasta_file = os.path.join(out_seq_dir, name_fasta) + ".fasta"
    l_sample, sample_header, sample_sequence = parse_fasta(fasta_file)

    # get SNPs
    pos = []
    ref = []
    alt = []
    for e in range(len(l_ref_sequence)):
        
        if l_sample[e] != "N" and l_ref_sequence[e] != l_sample[e] :
            
            pos.append(e + 1)
            ref.append(l_ref_sequence[e])

            # If ambiguous position set N
            if l_sample[e] == "X":
                alt.append("N")
            else:
                alt.append(l_sample[e])
    
    # Define DataFrame
    df = pd.DataFrame()
    n_rows = len(pos)
    
    # Fill SNPs
    df["REGION"] = [ref_header] * n_rows
    df["POS"] = pos
    df["REF"] = ref
    df["ALT"] = alt
    df["REF_DP"] = [0] * n_rows
    df["REF_RV"] = [0] * n_rows
    df["REF_QUAL"] = [0] * n_rows
    df["ALT_DP"] = [600] * n_rows
    df["ALT_RV"] = [73] * n_rows
    df["ALT_QUAL"] = [65] * n_rows
    df["ALT_FREQ"] = [1] * n_rows
    df["TOTAL_DP"] = [600] * n_rows
    df["PVAL"] = [0] * n_rows
    df["PASS"] = ["TRUE"] * n_rows
    df["GFF_FEATURE"] = ["NA"] * n_rows
    df["REF_CODON"] = ["NA"] * n_rows
    df["REF_AA"] = ["NA"] * n_rows
    df["ALT_CODON"] = ["NA"] * n_rows
    df["ALT_AA"] = ["NA"] * n_rows

    # return tsv as df
    return df
