import os
import pandas as pd
import numpy as np
from pandarallel import pandarallel
import dataframe_image as dfi
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

    # Same number cov and tsv files
    if len(args.tsv_file) != len(args.cov_file):
        print("Not equal number of tsv and cov files provided")
        return 1
    
    # Check mutation directory
    elif not os.path.isdir(args.mut_dir):
        print("%s: No such file or directory" %args.mut_dir)
        return 1
    
    # Check reference genome
    elif not os.path.isfile(args.ref_genome):
        return 1
    
    # Check tsv files
    for file in args.tsv_file:
        if not os.path.isfile(file):
            print("%s: No such file or directory" %file)
            return 1
        else:
            return 0
    
    # Check cov files
    for file in args.cov_file:
        if not os.path.isfile(file):
            print("%s: No such file or directory" %file)
            return 1
        else:
            return 0

def parse_ref_genome(args):

    # reference fasta
    fasta_ref = open(args.ref_genome, "r")
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


def filter_tsv(args, file, drop_features=True):

    # Read tsv
    df = pd.read_csv(file, sep=args.file_sep)

    # Drop duplicates
    df.drop_duplicates(subset=["POS"], keep="first", inplace=True)

    # Select SNPs >= min_DP_freq
    df = df[df.TOTAL_DP >= args.min_DP_freq]

    # Specify to wich gen affect the SNV
    df["GEN"] = [pos_to_gen(pos) for pos in df.POS]

    # Round to 2 ALT_FREQ
    df.ALT_FREQ = df.ALT_FREQ.round(2)

    # Set REF_FREQ
    df["REF_FREQ"] = df.parallel_apply(lambda row: 1 - row.ALT_FREQ, axis = 1)

    # Drop non relevant features
    if drop_features:
        df.drop(columns=["REGION", "REF_RV", "REF_QUAL", "ALT_RV",
                        "ALT_QUAL", "PVAL", "PASS", "GFF_FEATURE"], inplace=True)

    df = df[["POS", "REF", "ALT", "TOTAL_DP", "REF_DP", "REF_FREQ",
            "ALT_DP", "ALT_FREQ", "REF_CODON", "REF_AA", "ALT_CODON",
            "ALT_AA", "GEN"]]
    return df

def parse_coverage_file(args):

    cov_d = {}
    for file in args.cov_file:
        name = file.replace(".cov", "")
        cov_d[name] = {}
        f = open(file, "r")
        for line in f:
            l = line.strip()
            position = int(l.split("\t")[1])
            coverage = int(l.split("\t")[2])
            cov_d[name][position] = coverage
        f.close()
    return cov_d

def get_alingment_stats(df, args, df_concat, first_sample_name, second_sample_name):

    text_stats = open(args.out_dir + "/" + "alingment_stats.txt", "w")
    df_first = pd.read_csv(first_sample_name, sep=",")
    df_second = pd.read_csv(second_sample_name, sep=",")
    DF = [df_first, df_second]

    df_t = df.T
    to_write = "Total number of SNVs: "+ str(len(df.columns)) + "\n"

    # Discrepant SNVs
    id_discrepant_pos = np.where(df_t[">" + args.tsv_file[0].replace(".tsv", "")] != df_t[">" + args.tsv_file[1].replace(".tsv", "")])
    df_not_equal = df_concat.iloc[:,id_discrepant_pos[0]]
    df_not_equal_ = df_not_equal.copy()
    
    for i in range(len(DF)):
        df_ = DF[i]
        for column in df_not_equal[1:]:
            if not sum(df_["POS"] == int(column)):
                continue

            if df_not_equal[column][i + 1] == "R" or df_[df_["POS"] == int(column)].REF.values[0] == df_not_equal[column][i + 1]:
                df_not_equal[column][i + 1] = df_not_equal[column][i + 1] + " (" + str(round(df_[df_["POS"] == int(column)].REF_FREQ.values[0], 2)) + ")"

            else:
                df_not_equal[column][i + 1] = df_not_equal[column][i + 1] + " (" + str(round(df_[df_["POS"] == int(column)].ALT_FREQ.values[0], 2)) + ")"

    
    df_not_equal.to_csv(args.out_dir + "/" + args.tsv_file[0].replace(".tsv", "") + "_" + args.tsv_file[1].replace(".tsv", "") + "_not_share.csv", sep=",")
    to_write += "    Total number of unique SNVs: " + str(df_not_equal.shape[1]) + "\n"

    only_first = []
    only_second = []
    bases = ["A", "C", "G", "T"]

    for position in df_not_equal_:
        # only first
        if df_not_equal_[position][1] in bases and df_not_equal_[position][1] != df_not_equal_[position][0] and df_not_equal_[position][0] == df_not_equal_[position][2]:
            only_first.append(position)
        elif df_not_equal_[position][2] in bases and df_not_equal_[position][2] != df_not_equal_[position][0] and df_not_equal_[position][0] == df_not_equal_[position][1]:
            only_second.append(position)

    # Only first SNVs
    df_only_first = df_not_equal_[only_first]
    for i in range(len(DF)):
        if i == 0:
            df_ = DF[i]
            for column in df_only_first[1:]:

                if df_only_first[column][i + 1] == "R" or df_[df_["POS"] == int(column)].REF.values[0] == df_only_first[column][i + 1]:
                    df_only_first[column][i + 1] = df_only_first[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].REF_FREQ.values[0]) + ")"

                else:
                    df_only_first[column][i + 1] = df_only_first[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].ALT_FREQ.values[0]) + ")"

    
    df_only_first.to_csv(args.out_dir + "/" + args.tsv_file[0].replace(".tsv", "") + "_only.csv", sep=",")
    to_write += "        Number of SNVs in " + args.tsv_file[0] + " : " + str(len(only_first)) + "\n"

    # Only second SNVs
    df_only_second = df_not_equal_[only_second]
    for i in range(len(DF)):
        if i == 1:
            df_ = DF[i]
            for column in df_only_second[1:]:

                if df_only_second[column][i + 1] == "R" or df_[df_["POS"] == int(column)].REF.values[0] == df_only_second[column][i + 1]:
                    df_only_second[column][i + 1] = df_only_second[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].REF_FREQ.values[0]) + ")"

                else:
                    df_only_second[column][i + 1] = df_only_second[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].ALT_FREQ.values[0]) + ")"

    df_only_second.to_csv(args.out_dir + "/" + args.tsv_file[1].replace(".tsv", "") + "_only.csv", sep=",")
    to_write += "        Number of SNVs in " + args.tsv_file[1] + " : " + str(len(only_second)) + "\n"
    
    # Shared SNVs
    id_equal_pos = np.where(df_t[">" + args.tsv_file[0].replace(".tsv", "")] == df_t[">" + args.tsv_file[1].replace(".tsv", "")])

    df_equal = df_concat.iloc[:,id_equal_pos[0]]
    for i in range(len(DF)):
        df_ = DF[i]
        for column in df_equal[1:]:

            if df_equal[column][i + 1] == "R" or df_[df_["POS"] == int(column)].REF.values[0] == df_equal[column][i + 1]:
                df_equal[column][i + 1] = df_equal[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].REF_FREQ.values[0]) + ")"

            else:
                df_equal[column][i + 1] = df_equal[column][i + 1] + " (" + str(df_[df_["POS"] == int(column)].ALT_FREQ.values[0]) + ")"

    
    df_equal.to_csv(args.out_dir + "/" + args.tsv_file[0].replace(".tsv", "") + "_" + args.tsv_file[1].replace(".tsv", "") + "_share.csv", sep=",")
    to_write += "    Total number of shared SNVs: " + str(df_equal.shape[1]) + "\n"

    text_stats.write(to_write)
    text_stats.close()

def get_general_stats(df, name_file, out_dir, min_homo_prop):
    
    out_file = open(out_dir + "/" + name_file + "_stats.txt", "w")

    # Total number of SNVs
    to_write = "Total number of SNVs: " + str(df.shape[0]) + "\n"
    
    # HOMO SNVs
    HOM_SNVs = df[df.ALT_FREQ > min_homo_prop]
    HOM_SNVs.LINEAGE = df.LINEAGE.apply(lambda y: np.nan if len(y)==0 else y)
    HOM_SNVs.fillna(0, inplace=True)
    to_write += "    Number of HOMO SNVs: " + str(HOM_SNVs.shape[0]) + "\n"
    HOM_SNVs.to_csv(out_dir + "/" + name_file + "_HOMO.csv", index=False, sep=",")
    
    # HOMO SNVs related to a certain lineage
    HOM_SNVs_lineage = HOM_SNVs[HOM_SNVs.LINEAGE != 0]
    to_write += "        Number of HOMO SNVs related to a lineage: " + str(HOM_SNVs_lineage.shape[0]) + "\n"

    # HOMO SNVs not related to a certain lineage
    HOM_SNVs_not_lineage = HOM_SNVs[HOM_SNVs.LINEAGE == 0]
    to_write += "        Number of HOMO SNVs not related to a lineage: " + str(HOM_SNVs_not_lineage.shape[0]) + "\n"

    # HTZ SNVs
    HTZ_SNVs = df[df.ALT_FREQ < min_homo_prop]
    to_write += "    Number of HTZ SNVs: " + str(HTZ_SNVs.shape[0]) + "\n"
    HTZ_SNVs.to_csv(out_dir + "/" + name_file + "_HTZ.csv", index=False, sep=",")
    
    # INDEL HTZ SNVs
    INDEL_HTZ_SNVs = HTZ_SNVs[HTZ_SNVs.ALT.str[0].isin(["+", "-"])]
    to_write += "        Number of INDELS in HTZ SNVs: " + str(INDEL_HTZ_SNVs.shape[0]) + "\n"
    
    # NOT INDEL HTZ SNVs
    NOT_INDEL_HTZ_SNVs = HTZ_SNVs[~(HTZ_SNVs.ALT.str[0].isin(["+", "-"]))]
    to_write += "        Number of NOT INDELS in HTZ SNVs: " + str(NOT_INDEL_HTZ_SNVs.shape[0]) + "\n"

    out_file.write(to_write)
    out_file.close()

##############################################
def parse_aln(aln_file):
    sequences = []
    ids = []
    with open(aln_file, "r") as f:
        for line in f:
            if line[0] == ">":
                if len(ids):
                    sequences.append(seq)
                ids.append(line.strip())
                seq = ""
            else:
                seq += line.strip().upper()
        sequences.append(seq)
    
    return (ids, sequences)

def find_snps(reference_seq,input_seqs, rest_ids):
    non_amb = ["A","T","G","C", "L", "R"]
    # Key id, value snps
    snp_dict = {}
    # key snps value counter
    snp_counter = {}
    # key aln pos, value ref pos
    aln_ref = {}

    for index in range(len(input_seqs)):
        # sequence
        query_seq = input_seqs[index]
        ref_pos = 0
        # list to store the snps
        snps =[]

        for i in range(len(query_seq)):
            bases = [query_seq[i],reference_seq[i]]
            # consider reference position
            if reference_seq[i] in non_amb:
                ref_pos += 1
            if bases[0] != bases[1]:
                # consider only non gap positions
                if bases[0] in non_amb and bases[1] in non_amb:
                    snp = "%s%s%s" %(bases[1], str(ref_pos), bases[0])
                    if not i + 1 in aln_ref:
                        aln_ref[i + 1] = ref_pos
                    if snp not in snp_counter:
                        snp_counter[snp] = 1
                    else:
                        snp_counter[snp]+=1
                    snps.append(snp)
        snp_dict[rest_ids[index]] = snps

    return snp_dict, snp_counter, aln_ref

# create codified DF
def aln2matrix(ids, aln_ref, sequences, cov_d, args):
    df = pd.DataFrame()
    df.index = ids

    for position in aln_ref:
        ref_position = aln_ref[position]
        base_positions = []
        for i in range(len(sequences)):
            sequence = sequences[i]
            id = ids[i][1:]
            # Add letter codification
            if sequence[position - 1] == "N":
                base_positions.append("")
            elif id != "NC_045512.2" and cov_d[id][position] < args.min_DP_freq:
                base_positions.append("")
            else:
                base_positions.append(sequence[position - 1])
        df[ref_position] = base_positions
    df = df[sorted(df.columns)]
    return(df)

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

            elif cell not in ["A", "C", "G", "T"]:
                l_colors += [highlight + colours[-1] + ";"]

            elif cell != ref_cell:
                if cell == "A":
                    l_colors += [highlight + colours[0] + ";"]
                elif cell == "C":
                    l_colors += [highlight + colours[1] + ";"]
                elif cell == "G":
                    l_colors += [highlight + colours[2] + ";"]
                elif cell == "T":
                    l_colors += [highlight + colours[3] + ";"]
                else:
                    l_colors += [highlight + colours[0] + ";"]
            else:
                l_colors += [highlight + "#cbaca4" + ";"]

    return l_colors

##############################################################

def plot_SNPs(args, df, mutations):

    # Dictonary to store SNPs-location
    l_SNPs = {}

    for variant in mutations:

        l_SNPs[variant] = []
        mut_dict = mutations[variant]

        for i in range(len(df.columns)):
            pos = list(df.columns)[i]
            if pos in mut_dict:
                l_SNPs[variant].append(variant)
            else:
                l_SNPs[variant].append("")

    l_SNPs_df = pd.DataFrame()
    for variant in mutations:
        l_SNPs_df[variant] = l_SNPs[variant]
    l_SNPs_df.index = df.columns

    if len(args.get_SNVs):
        for variant in mutations:
            if variant not in args.get_SNVs:
                l_SNPs_df.drop(columns=[variant], inplace=True)

    df_concat = pd.concat([df, l_SNPs_df.T], sort=False)

    # Save to png
    dfi.export(df_concat.style.apply(color_df, axis = 0), args.out_dir + "/" + args.out_name + ".png", max_cols=-1)
    return df_concat


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
