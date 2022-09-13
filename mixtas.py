# imports
import dataframe_image as dfi
import subprocess
import os
import pandas as pd
import numpy as np
import utils
import warnings
warnings.filterwarnings("ignore")

def get_lineage(args, name_tsv, mutations):

    # create out/tsv dir
    dir_name_tsv = os.path.join(args.out_dir, name_tsv)
    utils.check_create_dir(dir_name_tsv)

    # create out/tsv/explode dir
    dir_name_tsv_explode = os.path.join(dir_name_tsv, "Explode")
    utils.check_create_dir(dir_name_tsv_explode)

    # Read tsv
    df = pd.read_csv(args.tsv_file, sep=args.file_sep)

    # Drop duplicates
    df.drop_duplicates(subset=["POS"], keep="first", inplace=True)

    # Discard indel positions
    df = df[df.ALT.str.len() == 1]

    # Select SNPs with a minimum DP
    df = df[df.TOTAL_DP > args.min_DP]

    # Set gen of each SNP
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

    # Order columns
    df = df[["POS", "REF", "ALT", "TOTAL_DP", "REF_DP", "REF_FREQ",
            "ALT_DP", "ALT_FREQ", "REF_CODON", "REF_AA", "ALT_CODON",
            "ALT_AA", "GEN"]]

    # Get SNPs lineage
    df["LINEAGE"] = utils.indetify_variants(df, mutations)

    # separate lineages label 
    df_explode = df.explode("LINEAGE")

    # Change [] ""
    df.LINEAGE = df.LINEAGE.apply(lambda y: np.nan if len(y)==0 else y)
    df.fillna("", inplace=True)

    # Store df
    df.to_csv("%s_lineage.csv" %(dir_name_tsv + "/" + name_tsv),
                index=False, sep=",")
    df_explode.to_csv("%s_exp_lineage.csv" %(dir_name_tsv_explode + "/" + name_tsv),
                index=False, sep=",")
    
    return df

def get_HTZ(df, args, name_tsv, mutations):

    # directories
    dir_name_tsv = os.path.join(args.out_dir, name_tsv)
    dir_name_tsv_explode = os.path.join(dir_name_tsv, "Explode")

    # HTZ SNVs
    HTZ_SNVs = df[(df.ALT_FREQ < args.min_HOM) & (df.ALT_FREQ > (1 - args.min_HOM))]
    HTZ_SNVs_explode = HTZ_SNVs.explode("LINEAGE")

    # Store df
    HTZ_SNVs.to_csv("%s_htz.csv" %(dir_name_tsv + "/" + name_tsv),
            index=False, sep=",")
    HTZ_SNVs_explode.to_csv("%s_exp_htz.csv" %(dir_name_tsv_explode + "/" + name_tsv),
            index=False, sep=",")

    # HOM SNVs
    HOM_SNVs = df[(df.ALT_FREQ > args.min_HOM) | (df.ALT_FREQ < (1 - args.min_HOM))]
    HOM_SNVs_explode = HOM_SNVs.explode("LINEAGE")

    # Store df
    HOM_SNVs.to_csv("%s_hom.csv" %(dir_name_tsv + "/" + name_tsv),
            index=False, sep=",")
    HOM_SNVs_explode.to_csv("%s_exp_hom.csv" %(dir_name_tsv_explode + "/" + name_tsv),
            index=False, sep=",")

    # HTZ stats
    name_stats_file = utils.get_HTZ_stats(df, HTZ_SNVs, HOM_SNVs, mutations, name_tsv, dir_name_tsv)

    # Infer if co-infection
    utils.infer_infection(name_stats_file, HOM_SNVs, mutations, name_tsv, dir_name_tsv)

    return HTZ_SNVs, HOM_SNVs

def get_alingment(args, script_dir, name_tsv, HTZ_SNVs, HOM_SNVs, df, mutations):

    # directories
    dir_name_tsv = os.path.join(args.out_dir, name_tsv)

    # Parse reference sequence
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")
    l_ref_sequence, header, ref_sequence = utils.parse_fasta(ref_genome)

    # List with Sample1 (<ALT_FREQ), Sample1+2 and Sample2 (>ALT_FREQ)
    genomes = [l_ref_sequence.copy(), l_ref_sequence.copy(), l_ref_sequence.copy()]
    sequences = [l_ref_sequence.copy(), l_ref_sequence.copy(), l_ref_sequence.copy()]

    for i in range(len(genomes)):

        # genome
        genome = genomes[i]
        # sequence
        sequence = sequences[i]

        # Introduce HOM positions in ref_genome
        if HOM_SNVs.shape[0]:
            for _, row in HOM_SNVs.iterrows():

                # position
                coordinate = row["POS"] - 1

                # If SNV
                if row["REF_FREQ"] < row["ALT_FREQ"]:
                    genome[coordinate] = row["ALT"]
                    sequence[coordinate] = row["ALT"]

        
        # Introduce HTZ positions in ref_genome
        if HTZ_SNVs.shape[0]:
            for _, row in HTZ_SNVs.iterrows():
                
                # positionhom
                coordinate = row["POS"] - 1

                # If Sample1
                if i == 0:

                    if row["REF_FREQ"] > row["ALT_FREQ"]:
                        # If % high HTZ < 0.55 set X 
                        if row["REF_FREQ"] <= args.ambiguity:
                            genome[coordinate] = "X" + " (" + str(row["ALT_FREQ"]) + ")"
                            sequence[coordinate] = "X"
                        else:
                            genome[coordinate] = row["ALT"] + " (" + str(row["ALT_FREQ"]) + ")"
                            sequence[coordinate] = row["ALT"]
                    else:
                        if row["ALT_FREQ"] <= args.ambiguity:
                            genome[coordinate] = "X" + " (" + str(row["REF_FREQ"]) + ")"
                            sequence[coordinate] = "X"
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
                        # If % high HTZ < 0.55 set X 
                        if row["REF_FREQ"] <= args.ambiguity:
                            genome[coordinate] = "X" + " (" + str(row["REF_FREQ"]) + ")"
                            sequence[coordinate] = "X"
                        else:
                            genome[coordinate] = row["REF"] + " (" + str(row["REF_FREQ"]) + ")"
                            sequence[coordinate] = row["REF"]
                    else:
                        if row["ALT_FREQ"] <= args.ambiguity:
                            genome[coordinate] = "X" + " (" + str(row["ALT_FREQ"]) + ")"
                            sequence[coordinate] = "X"
                        else:
                            genome[coordinate] = row["ALT"] + " (" + str(row["ALT_FREQ"]) + ")"
                            sequence[coordinate] = row["ALT"]
    
    # Store sequences to fasta files
    # out sequences dir
    out_seq_dir = os.path.join(dir_name_tsv, "Sequences")
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
    
    aln2df(args, name_tsv, dir_name_tsv, genomes, l_ref_sequence, df, mutations)

def aln2df(args, name_tsv, dir_name_tsv, genomes, l_ref_sequence, df, mutations):

    # out_aln_dir
    out_aln_dir = dir_name_tsv + "/ALN"
    utils.check_create_dir(out_aln_dir)

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

    # include lineages to df
    df_concat = utils.include_lineages(args, df_aln_SNV, mutations)

    # Store df
    df_concat.to_csv("%s_aln.csv" %(out_aln_dir + "/" + name_tsv), sep=",")

    # Color Dataframe
    dfi.export(df_concat.style.apply(utils.color_df, axis = 0),
                    "%s_aln.png" %(out_aln_dir + "/" + name_tsv),
                    max_cols=-1)
    
    # Select only HTZ positions in Alingment
    df_concat_t = df_concat.T
    df_concat_HTZ = df_concat_t[(df_concat_t["Sample1"] != df_concat_t["Sample2"]) &
                    (df_concat_t["Sample1"] != "X") & (df_concat_t["Sample2"] != "X")]
    df_concat_HTZ = df_concat_HTZ.T

    # Store df
    df_concat_HTZ.to_csv("%s_HTZ_aln.csv" %(out_aln_dir + "/" + name_tsv), sep=",")

    # Color Dataframe
    dfi.export(df_concat_HTZ.style.apply(utils.color_df, axis = 0),
                    "%s_HTZ_aln.png" %(out_aln_dir + "/" + name_tsv),
                    max_cols=-1)

def compare_episode(args, name_tsv, script_dir):
    
    # GENERATE ALINGMENT
    # out episode dir
    dir_name_tsv = os.path.join(args.out_dir, name_tsv)
    out_epi_dir = os.path.join(dir_name_tsv, "Episode")
    out_seq_dir = os.path.join(dir_name_tsv, "Sequences")
    utils.check_create_dir(out_epi_dir)

    # Construct alignment
    alin_episode = open(os.path.join(out_epi_dir, "episode_aln.fasta"), "w")

    # Parse reference sequence
    ref_genome = os.path.join(script_dir, "COVID_ref.fasta")
    l_ref_sequence, header_ref, ref_sequence = utils.parse_fasta(ref_genome)

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

    # Align samples with mafft
    aln_name = os.path.join(out_epi_dir, "episode_aln.aln")
    try:
        subprocess.call("mafft --quiet --maxiterate 100 %s  > %s" %(os.path.join(out_epi_dir, "episode_aln.fasta"),
                        aln_name),
                        shell=True)
    except:
        print("MAFFT aligner is not installed")
        exit(1)

    # snipit
    if args.snipit:
        subprocess.run(["snipit", aln_name, "-f", "pdf",
                            "--flip-vertical", "-o",
                            aln_name.replace(".aln", "")])
    
    utils.mini_compare(aln_name)
