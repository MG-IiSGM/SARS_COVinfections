import os, multiprocessing
import pandas as pd
import numpy as np
from pandarallel import pandarallel
import matplotlib.pyplot as plt
import re, subprocess, gzip

nproc = multiprocessing.cpu_count()
pandarallel.initialize(nb_workers=nproc, verbose=0)

def extract_read_list(input_dir):
    """
    Search files in a directory sort by name and extract comon name of R1 and R2
    with extract_sample() function
    190615 - Limit only parent folder, not subdirectories
    """
    input_dir = os.path.abspath(input_dir)
    all_fasta = []
    r1_list = []
    r2_list = []
    for root, _, files in os.walk(input_dir):
        if root == input_dir:  # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                is_fasta = re.match(r'.*\.f(ast)*[aq](\.gz)*', filename)
                if is_fasta:
                    all_fasta.append(filename)
    all_fasta = sorted(all_fasta)
    if len(all_fasta) % 2 == 0:
        for index, fasta_file in enumerate(all_fasta):
            if index % 2 == 0:
                r1_list.append(fasta_file)
            elif index % 1 == 0:
                r2_list.append(fasta_file)
 
    r1_list = sorted(r1_list)
    r2_list = sorted(r2_list)

    return r1_list, r2_list

def extract_sample(R1_file, R2_file):
    """
    Extract sample from R1, R2 files.
    """
    basename_R1 = os.path.basename(R1_file)
    basename_R2 = os.path.basename(R2_file)

    sample_name_R = os.path.commonprefix([basename_R1, basename_R2])

    long_suffix = re.search('_S.*', sample_name_R)
    short_suffix = re.search('_R.*', sample_name_R)
    bar_suffix = re.search('_$', sample_name_R)
    dot_suffix = re.search('.R$', sample_name_R)

    if long_suffix:
        match = long_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif short_suffix:
        match = short_suffix.group()
        sample_name = sample_name_R.split(match)[0]
    elif bar_suffix:
        match = bar_suffix.group()
        sample_name = sample_name_R.rstrip("_")
    elif dot_suffix:
        match = dot_suffix.group()
        sample_name = sample_name_R.rstrip(".R")
    else:
        sample_name = sample_name_R

    return sample_name

def fastqc_quality(r1, r2, output_dir, threads=8):
    check_create_dir(output_dir)

    cmd = ['fastqc', r1, r2, '-o', output_dir, '--threads', str(threads)]

    subprocess.call(cmd, shell=True)

def fastp_trimming(r1, r2, sample, output_dir, threads=6, min_qual=20, window_size=10, min_len=35):
    check_create_dir(output_dir)

    output_trimmed_r1 = os.path.join(
        output_dir, sample + ".trimmed_R1.fastq.gz")
    output_trimmed_r2 = os.path.join(
        output_dir, sample + ".trimmed_R2.fastq.gz")

    html_dir = os.path.join(output_dir, 'html')
    json_dir = os.path.join(output_dir, 'json')

    check_create_dir(html_dir)
    check_create_dir(json_dir)

    html_file = os.path.join(html_dir, sample + '_fastp.html')
    json_file = os.path.join(json_dir, sample + '_fastp.json')

    cmd = ['fastp',
           '--in1', r1,
           '--in2', r2,
           '--out1', output_trimmed_r1,
           '--out2', output_trimmed_r2,
           '--detect_adapter_for_pe',
           '--cut_tail',
           '--cut_window_size', str(window_size),
           '--cut_mean_quality', str(min_qual),
           '--length_required', str(min_len),
           '--json', json_file,
           '--html', html_file,
           '--thread', str(threads)]

    subprocess.call(cmd, shell=True)

def bwa_mapping(r1, r2, reference, sample, output_dir, threads=8):
    """
    #Store output in a file when it is outputted in stdout
    https://stackoverflow.com/questions/4965159/how-to-redirect-output-with-subprocess-in-python
    """
    r1 = os.path.abspath(r1)
    r2 = os.path.abspath(r2)
    reference = os.path.abspath(reference)

    sample_name = sample + ".sam"
    output_file = os.path.join(output_dir, sample_name)

    check_create_dir(output_dir)
    
    cmd_index = ["bwa", "index", reference]
    subprocess.call(cmd_index, shell=True)
    
    cmd_map = ["bwa", "mem", "-Y", "-M", "-t", str(threads), "-o", output_file, reference, r1, r2]
    subprocess.call(cmd_map, shell=True)

def sam_to_index_bam(sample, output_dir, r1, threads):
    # input_sam_path = os.path.abspath(input_sam)
    # if output_bam == "inputdir":
    #     output_bam = os.path.dirname(input_sam_path)
    # else:
    #     output_bam = output_bam

    sample_name = sample + ".sam"
    input_sam_path = os.path.join(output_dir, sample_name)

    input_name = (".").join(os.path.basename(input_sam_path).split(".")[:-1])

    output_bam_name = input_name + ".bam"
    output_bam_path = os.path.join(output_dir, output_bam_name)

    output_sorted_name = input_name + ".sorted.bam"
    output_sorted_path = os.path.join(output_dir, output_sorted_name)

    output_bg_sorted_name = input_name + ".rg.sorted.bam"
    output_bg_sorted_path = os.path.join(output_dir, output_bg_sorted_name)

    cmd_view = ["samtools", "view", "-Sb", input_sam_path, "--threads", str(threads), "-o", output_bam_path,]
    subprocess.call(cmd_view, shell=True)

    check_remove_file(input_sam_path)
    
    cmd_sort = ["samtools", "sort", output_bam_path, "-o", output_sorted_path]
    subprocess.call(cmd_sort, shell=True)

    check_remove_file(output_bam_path)

    add_SG(sample, output_sorted_path, output_bg_sorted_path, r1)

    check_remove_file(output_sorted_path)

def add_SG(sample, input_bam, output_bg_sorted, r1):
    """
    @MN00227:45:000H255J3:1:11102:21214:1110 1:N:0:18
    @NS500454:48:HKG57BGXX:1:11101:17089:1032 2:N:0:TCCTGAGC+TCTTACGC
    @NS500454:27:HJJ32BGXX:1:11101:12392:1099 1:N:0:2
    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:
    <is filtered>:<control number>:<sample number | barcode1'+barcode2'>
    ID = Read group identifier {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} 
    PU = Platform Unit #optional
    SM = Sample
    PL = Platform/technology used to produce the read (ILLUMINA, SOLID, LS454, HELICOS and PACBIO)
    LB = DNA preparation library identifier
    """

    with gzip.open(r1) as f:
        first_line = f.readline().strip().decode()
    #print(first_line)
    first_line_list = first_line.split(":")
    if len(first_line_list) > 4:
        rg_id = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
        rg_pu = ".".join([first_line_list[2],first_line_list[3],first_line_list[-1]])
    else:
        rg_id = first_line_list[0]
        rg_pu = first_line_list[0]
    rg_sm = sample
    rg_pl = "ILLUMINA"
    rg_lb = "lib_" + sample

    rg_id_param = "RGID=" + rg_id
    rg_pu_param = "RGPU=" + rg_pu
    rg_sm_param = "RGSM=" + rg_sm
    rg_pl_param = "RGPL=" + rg_pl
    rg_lb_param = "RGLB=" + rg_lb

    input_param = "INPUT=" + input_bam
    output_param = "OUTPUT=" + output_bg_sorted


    cmd = ["picard", "AddOrReplaceReadGroups", 
    input_param, output_param, rg_id_param, rg_lb_param, rg_pl_param, rg_pu_param, rg_sm_param,
    "SORT_ORDER=coordinate"]
    subprocess.call(cmd, shell=True)

def picard_markdup(input_bam):
    
    input_bam = os.path.abspath(input_bam)
    
    path_file_name = input_bam.split(".")[0]
    file_name = input_bam.split("/")[-1]
    output_markdup = path_file_name + ".rg.markdup.bam"
    output_markdup_sorted = path_file_name + ".rg.markdup.sorted.bam"

    output_dir = ('/').join(input_bam.split('/')[0:-1])
    stat_output_dir = os.path.join(output_dir, "Stats")
    stat_output_file = file_name + ".markdup.metrics.txt"
    stat_output_full = os.path.join(stat_output_dir, stat_output_file)

    check_create_dir(stat_output_dir)

    cmd_markdup = ["picard", "MarkDuplicates", "-I", input_bam, "-O", output_markdup, "-M", stat_output_full]
    subprocess.call(cmd_markdup, shell=True)
    
    #samtools sort: samtools sort $output_dir/$sample".sorted.bam" -o $output_dir/$sample".sorted.bam"
    cmd_sort = ["samtools", "sort", output_markdup, "-o", output_markdup_sorted]
    subprocess.call(cmd_sort, shell=True)

    check_remove_file(input_bam)
    check_remove_file(output_markdup)

def ivar_trim(input_bam, primers_file, sample, min_length=30, min_quality=20, sliding_window_width=4):
    """
    Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]
        Input Options    Description
           -i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           -b    (Required) BED file with primer sequences and positions
           -m    Minimum length of read to retain after trimming (Default: 30)
           -q    Minimum quality threshold for sliding window to pass (Default: 20)
           -s    Width of sliding window (Default: 4)
           -e    Include reads with no primers. By default, reads with no primers are excluded
        Output Options   Description
           -p    (Required) Prefix for the output BAM file
    """
    
    input_bam = os.path.abspath(input_bam)
    input_bai = input_bam + ".bai"
    primers_file = os.path.abspath(primers_file)

    prefix = input_bam.split('.')[0] + ".rg.markdup.trimmed"
    output_trimmed_bam = prefix + ".bam"
    output_trimmed_sorted_bam = input_bam.split('.')[0] + ".rg.markdup.trimmed.sorted.bam"
    
    cmd = ["ivar", "trim", "-i", input_bam, "-b", primers_file, "-p", prefix, "-m", str(min_length), "-q", str(min_quality), "-s", str(sliding_window_width), "-e"]
    subprocess.call(cmd, shell=True)

    check_remove_file(input_bam)

    cmd_sort = ["samtools", "sort", output_trimmed_bam, "-o", output_trimmed_sorted_bam]
    subprocess.call(cmd_sort, shell=True)

    check_remove_file(output_trimmed_bam)

    cmd_index = ["samtools", "index", output_trimmed_sorted_bam]
    subprocess.call(cmd_index, shell=True)

    check_remove_file(input_bai)

def ivar_variants(reference, input_bam, output_variant, sample, min_quality=20, min_frequency_threshold=0.8, min_depth=20):
    """
    Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference [<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] [-r <reference-fasta>] [-g GFF file]
        Note : samtools mpileup output must be piped into ivar variants
        Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
           -m    Minimum read depth to call variants (Default: 0)
           -r    Reference file used for alignment. This is used to translate the nucleotide sequences and identify intra host single nucleotide variants
           -g    A GFF file in the GFF3 format can be supplied to specify coordinates of open reading frames (ORFs). In absence of GFF file, amino acid translation will not be done.
        Output Options   Description
           -p    (Required) Prefix for the output tsv variant file
    """
    ivar_folder = os.path.join(output_variant, 'ivar_raw')
    check_create_dir(ivar_folder)
    prefix = ivar_folder + '/' + sample

    input = {'reference' : reference,
            'input_bam': input_bam,
            'prefix' : prefix,
            'min_quality': str(min_quality),
            'min_frequency_threshold': str(min_frequency_threshold),
            'min_depth': str(min_depth)}


    cmd = "samtools mpileup -aa -A -d 0 -B -Q 0 --reference {reference} {input_bam} | \
        ivar variants -p {prefix} -q {min_quality} -t {min_frequency_threshold} -m {min_depth} -r {reference}".format(**input)

    subprocess.call(cmd, Shell=True)

def ivar_consensus(input_bam, output_consensus, sample, min_quality=20, min_frequency_threshold=0.8, min_depth=20, uncovered_character='N'):
    """
    ivar consensus
        Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 
        Note : samtools mpileup output must be piped into ivar consensus
        Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call consensus. (Default: 0)
                 Frequently used thresholds | Description
                 ---------------------------|------------
                                          0 | Majority or most common base
                                        0.2 | Bases that make up atleast 20% of the depth at a position
                                        0.5 | Strict or bases that make up atleast 50% of the depth at a position
                                        0.9 | Strict or bases that make up atleast 90% of the depth at a position
                                          1 | Identical or bases that make up 100% of the depth at a position. Will have highest ambiguities
           -m    Minimum depth to call consensus(Default: 10)
           -k    If '-k' flag is added, regions with depth less than minimum depth will not be added to the consensus sequence. Using '-k' will override any option specified using -n 
           -n    (N/-) Character to print in regions with less than minimum coverage(Default: N)
        Output Options   Description
           -p    (Required) Prefix for the output fasta file and quality file
    """

    prefix = output_consensus + '/' + sample

    input = {'input_bam': input_bam,
            'prefix' : prefix,
            'min_quality': str(min_quality),
            'min_frequency_threshold': str(min_frequency_threshold),
            'min_depth': str(min_depth),
            'uncovered_character': uncovered_character}

    cmd = "samtools mpileup -aa -A -d 0 -B -Q 0  {input_bam} | \
        ivar consensus -p {prefix} -q {min_quality} -t {min_frequency_threshold} -m {min_depth} -n {uncovered_character}".format(**input)

    subprocess.call(cmd, Shell=True)

def create_bamstat(input_bam, output_dir, sample, threads=8):
    output_file = os.path.join(output_dir, sample + ".bamstats")
    cmd = "samtools flagstat --threads {} {} > {}".format(str(threads), input_bam, output_file)
    subprocess.call(cmd, Shell=True)

def create_coverage(input_bam, output_dir, sample):
    output_file = os.path.join(output_dir, sample + ".cov")
    cmd = "samtools depth -aa {} > {}".format(input_bam, output_file)
    subprocess.call(cmd, Shell=True)

# check if directory exists
def check_create_dir(path):

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """
    if os.path.exists(file_name):
        os.remove(file_name)

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
def check_argmunets(args):

    # Check mutation directory
    if not os.path.isdir(args.input_dir):
        print("%s: No such file or directory" %args.input_dir)
        return 1
    
    # Check reference genome
    if not os.path.isfile(args.reference):
        print("%s: No such file or directory" %args.reference)
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
