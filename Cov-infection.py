# imports
import argparse, os, utils, misc, sys

# parse arguments (-B)
parser = argparse.ArgumentParser(description="Script to infer a potential co-infection")

# INPUT parameters
parser.add_argument("-i", "--input_dir", metavar="input_directory",
                    type=str, required=True, help="Input directory containing all fastq files")
parser.add_argument("-r", "--reference", metavar="reference",
                    type=str, required=True, help="Reference genome to map fastq files")
parser.add_argument("-o", "--out_dir", help="Output directory", type=str, 
                    required=False, default=".")
parser.add_argument('-t', '--threads', type=str, dest="threads",
                                  required=False, default=4, help='Threads to use')
parser.add_argument('-p', '--primers', type=str, default='/home/laura/DATABASES/Anotacion/COVID/primers/nCoV-2019.bed',
                                 required=False, help='Bed file including primers to trim')

# COINFECTION ARGUMENTS
parser.add_argument("--min_DP", help="minimum frequency (depth) to accept a SNP", default=10,
                    type=int)
parser.add_argument("--min_HOM", help="minimum proportion for homocygosis",
                    default=0.85, type=float)
parser.add_argument("--ambiguity", help="min value to segregate", default=0.55, type=float)
parser.add_argument("--max_mean_htz", help="maximum htz proportion",
                    default=0.75, type=float)
parser.add_argument("--max_std_htz", help="maximum deviation of htz proportion",
                    default=0.08, type=float)
parser.add_argument("--max_extra_std_htz", help="maximum extra deviation of htz proportion",
                    default=0.015, type=float)
parser.add_argument("--SNPs_out", help="percentage of SNPs out mean htz +- (std htz + extra_std)",
                    default=0.3, type=float)

parser.add_argument("--pangolin", help="pangolin annotation", 
                    action="store_true")
parser.add_argument("--snipit", help="snipit analysis", 
                    action="store_true")

# main
def main():

    # Script directory
    abs_path = os.path.abspath(sys.argv[0])
    script_dir = os.path.dirname(abs_path)

    # check arguments
    args = parser.parse_args()

    if utils.check_argmunets(args):
        exit(1)

    else:

        # Check out dir
        utils.check_create_dir(args.out_dir)

        # Output subdirectories
        out_qc_dir = os.path.join(args.out_dir, "Quality")
        out_qc_pre_dir = os.path.join(out_qc_dir, "raw")
        out_qc_post_dir = os.path.join(out_qc_dir, "processed")
        out_trim_dir = os.path.join(args.out_dir, "Trimmed")
        out_map_dir = os.path.join(args.out_dir, "Bam")
        out_variant_dir = os.path.join(args.out_dir, "Variants")
        out_consensus_dir = os.path.join(args.out_dir, "Consensus")
        out_stats_dir = os.path.join(args.out_dir, "Stats")
        out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")
        out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")


        # Obtain all R1 and R2 from folder
        r1, r2 = utils.extract_read_list(args.input_dir)
    
    # Loop for each fastq
    for r1_file, r2_file in zip(r1, r2):

        sample = utils.extract_sample(r1_file, r2_file)
        print(sample)

        # QUALITY CHECK in RAW with fastqc
        ######################################################
        utils.check_create_dir(out_qc_dir)
        utils.fastqc_quality(r1_file, r2_file,
                                out_qc_pre_dir, args.threads)
        
        # QUALITY TRIMMING AND ADAPTER REMOVAL WITH fastp
        ###################################################
        out_trim_name_r1 = sample + ".trimmed_R1.fastq.gz"
        out_trim_name_r2 = sample + ".trimmed_R2.fastq.gz"
        output_trimming_file_r1 = os.path.join(
            out_trim_dir, out_trim_name_r1)
        output_trimming_file_r2 = os.path.join(
            out_trim_dir, out_trim_name_r2)
        utils.fastp_trimming(r1_file, r2_file, sample, out_trim_dir,
                threads=args.threads, min_qual=20, window_size=10, min_len=35)
        
        # QUALITY CHECK in TRIMMED with fastqc
        ######################################################
        utils.fastqc_quality(output_trimming_file_r1, output_trimming_file_r2,
         out_qc_post_dir, args.threads)
        
        # MAPPING WITH BWA - SAM TO SORTED BAM - ADD HEADER SG
        #####################################################
        out_map_name = sample + ".rg.sorted.bam"
        output_map_file = os.path.join(out_map_dir, out_map_name)
        utils.bwa_mapping(output_trimming_file_r1, output_trimming_file_r2,
            args.reference, sample, out_map_dir, threads=args.threads)
        utils.sam_to_index_bam(
            sample, out_map_dir, output_trimming_file_r1, threads=args.threads)
        
         #MARK DUPLICATES WITH PICARDTOOLS ###################
        #####################################################
        out_markdup_name = sample + ".rg.markdup.sorted.bam"
        output_markdup_file = os.path.join(
            out_map_dir, out_markdup_name)
        utils.picard_markdup(output_map_file)

        #TRIM PRIMERS WITH ivar trim ########################
        #####################################################
        utils.ivar_trim(output_markdup_file, args.primers, sample,
                min_length=30, min_quality=20, sliding_window_width=4)
        
        #VARIANT CALLING WTIH ivar variants##################
        #####################################################
        utils.check_create_dir(out_variant_dir)
        out_ivar_variant_name = sample + ".tsv"
        out_ivar_variant_file = os.path.join(
            out_variant_dir, out_ivar_variant_name)
        out_markdup_trimmed_name = sample + ".rg.markdup.trimmed.sorted.bam"
        output_markdup_trimmed_file = os.path.join(
                out_map_dir, out_markdup_trimmed_name)
        utils.ivar_variants(args.reference, output_markdup_trimmed_file, out_variant_dir, sample,
                min_quality=0, min_frequency_threshold=0.1, min_depth=0)
        
        #CREATE CONSENSUS with ivar consensus##################
        #######################################################
        utils.check_create_dir(out_consensus_dir)
        out_ivar_consensus_name = sample + ".fa"
        out_ivar_consensus_file = os.path.join(
                out_consensus_dir, out_ivar_consensus_name)
        utils.ivar_consensus(output_markdup_trimmed_file, out_consensus_dir, sample,
                    min_quality=20, min_frequency_threshold=0.8, min_depth=10, uncovered_character='N')
        utils.replace_consensus_header(out_ivar_consensus_file)


        ########################CREATE STATS ###################
        ########################################################
        utils.check_create_dir(out_stats_dir)
        utils.check_create_dir(out_stats_bamstats_dir)
        utils.check_create_dir(out_stats_coverage_dir)
        out_coverage_name = sample + ".cov"
        out_coverage_file = os.path.join(
            out_stats_coverage_dir, out_coverage_name)
        utils.create_bamstat(output_markdup_trimmed_file,
                out_stats_bamstats_dir, sample, threads=args.threads)
        utils.create_coverage(output_markdup_trimmed_file,
                                out_stats_coverage_dir, sample)


        # parse variant calling file
        df = misc.parse_vcf(args, out_ivar_variant_file)

        # get alingment
        misc.get_alingment(args, out_consensus_dir, sample, df, out_coverage_file)

        # Get stats
        utils.quality_control(df, args, sample, out_consensus_dir)
    
    # If all OK
    exit(0)

if __name__ == "__main__":
    main()
