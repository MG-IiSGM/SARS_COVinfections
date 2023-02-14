# SARS_COVinfections

A co-infection detection pipeline that focuses on SNP call frequency. 

In each specimen, the pipeline identifies genome positions where alternate allele frequency is ≥85% (homozygous SNPs); positions where more than one allele co-exists and major allele frequency is between 15 and 85% (heterozygous SNPs); and positions where the alternative allele is in a proportion of ≤15% and is considered to be background sequencing noise and therefore ruled out.

To run the pipeline. First you need to install the conda environment via `SARS_COVinfections.yml`. Once installed, run the code:

```{bash, eval=FALSE}
conda activate covidma

python Cov-infection.py -i "FASTQ_DIR" \
  -r "REF_GENOME" -o "OUTPUT_DIR" \
  -t "THREADS" -p "PRIMERS_Artic"
```

Additional options:

```
--min_DP: minimum frequency (depth) to accept a SNP.
--min_HOM: minimum proportion for homocygosis.
--ambiguity: min confident proportion to segregate.
--pangolin: pangolin annotation
--snipit: snipit visualization of SNPs
```
To see all available options:

```
python Cov-infection.py -h
```

The output directories created are:

* Bam: Bam file with reads mapped to reference.
* Consensus: 
            * ivar consensus fasta
            * FASTQ_name folder: 
                                  * ALN: Visual representation of HTZ positions.
                                  * Sequences: Mayority (2) and Minority (1) sequences.
                                  * Stats: Co-infection stats.
* Quality: Fastq quality.
* Stats: Coverage and Bam stats.
* Trimmed: Fastq trimmed.
* Variants: tsv file with SNVs
