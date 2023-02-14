# SARS_COVinfections

A *co-infection detection pipeline* focusing on SNP call frequencies. 

In each sample, the pipeline identifies genomic positions where the alternative allele frequency is ≥85% (homozygous SNPs); positions where more than one allele coexists and the major allele frequency is between 15 and 85% (heterozygous SNPs); and positions where the alternative allele is ≤15% and is considered background sequencing noise and therefore excluded.

To run the pipeline. First you need to install the conda environment via `SARS_COVinfections.yml`. Once installed, run the code:

## 1. Activate conda environment

```{bash, eval=FALSE}
conda activate covidma
```

## 2. Run the python pipeline

```

python Cov-infection.py -i "FASTQ_DIR" \
  -r "REF_GENOME" -o "OUTPUT_DIR" \
  -t "THREADS" -p "PRIMERS_Artic"
```

### Mandatory arguments:

```
-i: input directory of paired fastq files.
-r: reference genome in fasta format.
-o: output directory. If it does not exist, it will be created.
-t: Number of threads.
-p: Primers used for sequencing.
```


### Optional arguments:

```
--min_DP: minimum frequency (depth) to accept a SNP.
--min_HOM: minimum proportion for homocygosis.
--ambiguity: minimum confident fraction to segregate.
--pangolin: pangolin annotation.
--snipit: snipit visualisation of SNPs.
```
To see all available options:

```
python Cov-infection.py -h
```

## Output generated

The output directories created are

* **Bam**: Bam file with reads mapped to reference.

* **Consensus**: 
  * ivar consensus fasta
  * Co-infection output** FASTQ_name folder: 
    * ALN: Visual representation of HTZ positions.
    * Sequences: Majority (2) and Minority (1) sequences.
    * Stats: Co-infection statistics.

* **Quality**: Fastq quality.

* **Stats**: Coverage and bam stats.

* **Trimmed**: Fastq trimmed.

* **Variants**: tsv file with SNVs.
