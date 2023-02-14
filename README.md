# SARS_COVinfections

A co-infection detection pipeline that focuses on SNP call frequency. 

In each specimen, the pipeline identifies genome positions where alternate allele frequency is ≥85% (homozygous SNPs); positions where more than one allele co-exists and major allele frequency is between 15 and 85% (heterozygous SNPs); and positions where the alternative allele is in a proportion of ≤15% and is considered to be background sequencing noise and therefore ruled out.

To run the pipeline. You first need to install the conda environment thourgh `SARS_COVinfections.yml`. Once installed, to run the code:

```
conda activate covidma

python Cov-infection.py -i "FASTQ_DIR" \
  -r "REF_GENOME" -o "OUTPUT_DIR" \
  -t "THREADS" -p "PRIMERS_Artic"
```
