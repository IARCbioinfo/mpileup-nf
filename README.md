# mpileup-nf
Nextflow pipeline for coverage computation with samtools mpileup (bed parallelization) 

Command line example:

```
nextflow run iarcbioinfo/mpileup-nf --bam_folder BAM/ --fasta_ref ref.fasta  --bedfile bedfile.bed --nsplit 100 --out_table  output.txt
```

Output tab-delimited file looks like:
```
  chr pos K211007X K211008Y K211010Z
  chr1 801957 314 252 287
  chr1 866511 212 128 169
```
with samples in columns and positions in lines.
