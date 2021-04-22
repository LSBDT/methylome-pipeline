
## requirements

### tools
bedtools
datamash
meme-5.0.1 (centrimo,fasta-get-markov)
R
R library methylKit
R library graphics
R library genomation

### files
genome reference hg19.fa
chrom.sizes file
refGene bed file
CpG island bed file
motif file (meme format)
gene name motif list file
motif sequence file

## usege

./pipeline.sh {methylkit input directory} {transcript input directory}

input directory contains methylkit input data

### input files

{methylkit input directory}/{celltype}/{samplename}.gz

{transcript input directory}/{celltype}/{samplename}.tab
column1: genename
column2: fpkm

## test
./pipeline.sh testdata/input/methyl testdata/input/transcript
