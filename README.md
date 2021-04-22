# The pipeline for building a dataset of differentially methylated regions and pioneering TFs

The pipeline build a dataset of differentially methylated regions in the human genome and predicting the pioneering transcription factors that cause the demethylation and differentiation. The pipeline uses IHEC datasets and IMAGE TFBS data.


## requirements

### tools

- bedtools
- datamash
- meme-5.0.1 (centrimo,fasta-get-markov)
- R
- R library methylKit
- R library graphics
- R library genomation


### files

- genome reference hg19.fa
- chrom.sizes file
- refGene bed file
- CpG island bed file
- motif file (meme format)
- gene name motif list file
- motif sequence file

## usege

```./pipeline.sh {methylkit input directory} {transcript input directory}```

input directory contains methylkit input data

### input files

```{methylkit input directory}/{celltype}/{samplename}.gz```

```{transcript input directory}/{celltype}/{samplename}.tab```

The '.tab' file is a tab-delimited text whose columns are
- column1: genename
- column2: fpkm

## for test run

```./pipeline.sh testdata/input/methyl testdata/input/transcript```


## Datasets used in the report

The input files to create the original dataset used in the report by Miyajima and Noguchi et al. are available at the following site:

https://genomec.gsc.riken.jp/gerg/owncloud/index.php/s/o8N21uW976my2t8


