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


### reference files

- genome reference hg19.fa, hg19.fa.fai
- chrom.sizes file
- refGene bed file
- CpG island bed file
- motif file (meme format)
- gene name motif list file
- motif sequence file

## usege

```./pipeline.sh -d {reference directory} {reference cell name} {methylkit input directory} {transcript input directory}```

input directory contains methylkit input data

### input files

```{methylkit input directory}/{celltype}/{samplename}.gz```

The files are gzipped tab-delimited texts with the following columns.

- column1: chr
- column2: start
- column3: end
- column4: total coverage
- column5: methyl coverage
- column6: methyl %

```{transcript input directory}/{celltype}/{samplename}.tab```

The files are tab-delimited texts with the following columns.

- column1: genename
- column2: fpkm

### output files

```{workdir}/centrimo/summary/{reference_cell}.{meth,demeth}.ave.fpkm.motif.pval_peak_cons.tab```

- meth: methylated in reference_cell
- demeth: demethylated in reference_cell

- rows: motifs
- columns: p-value, C-value, concentration, fpkm


## for test run

~~~~
cd testdata/work
../../src/pipeline.sh -d ../db celltype1 ../input/methyl ../input/transcript/
~~~~


## Datasets used in the report

The input files to create the original dataset used in the report by Miyajima and Noguchi et al. are available at the following site:

https://genomec.gsc.riken.jp/gerg/owncloud/index.php/s/o8N21uW976my2t8

