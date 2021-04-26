#!/bin/bash
#########################
cell=$1
indir=$2
outdir=average
#########################
mkdir -p ${outdir}

input_file=${indir}/${cell}/*.cov.meth.gz
count=`ls ${input_file}|wc -l`
zcat ${input_file}|cut -f1-5 |sort -k1,1 -k2,2n > ${outdir}/${cell}.sum.bed
bedtools merge -d -1 -i ${outdir}/${cell}.sum.bed -c 4,5 -o sum,sum |
    awk 'BEGIN{OFS="\t";OFMT="%.2f"}{meth=100*$5/$4;print $1"."$2+1,$1,$2+1,"F",int($4/'$count'),meth,100-meth}' |
    gzip -c > ${outdir}/${cell}.txt.gz
rm  ${outdir}/${cell}.sum.bed
exit 0
