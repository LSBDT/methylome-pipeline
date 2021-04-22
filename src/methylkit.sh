#!/bin/sh
#run methylKit
#################################################
wkdir=$1
cell1=$2   #reference celltype
cell2=$3   #比較対象のサンプル
interval=$4
q=$5
diff=$6
locount=$7
rscript=$8
refgene=$9
refcpg=${10}
#################################################
outdir=${wkdir}/methylkit
inputdir=${wkdir}/average/sum
cd ${wkdir}
mkdir -p ${outdir} ${outdir}/QC $outdir/Rdata

sample1=${inputdir}/$cell1.txt.gz
sample2=${inputdir}/$cell2.txt.gz
out=${cell1}....${cell2}

zcat ${sample1} > ${sample1}.${cell2}.txt
zcat ${sample2} > ${sample2}.${cell1}.txt
R --slave --args $wkdir $interval $outdir $diff $q ${sample1}.${cell2}.txt ${sample2}.${cell1}.txt $out ${cell1} ${cell2} $locount $refgene $refcpg < $rscript 1>&2
gzip -f ${outdir}/${out}.txt
rm ${sample1}.${cell2}.txt ${sample2}.${cell1}.txt
