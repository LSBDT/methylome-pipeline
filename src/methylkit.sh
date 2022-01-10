#!/bin/bash
#################################################
cell1=$1   #reference celltype
cell2=$2   #sample
interval=$3
q=$4
diff=$5
locount=$6
rscript=$7
genomever=$8
refgene=$9
refcpg=$10
#################################################
outdir=methylkit
inputdir=average
mkdir -p ${outdir} ${outdir}/QC

sample1=${inputdir}/$cell1.txt.gz
sample2=${inputdir}/$cell2.txt.gz
out=${cell1}....${cell2}

zcat ${sample1} > ${sample1}.${cell2}.txt
zcat ${sample2} > ${sample2}.${cell1}.txt
echo " --vanilla --args $PWD $interval $outdir $diff $q ${sample1}.${cell2}.txt ${sample2}.${cell1}.txt $out ${cell1} ${cell2} $locount $genomever $refgene $refcpg < $rscript 1>&2"
R --vanilla --args $PWD $interval $outdir $diff $q ${sample1}.${cell2}.txt ${sample2}.${cell1}.txt $out ${cell1} ${cell2} $locount $genomever $refgene $refcpg < $rscript 1>&2
gzip -f ${outdir}/${out}.txt
rm ${sample1}.${cell2}.txt ${sample2}.${cell1}.txt

exit 0
