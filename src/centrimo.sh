#!/bin/bash
#################################################
cell1=$1
cell2=$2
genome=$3
chrom=$4
motif=$5
#################################################
defoutdir=centrimo
inputdir=methylkit
mkdir -p $defoutdir

infile=${inputdir}/${cell1}....${cell2}.txt.gz
outdir=${defoutdir}/${cell1}....${cell2}
mkdir -p ${outdir}

zcat ${infile}|cut -f2-|grep ^chr | sort -k1,1 -k2,2n >  ${outdir}/tmp.txt
cat ${outdir}/tmp.txt | awk 'BEGIN{OFS="\t"}{if($7>0) print $1,$2,$3,$1":"$2"-"$3}' > ${outdir}/${cell1}_demeth.bed
cat ${outdir}/tmp.txt | awk 'BEGIN{OFS="\t"}{if($7<0) print $1,$2,$3,$1":"$2"-"$3}' > ${outdir}/${cell1}_meth.bed
rm ${outdir}/tmp.txt

fn1(){
#    if [[ $(cat ${outdir}/${cell1}_$1.bed | wc -l) -ge 100 ]];then
        cat ${outdir}/${cell1}_$1.bed |
            bedtools flank -i /dev/stdin -g ${chrom} -b 1000 |
            awk 'OFS="\t"{if(NR%2==1){a1=$1;a2=$2}else{print a1,a2,$3,$4}}'|
            bedtools getfasta -name -fi ${genome%.gz} \
            -bed /dev/stdin -fo ${outdir}/${cell1}_$1.fasta
        fasta-get-markov ${outdir}/${cell1}_$1.fasta -m 1 > ${outdir}/bg_$1.txt
        centrimo -oc ${outdir} --bfile ${outdir}/bg_$1.txt ${outdir}/${cell1}_$1.fasta ${motif}
        mv ${outdir}/centrimo.html   ${outdir}/centrimo_$1.html
        mv ${outdir}/centrimo.txt    ${outdir}/centrimo_$1.txt
        mv ${outdir}/site_counts.txt ${outdir}/site_counts_$1.txt
#        rm ${outdir}/${cell1}_$1.fasta
#    fi
}
fn1 meth
fn1 demeth

exit 0
