#!/bin/bash
######################################################################################
indir=centrimo/summary ### inputdir
outdir=summary         ### outpudir
transcript=transcript          ### trascript data directory
ref=$1
######################################################################################
mkdir -p $outdir

type=$2.$ref
centrimo=$indir/$type.motif.pval.tab
centrimo2=$indir/$type.motif.pval_peak_cons.tab
out=${outdir}/$type.ave.fpkm.motif
out2=${out}.pval_peak_cons
tmp=${outdir}/$type.tmp

cat $transcript/ave.fpkm.tab |datamash transpose|sort -k1,1 > $tmp.1
cat $centrimo|datamash transpose|cut -f1|grep -v "^$"|sort|join -t$'\t' -a1 /dev/stdin $tmp.1 > $tmp.22
col=`head -1 $tmp.22|awk -F"\t" '{print NF}'`
cat $centrimo|cut -f1,2 | sort -k2,2 > $tmp.2
cat $tmp.22|awk 'OFS="\t"{tmp=$1;for(i=2;i<='$col';i++){if($i==""){tmp=tmp"\tNA"}else{tmp=tmp"\t"$i}};print tmp}'|datamash transpose|
    join -t$'\t' -1 2 $tmp.2 /dev/stdin|cut -f2,5-|sort -k1,1 > $out.tab
cat $tmp.1 | grep  -e "00geneid" -e $ref$'\t' |datamash transpose |join -t$'\t' -1 2 $tmp.2 /dev/stdin|cut -f2- > $out.ref
col2=`paste $centrimo2 $out.tab |head -1|awk '{print NF}'`
cat $centrimo2 |sort -k1,1 |paste /dev/stdin $out.tab|paste /dev/stdin $out.ref|cut -f1-${col2},$((col2 + 2)) /dev/stdin > $out2.tab ## output: centrimo(pval_peak_cons), transcript, transcript (reference)
rm $tmp.[123] $tmp.22

exit 0
