#!/bin/bash
################################################################################################
wkdir=$1
category=$2
motifinput=$3
motifseq=$4
th=$5
peakth=$6
indir=$wkdir/centrimo/                       ## input dir
outdir=$wkdir/centrimo/summary         ## output dir
################################################################################################
cd $wkdir
mkdir -p $outdir

if [ ! -f $wkdir/image/motif.list ];then
  cat $motifinput | cut -f1 | sort |uniq > $wkdir/image/gene.list
  cat $motifinput | cut -f2 | sort |uniq > $wkdir/image/motif.list
  cat $motifinput | awk 'OFS="\t"{print $2,$1}'|sort -k1,1 > $wkdir/image/motif.gene.list
fi

fn2(){
  fn(){
    mtype=$1
    out=$outdir/$category.$mtype.motif.${th}_$peakth.tab
    out10=$outdir/$category.$mtype.motif.pval.tab
    out11=$outdir/$category.$mtype.motif.peak.tab
    out12=$outdir/$category.$mtype.motif.cons.tab
    header=$outdir/$category.$mtype.header.tab
    cat $motifinput | awk 'OFS="\t"{print $2,$1,$3}'|sort -k1,1 > $out
    cp -fp $out $out10
    cp -fp $out $out11
    cp -fp $out $out12
    echo -ne "00motif\t00geneid\t01type\t" > $header
    for f in `ls $indir/$category....*/centrimo_$mtype.tsv`
    do
      if [ -s $f ];then
        type=`dirname $f|sed 's|.*/||'|sed 's/'$category'\.\.\.\.//'`
        echo -ne $type"\t" >> $header
        ## calc. ratio of -80~-100pb/-120~-100pb
        cat ${f%/*}/site_counts_$mtype.txt|awk 'OFS="\t"{
          if($1~"DB"){printf $4"\t"}
          else if($1>=-120 && $1<-100){t+=$2;}
          else if($1>=-100 && $1<-80){t2+=$2;t3+=$2;}
          else if($1>=-80  && $1<100){t3+=$2;}
          else{if($1==100){print "-100..-80/-120..-100",t2/20,t3}t=0;t2=0;t3=0;}
        }'|sort -k1,1 > ${f%/*}/site_counts_$mtype.peak.tab

        cat $f|tail +2|grep -v ^#|grep -v ^$|cut -f2-|sort -k1,1|join -t$'\t' /dev/stdin ${f%/*}/site_counts_$mtype.peak.tab|sort -k7,7n > ${f%/*}/motif.site_counts_$mtype.peak.tab
        cat $f|tail +2|grep -v ^#|grep -v ^$|awk 'OFS="\t"{print $2,$7,$12}'|sort -k1,1 > $out10.tmp2
        join -t$'\t' $out10.tmp2 ${f%/*}/site_counts_$mtype.peak.tab|awk 'BEGIN{OFS="\t"}{print $1,$2}'|sort -k1,1 > $out10.tmp1
        join -t$'\t' $out10.tmp2 ${f%/*}/site_counts_$mtype.peak.tab|awk 'BEGIN{OFS="\t"}{print $1,$5}'|sort -k1,1 > $out11.tmp1
        join -t$'\t' $out10.tmp2 ${f%/*}/site_counts_$mtype.peak.tab|awk 'BEGIN{OFS="\t"}{print $1,$6/$3}'|sort -k1,1 > $out12.tmp1
        join -t$'\t' -a 1 $wkdir/image/motif.list $out10.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out10.tmp2
        join -t$'\t' -a 1 $wkdir/image/motif.list $out11.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out11.tmp2
        join -t$'\t' -a 1 $wkdir/image/motif.list $out12.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out12.tmp2
        join -t$'\t' $out10 $out10.tmp2 > $out10.tmp3
        join -t$'\t' $out11 $out11.tmp2 > $out11.tmp3
        join -t$'\t' $out12 $out12.tmp2 > $out12.tmp3
        mv -f $out10.tmp3 $out10
        mv -f $out11.tmp3 $out11
        mv -f $out12.tmp3 $out12
      fi
    done
    cat $header|sed 's/\t$/\n/'|cat /dev/stdin $out10 > $out10.tmp1; mv -f $out10.tmp1 $out10
    cat $header|sed 's/\t$/\n/'|cat /dev/stdin $out11 > $out11.tmp1; mv -f $out11.tmp1 $out11
    cat $header|sed 's/\t$/\n/'|cat /dev/stdin $out12 > $out12.tmp1; mv -f $out12.tmp1 $out12
    paste $out10 $out11 $out12|join -t$'\t' $motifseq /dev/stdin  > $outdir/$category.$mtype.motif.pval_peak_cons.tab
    rm $out10.tmp[1-4] $out11.tmp[1-4] $out12.tmp[1-4]

}
fn meth
fn demeth

