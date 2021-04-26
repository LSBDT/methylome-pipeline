#!/bin/bash
#################################
refcell=$1
motifinput=$2
motifseq=$3
th=$4
peakth=$5
indir=centrimo/
outdir=centrimo/summary
#################################
mkdir -p $outdir

if [[ ! -f image/motif.list ]];then
  mkdir -p image
  cat $motifinput | cut -f1 | sort |uniq > image/gene.list
  cat $motifinput | cut -f2 | sort |uniq > image/motif.list
  cat $motifinput | awk 'OFS="\t"{print $2,$1}'|sort -k1,1 > image/motif.gene.list
fi

fn(){
    mtype=$1
    out=$outdir/${refcell}.$mtype.motif.${th}_$peakth.tab
    out10=$outdir/${refcell}.$mtype.motif.pval.tab
    out11=$outdir/${refcell}.$mtype.motif.peak.tab
    out12=$outdir/${refcell}.$mtype.motif.cons.tab
    header=$outdir/${refcell}.$mtype.header.tab
    cat $motifinput | awk -F"\t" 'OFS="\t"{print $2,$1,$3}'|sort -k1,1 > $out
    cp -fp $out $out10
    cp -fp $out $out11
    cp -fp $out $out12
    echo -ne "00motif\t00geneid\t01type\t" > $header
    for f in `ls $indir/${refcell}....*/centrimo_$mtype.txt`;do
      if [[ -s $f ]];then
        type=`dirname $f|sed 's|.*/||'|sed 's/'${refcell}'\.\.\.\.//'`
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
        join -t$'\t' -a 1 image/motif.list $out10.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out10.tmp2
        join -t$'\t' -a 1 image/motif.list $out11.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out11.tmp2
        join -t$'\t' -a 1 image/motif.list $out12.tmp1|awk 'OFS="\t"{if($2==""){print $1,0}else{print $1,$2}}' > $out12.tmp2
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
    paste $out10 $out11 $out12|join -t$'\t' $motifseq /dev/stdin  > $outdir/${refcell}.$mtype.motif.pval_peak_cons.tab
    rm $out10.tmp[1-2] $out11.tmp[1-2] $out12.tmp[1-2]
}
fn meth
fn demeth

exit 0
