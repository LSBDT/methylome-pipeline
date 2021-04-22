#!/bin/bash
#######################################################################
wkdir=$1
th=$2
inputdir=$3
genelist=${wkdir}/image/gene.list
outdir=${wkdir}/transcript/summary
#######################################################################
mkdir -p ${wkdir} $outdir
cd ${wkdir}

fn(){
    count=2
    for f in ${wkdir}/$dd/*.tab ;do
        if [[ $count == 2 ]];then
            join -t$'\t' -a 1 $genelist $f|awk 'OFS="\t"{if($'$count'==""){print $1,"0"}else{print $0}}' > tmp2
        else
            join -t$'\t' -a 1 tmp1 $f|awk 'OFS="\t"{if($'$count'==""){print $0,"0"}else{print $0}}' > tmp2
        fi
        mv -f tmp2 tmp1
        ((count+=1))
    done
    cat tmp1 | awk 'BEGIN{OFS="\t"}{
    ave=0
    for(i=2;i<=NF;i++){
      ave += $i
    }
    print $1,ave/(NF-1)}' > $outdir/$dd.ave.tab
    rm tmp1
}

for d in ${inptdir}/* ;do
  dd=`echo ${d##*/}`
  fn
done

fn2(){
    if [[ $count == 1 ]];then
        echo -e "00geneid\t"$ff > header1
        cat $f > tmp3
    else
        cat header1|awk 'OFS="\t"{print $0,"'"$ff"'"}' > header2
        mv -f header2 header1
        join -t$'\t' tmp3 $f > tmp4
        mv -f tmp4 tmp3
    fi
}

count=1
for f in $outdir/*.ave.tab ;do
    ff=`echo ${f##*/}|sed 's/\.ave\.tab$//'`
    fn2
    ((count+=1)
done
cat header1 tmp3 > $outdir/ave.fpkm.tab
rm header1 tmp[3]

