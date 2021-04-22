#!/bin/bash
export LC_ALL=C
while getopts w:i:q:d:l:r:c:g:s:m: OPT ;do
  case ${OPT} in
    w) wkdir=${OPTARG} ;;
    i) interval=${OPTARG} ;;
    q) q-value=${OPTARG} ;;
    d) diff=${OPTARG} ;;
    l) locount=${OPTARG} ;;
    r) refgene=${OPTARG} ;;
    c) refcpg=${OPTARG} ;;
    g) genome=${OPTARG} ;;
    s) chrom_sizes=${OPTARG} ;;
    m) motif=${OPTARG} ;;
  esac
done
shift $((OPTIND - 1))

refcell=$1
inputdir=$2

[[ -z ${wkdir} ]] && wkdir=$PWD
[[ -z ${interval} ]] && interval=200
[[ -z ${q_value} ]] && q_value=0.0001
[[ -z ${diff} ]] && diff=50
[[ -z ${locount} ]] && locount=100
echo $1 $2 $wkdir

srcdir=`echo $PWD/$0|sed 's|/[^/]*$||'`

cd ${wkdir}

for dir in `ls -d ${inputdir}/*` ;do
    bash ${srcdir}/average.merge.sh ${wkdir} ${dir##*/} ${inputdir}
done

for file in ${wkdir}/average/sum/* ;do
    if [[ file != ${wkdir}/average/sum/${refcell}.txt.gz ]];then
        cell2=`echo ${file##*/}|sed 's/\.txt\.gz//'`
        bash ${srcdir}/methylkit.sh ${wkdir} ${refcell} ${cell2} ${interval} ${q_value} ${diff} ${locount} ${srcdir}/methylkit.R ${refgene} ${refcpg}
        bash ${srcdir}/centrimo.sh ${wkdir} ${refcell} ${cell2} ${genome} ${chrom_sizes} ${motif}
    fi
done

bash ${srcdir}/centrimo.summary.sh ${wkdir} ${motifinput} ${motifseq} -230.259 1.2

#
bash ${srcdir}/transcript.summary.sh ${wkdir} 1


$outdir/$category.$mtype.motif.pval_peak_cons.tab
for file in ${wkdir}/centrimo/summary/* ;do
    if [[ file != ${wkdir}/centrimo/summary/${refcell}.txt.gz ]];then
        cell2=`echo ${file##*/}|sed 's/\.txt\.gz//'`
        bash ${srcdir}/centrimo.transcript ${wkdir} demeth es_cell..nih_roadmap
        bash ${srcdir}/centrimo.transcript ${wkdir} meth es_cell..nih_roadmap
    fi
done
