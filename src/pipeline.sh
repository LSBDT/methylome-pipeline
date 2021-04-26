#!/bin/bash

#test command
#src/pipeline.sh -w /mnt/d/riken/methyl/git/testdata/work celltype1 /mnt/d/riken/methyl/git/testdata/input/methyl
#pipeline.sh -w /mnt/d/riken/methyl/git/testdata/work celltype1 /mnt/d/riken/methyl/git/testdata/input/methyl /mnt/d/riken/methyl/git/testdata/input/transcript/
#  export PATH=$PATH:/home/ngc/miniconda3/envs/bedtools/bin/
# ../../pipeline.sh -d ../db celltype1 ../input/methyl ../input/transcript/


srcdir=$(cd $(dirname $(readlink -f $0 || echo $0));pwd -P)
usage(){
    cat <<EOF
Usage:
    $(basename ${0}) [<options>] [-r refGene bed file] [-c CpG bed file] [methylation data directory] [transcript data directory]

Options:

    -w    working directory [def. current directory]
    -b    reference data directory
    -i    interval(bp) [def. 200]
    -q    q-value threshold [def. 0.0001]
    -d    differential methylation % threshold [def. 50]
    -l    methylkit count threshold [def. 100]
    -r    refGene bed file
    -p    CpG bed file
    -g    genome fasta file
    -c    genome chrom sizes file
    -m    motif list file
    -n    genename motif list file
    -s    motif sequence file
    -h    print this
EOF
    exit 1
}

[[ $# = 0 ]] && usage

export LC_ALL=C
while getopts w:d:i:q:d:l:r:p:g:c:m:n:s:h OPT ;do
  case ${OPT} in
    w) wkdir=${OPTARG} ;;
    d) datadir=${OPTARG} ;;
    i) interval=${OPTARG} ;;
    q) q-value=${OPTARG} ;;
    d) diff=${OPTARG} ;;
    l) locount=${OPTARG} ;;
    r) refgene=${OPTARG} ;;
    p) refcpg=${OPTARG} ;;
    g) genome=${OPTARG} ;;
    c) chrom_sizes=${OPTARG} ;;
    m) motif_list=${OPTARG} ;;
    n) gene_motif=${OPTARG} ;;
    s) motifseq=${OPTARG} ;;
    h) usage ;;
  esac
done
shift $((OPTIND - 1))
[[ -z ${wkdir} ]] && wkdir=$PWD
cd ${wkdir}

bedtools &> /dev/null || { echo "bedtools not found."; exit 1; }

check_file(){
    [[ -f $1 ]] && return
    echo "$2 file $1 dose not exist." 1>&2; exit 1;
}
check_dir(){
    [[ -d $1 ]] && return
    echo "$2 directory $1 dose not exist." 1>&2; exit 1;
}
refcell=$1
inputdir=$2
transcriptdir=$3
: ${refcell:?} || echo "no reference cell" 1>&2
check_dir ${inputdir} "methylation data directory"
check_dir ${transcriptdir} "transcript data directory"

[[ -z ${interval} ]] && interval=200
[[ -z ${q_value} ]] && q_value=0.0001
[[ -z ${diff} ]] && diff=50
[[ -z ${locount} ]] && locount=100

if [[ -n ${datadir} ]];then
    refgene=${datadir}/hg19.refGene.bed
    refcpg=${datadir}/cpgi.hg19.bed
    genome=${datadir}/hg19.fa.gz
    chrom_sizes=${datadir}/hg19.chrom.sizes
    motif_list=${datadir}/image.meme
    gene_motif=${datadir}/Genename_Motif.txt
    motifseq=${datadir}/motif.seq
else
    check_file ${refgene} "[-r refGene bed file]"
    check_file ${refcpg}  "[-p CpG bed file]"
    check_file ${genome}  "[-g genome fasta file]"
    check_file ${chrom_sizes} "[-c genome chrom sizes file]"
    check_file ${motif_list}   "[-m motif list file]"
    check_file ${gene_motif}  "[-n genename motif list file]"
    check_file ${motifseq}    "[-s motif sequence file]"
fi
#test(){
echo "methylation data prep"
for dir in `ls -d ${inputdir}/*` ;do
    [[ -d ${dir} ]] || { echo "${inputdir}/* dose not exist." 1>&2; exit 1; }
    bash ${srcdir}/average.merge.sh ${dir##*/} ${inputdir}
done
[[ $? = 0 ]] || exit 1

for file in average/*.txt.gz ;do
    if [[ ! ${file} = average/${refcell}.txt.gz ]];then
        cell2=`echo ${file##*/}|sed 's/\.txt\.gz//'`
        echo "run methylKit ${refcell} ${cell2}"
        bash ${srcdir}/methylkit.sh ${refcell} ${cell2} ${interval} ${q_value} ${diff} ${locount} ${srcdir}/methylkit.R ${refgene} ${refcpg}
        echo "run centrimo ${refcell} ${cell2}"
        bash ${srcdir}/centrimo.sh ${refcell} ${cell2} ${genome} ${chrom_sizes} ${motif_list}
    fi
done
#}
echo "centrimo summary"
bash ${srcdir}/centrimo.summary.sh ${refcell} ${gene_motif} ${motifseq} -230.259 1.2
[[ $? = 0 ]] || exit 1

echo "transcript summary"
bash ${srcdir}/transcript.summary.sh ${transcriptdir}
[[ $? = 0 ]] || exit 1

echo "centrimo transcript summary"
bash ${srcdir}/centrimo.transcript.sh demeth ${refcell}
bash ${srcdir}/centrimo.transcript.sh meth ${refcell}
