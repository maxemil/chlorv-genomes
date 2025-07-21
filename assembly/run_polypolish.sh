#!/bin/bash
set -e

## subs
usage(){
cat <<EOF
Usage:
  run_polypolish.sh -g FASTA FASTQ1 [FASTQ2 ... ]
run polypolish on genome assebly.
  -g  genome fasta file
  -t number of threads for BWA
  FASTQ fastq files with Illumina data
  -h  show this help
EOF
exit 0;
}

## prep
[[ $# -eq 0 ]] && usage;

# Execute getopt
ARGS=`getopt --name "run_polypolish.sh" \
    --options "g:t:h" \
    -- "$@"`
#Bad arguments
[ $? -ne 0 ] && exit 1;

# A little magic
eval set -- "$ARGS"

CPU=10

# Now go through all the options
while [ : ]; do
    case "$1" in
        -g)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            GENOME="$2";
            shift 2;;
        -t)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            CPU="$2";
            shift 2;;
        -h)
	    usage && exit 0;;
        --)
            shift
            break;;
        *)
            echo "$1: Unknown option" 1>&2 && exit 1;;
    esac
done

FASTQ=( $@ )
EXT="${GENOME##*.}"
BASE=${GENOME%%.$EXT}


bwa index $GENOME 2>> $BASE.polypolish_bwa.log

for((i=0; i < ${#FASTQ[@]}; i+=2))
do
    p=$((i/2+1))
    echo "processing pair $p" 
    fq1=( "${FASTQ[@]:i}" )
    fq2=( "${FASTQ[@]:i+1}" )
    bwa mem -t $CPU -a $GENOME $fq1 > "$BASE"_"$p"-1.sam 2>> $BASE.polypolish_bwa.log
    bwa mem -t $CPU -a $GENOME $fq2 > "$BASE"_"$p"-2.sam 2>> $BASE.polypolish_bwa.log
    polypolish filter --in1 "$BASE"_"$p"-1.sam --in2 "$BASE"_"$p"-2.sam --out1 "$BASE"_"$p"-1.filtered.sam --out2 "$BASE"_"$p"-2.filtered.sam 2> $BASE.polypolish.log
done
echo "polishing assembly" 

polypolish polish --debug $BASE.polypolish.tsv $GENOME "$BASE"_*.filtered.sam 2>> $BASE.polypolish.log > $BASE.polypolish.$EXT
