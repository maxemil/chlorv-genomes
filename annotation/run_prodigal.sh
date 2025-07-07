#!/bin/bash

## subs
usage(){
cat <<EOF
Usage:
  run_prodigal.sh FASTA 
run prodigal on genome.
  -p  mode, either meta or single (default single).
  -g  run prodigal-gv instead of vanilla prodigal.
  -h  show this help
EOF
exit 0;
}

## prep
[[ $# -eq 0 ]] && usage;

# Execute getopt
ARGS=`getopt --name "run_prodigal.sh" \
    --options "p:gh" \
    -- "$@"`

#Bad arguments
[ $? -ne 0 ] && exit 1;

# A little magic
eval set -- "$ARGS"

MODE='single'
CMD='prodigal'
# Now go through all the options
while true; do
    case "$1" in
        -p)
            [ ! -n "$2" ] && (echo "$1: value required" 1>&2 && exit 1);
            MODE=$2;
            shift 2;;
        -g)
            CMD="prodigal-gv";
            shift;;
        -h)
	    usage && exit 0;;
        --)
            shift
            break;;
        *)
            echo "$1: Unknown option" 1>&2 && exit 1;;
    esac
done

GENOME=$1
EXT="${GENOME##*.}"
BASE=${GENOME%%.$EXT}

echo "$CMD -i $GENOME -d $BASE.ffn -a $BASE.faa -o $BASE.gff -f gff -p $MODE"
$CMD -i $GENOME -d $BASE.ffn -a $BASE.faa -o $BASE.gff -f gff -p $MODE 