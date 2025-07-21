BASE=ChlorV-1
seqkit seq --min-qual ${params.min_read_quality} \
            --min-len ${params.min_read_length} \
            --out-file ${draft_name}.l${params.min_read_length}_q${params.min_read_quality}.${params.input_format} \
            fastq/*

flye -o "$BASE"-flye -t ${task.cpus} --meta --iterations 3 --nano-hq ${fastq} &> "$BASE"-flye.log

flye-post.py ${flye_dir} -m ${params.min_contig_length} -p "$BASE"
cp ${flye_dir}/"$BASE".fasta "$BASE".flye.fasta

run_medaka_parallel.py -a ${assembly} -f ${fastq} -p "$BASE"-medaka/"$BASE" -t ${task.cpus} -m ${params.medaka_model}
cp consensus.fasta "$BASE".medaka.fasta

# check assembly and identify giant virus genome
bfg "$BASE"_001 "$BASE".medaka.fasta > $BASE.fna

#if Illumina data is available:
run_polypolish.sh -g "$BASE".fna -t 40 "$BASE".R1.fastq.gz "$BASE".R2.fastq.gz