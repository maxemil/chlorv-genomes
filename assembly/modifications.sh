
dorado -vv 2> $out_dir/$run.dorado.log
dorado basecaller -v "$model",$mod $pod5_dir --models-directory /data/models/ $dorado_options 2>> $out_dir/$run.dorado.log | samtools view -F 3584 -u -S -@ 20 - | samtools sort -@ 20 -o $out_dir/"$run"_mod.bam
samtools index $out_dir/"$run"_mod.bam

BASE='ChlorV-1'
MAPPED_MOD='ChlorV-1_mapped.bam'

~/software/nanopore-basecalling/align_modifications.sh -m ChlorV-1_mod.bam \
            -g $BASE.fna \
            -a ChlorV-1.medaka.fasta \
            -c ChlorV-1..001 \
            -o $MAPPED_MOD

samtools view -b -@ 30 $MAPPED_MOD $BASE > $BASE.bam
samtools index $BASE.bam
modkit pileup $BASE.bam $BASE.bed --ref $BASE.fna -t 40 --log-filepath $BASE.modkit.log

# this part was inspired by https://github.com/AlexdeMendoza/6mA_evolution
samtools faidx $BASE.fna
bedtools slop -l 1 -r 8 -s -i $BASE.bed -g $BASE.fna.fai | bedtools merge -i - > $BASE.plus3.bed
bedtools getfasta -fi $BASE.fna -fo $BASE.plus3.fasta -s -bed $BASE.plus3.bed
seqkit fx2tab $BASE.plus3.fasta | cut -f 2 > $BASE.plus3.fx2tab
paste $BASE.bed $BASE.plus3.fx2tab > $BASE.6mA.context.bed

modkit motif search -i $BASE.bed -r $BASE.fna -o ./"$BASE".motif_search.tsv --threads 30
# modkit motif bed $BASE.fna CANNNNNNTG 1 > "$BASE"_motif.bed
modkit motif bed $BASE.fna TGCA 1 > "$BASE"_motif.bed