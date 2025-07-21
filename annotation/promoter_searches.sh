GENOME='BigV-3'

# get late motif
get_upstream_noncoding_region.py upstream -f ../"$GENOME".fna -g ../ORF_calling/"$GENOME".gff -o "$GENOME".upstream40.fasta -l 40
bfg --search-sequences TCTA $GENOME.upstream40.fasta > $GENOME.upstream40.TCTA.fasta
/home/meme/bin/meme -oc "$GENOME"_upstream40 -dna -minw 2 -maxw 9 $GENOME.upstream40.TCTA.fasta

# get early motif
upl='65'
get_upstream_noncoding_region.py upstream -f ../"$GENOME".fna -g ../ORF_calling/"$GENOME".gff -o "$GENOME".upstream"$upl".fasta -l $upl
bfg -v --search-sequences TCTA $GENOME.upstream"$upl".fasta > $GENOME.upstream"$upl".noTCTA.fasta
/home/meme/bin/meme -oc "$GENOME"_upstream"$upl" -dna -minw 2 -maxw 11 -allw $GENOME.upstream"$upl".noTCTA.fasta

# combine motifs
cat <(/home/meme/libexec/meme-5.5.0/meme-get-motif -id MEME-1 -ia "$GENOME"_upstream"$upl"/meme.txt)\
    <(/home/meme/libexec/meme-5.5.0/meme-get-motif -id MEME-1 -ia "$GENOME"_upstream40/meme.txt) |
    /home/meme/libexec/meme-5.5.0/meme-get-motif -all > early_late_motifs.txt

# find all individual motif occurrences
/home/meme/bin/centrimo --oc "$GENOME"_upstream50_centrimo "$GENOME".upstream50.fasta early_late_motifs.txt
/home/meme/bin/fimo --oc "$GENOME"_upstream50_fimo early_late_motifs.txt "$GENOME".upstream50.fasta