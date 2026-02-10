# IDs_NCLDVs_Bigna_add.bed
# scaffold_65     152116  152397  scaffold_65_30_3
# scaffold_39     109412  110713  scaffold_39_24_1
# scaffold_11     483024  483449  scaffold_11_121_2
# scaffold_6      1164847 1165521 scaffold_6_Genewise
# scaffold_87     36433   37107   scaffold_87_Genewise
# scaffold_10     402087  403472  scaffold_10_139_3

seqkit subseq --bed IDs_NCLDVs_Bigna_add.bed \
    Bigna1_masked_nuclear_scaffolds.fasta > IDs_NCLDVs_Bigna_add.fasta

for f in *.faa; 
do 
    base=${f%%.faa}
    mmseqs easy-search \
            --threads 40 -s 7 -e 1e-3 \
            IDs_NCLDVs_Bigna_add.fasta \
            $f \
            "$base"_add.o6 \
            /data/tmp \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
    mmseqs easy-search \
            --threads 40 -s 7 -e 1e-3 \
            Bigna1_filtered_proteins.fasta \
            $f \
            "$base"_prots.o6 \
            /data/tmp \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
done