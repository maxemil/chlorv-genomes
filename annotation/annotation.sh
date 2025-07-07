# ORF Calling

GENOME='BigV-4'
mkdir ORF_calling
cd ORF_calling
ln -s $(realpath ../$GENOME.fna $GENOME.fna)
run_prodigal.sh -p single -g $GENOME.fna
mkdir prodigal
mv $GENOME.faa $GENOME.gff $GENOME.ffn prodigal
/software/genemark_suite_linux_64/gmsuite/gmsn.pl --format GFF3 --output $GENOME.gm.gff --name $GENOME.gm --faa --clean --virus $GENOME.fna
mkdir genemark
mv $GENOME.gm* gms.log genemark
merge_gff.py -p prodigal/$GENOME.gff -g genemark/$GENOME.gm.gff -f $GENOME.fna -s '..'
gffread -g $GENOME.fna -x $GENOME.ffn -y $GENOME.faa $GENOME.merge.gff -s '..' -l 100
mv $GENOME.merge.gff $GENOME.gff
cd ..


# run sequence-based annotation

GENOME='BigV-4'
mkdir annotation
ln -s $(realpath ORF_calling/$GENOME.faa annotation)
cd annotation
seq-gc -Nbw 50 ../$GENOME.fna > $GENOME.gc.tsv
emapper.py -m mmseqs \
        --report_no_hits --data_dir /databases/eggnog-mapper \
        --temp_dir /data/tmp/ --cpu 20 -i $GENOME.faa \
        -o "$GENOME"
hmmsearch --cpu 30 --tblout $GENOME.gvogs.tblout \
        -E 0.0001 /databases/GVDB/GVOGs/gvog.complete.hmm $GENOME.faa
hmmsearch --cpu 30 --tblout $GENOME.pfam.tblout \
        -E 0.0001 /databases/pfamA/Pfam-A.hmm $GENOME.faa

tRNAscan-SE ../$GENOME.fna > $GENOME.trna.tsv
diamond blastp -p 30 --query $GENOME.faa --db /data/nr.dmnd --out $GENOME.nr.tax -f 102 --ultra-sensitive --tmpdir /data/tmp # --taxon-exclude 693272
parse_diamond_tax.py -i $GENOME.nr.tax -o $GENOME.nr.parsed_tax
mamba activate ips
/databases/interproscan_5.62_94.0/interproscan-5.62-94.0/interproscan.sh \
        -i $GENOME.faa -o $GENOME.ips.tsv -f tsv -T /data/tmp -exclappl ProSitePatterns
mamba deactivate

# avoid self-hits for
# CroV: --taxon-exclude 693272,3044800,1513235,3047716
# ChlorellaV: --taxon-exclude 2922419

# look for ncRNAs
cmscan --noali --cut_tc -g --nohmmonly --rfam --cpu 20 --tblout $GENOME.ncrna-genes.tblout /databases/db/ncRNA-genes ../$GENOME.fna
cmscan --noali --cut_tc -g --nohmmonly --rfam --cpu 20 --tblout $GENOME.ncrna-regions.tblout /databases/db/ncRNA-regions ../$GENOME.fna

combine_annotations.py -f $GENOME.faa -ips $GENOME.ips.tsv -g $GENOME.gvogs.tblout \
                    -ga /databases/GVDB/GVOGs/gvog.complete.annot.tsv \
                    -e $GENOME.emapper.annotations -t $GENOME.nr.parsed_tax \
                    -o $GENOME.combined.annotations #.xlsx --xlsx

combine_annotations.py -f $GENOME.faa -ips $GENOME.ips.tsv -g $GENOME.gvogs.tblout \
                    -ga /databases/GVDB/GVOGs/gvog.complete.annot.tsv \
                    -e $GENOME.emapper.annotations -t $GENOME.nr.parsed_tax \
                    -o $GENOME.combined.annotations.xlsx --xlsx