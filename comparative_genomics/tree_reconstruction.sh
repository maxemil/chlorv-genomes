git clone https://github.com/maxemil/gv-phylo.git
git clone https://github.com/maxemil/misc-scripts.git

nextflow run viral-evolution/gv-phylo.nf \
    --proteomes "proteomes/*" \
    --output_folder "ChlorVs_phylo" \
    --diamond_db /databases/gv-phylo/gvog_markers_v2.dmnd \
    --seeds "/databases/gv-phylo/marker_seeds_gvdb2/*.faa" \
    --gvdb_tsv /databases/GVDB-2.0/GVDB-2.0.csv \
    --taxa_colors /databases/gv-phylo/taxa_colors.txt \
    --selectors "Mimiviridae" \
    --subgroup "Imitervirales" \
    --treemode "iqtree-fast" \
    --marker_selection 'all' \
    -resume -w /data/work

# check individual markers and redo alignment if necessary
mafft --globalpair --maxiterate 1000 --thread 20 $faa > $marker.mafft
trimal -gt 0.1 -in $marker.mafft -out $marker.trimal    
iqtree -s $marker.trimal --ufboot 1000 -m MFP -mset LG -pre $marker.iqtree -nt 20

misc-scripts/concatenate.py \
    -out alignments_manual_concat/ChlorVs_phylo.aln \
    -t 3 \
    -sep '..' \
    alignments_manual/trimal/*

iqtree -s ChlorVs_phylo.aln \
    -m MFP \
    --mset LG,LG+C10,LG+C20,LG+C30,LG+C40,LG4X \
    --ufboot 1000 \
    --prefix ChlorVs_phylo.iq \
    -T 40 \
    -redo
