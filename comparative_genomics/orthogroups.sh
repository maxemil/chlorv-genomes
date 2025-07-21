cat Aliimimivirinae/* proteomes/* > Aliimimivirinae_wt_isolates.faa
# mmseqs easy-cluster Aliimimivirinae_wt_isolates.faa ChlorVs_cluster /data/tmp
orthofinder -f Aliimimivirinae_wt_isolates -t 40 -n ChlorVs_orthos -os -M msa


cat proteomes/* > Aliimimivirinae_only_isolates.faa
cat gff/* > Aliimimivirinae_only_isolates.gff
mmseqs easy-search --threads 40 -s 7.5 -e 1e-5 \
        Aliimimivirinae_only_isolates.faa \
        Aliimimivirinae_only_isolates.faa \
        Aliimimivirinae_only_isolates.o6 \
        /data/tmp