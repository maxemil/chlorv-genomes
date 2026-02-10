micromamba activate variant

GENOME=ChlorV-1

bcftools mpileup -Ou -f $GENOME.fna $GENOME.fna.sorted.bam | bcftools call -mv -Oz -o $GENOME.variants.vcf.gz
bcftools index $GENOME.variants.vcf.gz
bcftools filter -i 'QUAL>30' $GENOME.variants.vcf.gz -Oz -o $GENOME.filtered.vcf.gz