micromamba activate variant


bcftools mpileup -Ou -f ChlorV-1.fna ChlorV-1.fna.sorted.bam | bcftools call -mv -Oz -o ChlorV-1.variants.vcf.gz
bcftools index ChlorV-1.variants.vcf.gz
bcftools filter -i 'QUAL>30' ChlorV-1.variants.vcf.gz -Oz -o ChlorV-1.filtered.vcf.gz