#!/bin/bash

dorado -vv 2> $out_dir/$run.dorado.log
dorado duplex -v $model $pod5_dir $dorado_options 2>> $out_dir/$run.dorado.log > $out_dir/$run.bam

samtools view -O fastq -d dx:0 $out_dir/"$run".bam | pigz > $out_dir/"$run".simplex.fastq.gz
samtools view -O fastq -d dx:1 $out_dir/"$run".bam | pigz > $out_dir/"$run".duplex.fastq.gz

porechop -i $out_dir/"$run".simplex.fastq.gz -o $out_dir/"$run".simplex.trimmed.fastq.gz \
                --format fastq.gz --threads 20 &> $out_dir/"$run".simplex.porechop.log
porechop -i $out_dir/"$run".duplex.fastq.gz -o $out_dir/"$run".duplex.trimmed.fastq.gz \
                --format fastq.gz --threads 20 &> $out_dir/"$run".duplex.porechop.log