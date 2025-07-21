library(gggenomes)
library(ggplot2)
library(gggenes)   
library(data.table)
library(dplyr)

genes <- read_seqs('selected_methylases.faa')
annots <- fread('selected_methylases.annot', col.names=c('seqid', 'database', 'accession', 'description', 'tag', 'start', 'stop', 'evalue'))
annots <- annots %>% rename('end'='stop', 'seq_id'='seqid')
annots <- annots %>% mutate(across(evalue, as.double))

g = gggenomes(seqs=genes, genes=annots) +
    geom_seq(size=2) +
    geom_gene(size=6, shape=0, position=position_nudge(y=0.12), aes(fill=accession)) + 
    geom_gene_tag(aes(label=description), check_overlap=FALSE, position=position_nudge(y=0.2), size=9, angle=0) + 
    geom_seq_label() + 
    scale_x_continuous(name='position', breaks=c(0,100,200,300,365,400,403,500,589,600,700,800,895), limits=c(0,895))
ggsave("selected_methylases.pdf", limitsize=FALSE)
