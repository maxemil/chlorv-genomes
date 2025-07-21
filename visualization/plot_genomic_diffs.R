library(tidyverse)
library(gggenomes)
library(ggpubr)
library(ggrepel)

vseqs <- read_seqs(list.files('genomes', '*.fna$', full.names=TRUE), format='fasta')
vseqs <- vseqs %>% filter(seq_id %in% c('ChlorV-1', 'ChlorV-4', 'ChlorV-2', 'ChlorV-3'))

vgenes <- read_feats(list.files('gff', '*.gff$', full.names=TRUE), format='gff3')

prots = read_tsv('annotations/filtered.annotations', col_names=FALSE)

vgenes = merge(vgenes, prots, by.x='feat_id', by.y='X1', all.x=TRUE)

vlinks = read_tsv('Aliimimivirinae_only_isolates.filtered.links', col_names=FALSE)
colnames(vlinks) <- c('seq_id', 'seq_id2', 'start', 'end', 'start2', 'end2', 'strand')

vblast <- read_sublinks('Aliimimivirinae_only_isolates.o6') %>% filter(pident > 0.7)
vblast <- vblast %>% filter((feat_id %in% vgenes$feat_id) & (feat_id2 %in% vgenes$feat_id))
vblast <- vblast %>% filter(feat_id != vblast$feat_id2)

vgenes$match <- ifelse(vgenes$feat_id %in% vblast$feat_id, 'matched' , 'unmatched')

proteomics <- read_tsv('../01_assemblies/ChlorV-1/proteomics/ChlorV-1_proteomics_annotated.tsv')
vgenes$capsid <- ifelse(vgenes$feat_id %in% proteomics$Protein.Group, 'in_capsid' , 'not_found')

panorthologs <- read_tsv('Panorthologs.tsv', col_names=FALSE)
vgenes$match <- ifelse(vgenes$feat_id %in% panorthologs$X1, 'matched' , 'unmatched')

loci = list(
A = c('ChlorV-1..001', 'ChlorV-1..004', 'ChlorV-4..005', 'ChlorV-4..008', 'ChlorV-2..006', 'ChlorV-2..009', 'ChlorV-3..003', 'ChlorV-3..008'),
B = c('ChlorV-1..006', 'ChlorV-1..012', 'ChlorV-4..010', 'ChlorV-4..015', 'ChlorV-4..018', 'ChlorV-2..011', 'ChlorV-2..017', 'ChlorV-2..021', 'ChlorV-3..010', 'ChlorV-3..015'),
D = c('ChlorV-1..080', 'ChlorV-1..086', 'ChlorV-4..078', 'ChlorV-4..081', 'ChlorV-2..081', 'ChlorV-2..084', 'ChlorV-3..068', 'ChlorV-3..071'),
E = c('ChlorV-1..110', 'ChlorV-1..114', 'ChlorV-4..107', 'ChlorV-4..109', 'ChlorV-2..114', 'ChlorV-2..116', 'ChlorV-3..100', 'ChlorV-3..105'),
F = c('ChlorV-1..311', 'ChlorV-1..317', 'ChlorV-4..308', 'ChlorV-4..314', 'ChlorV-2..324', 'ChlorV-2..327', 'ChlorV-2..331', 'ChlorV-3..306', 'ChlorV-3..313'),
G = c('ChlorV-1..370', 'ChlorV-1..373', 'ChlorV-1..376', 'ChlorV-4..365', 'ChlorV-4..368', 'ChlorV-2..380', 'ChlorV-2..383', 'ChlorV-3..360', 'ChlorV-3..367'),
locus2 = c('ChlorV-1..041', 'ChlorV-1..047', 'ChlorV-1..054', 'ChlorV-4..043', 'ChlorV-4..045', 'ChlorV-4..049', 'ChlorV-2..046', 'ChlorV-2..052', 'ChlorV-3..038', 'ChlorV-3..044'),
locus4 = c('ChlorV-1..170', 'ChlorV-1..172', 'ChlorV-4..163', 'ChlorV-4..164', 'ChlorV-2..173', 'ChlorV-2..174', 'ChlorV-3..167', 'ChlorV-3..168'),
locus6 = c('ChlorV-1..398', 'ChlorV-1..399', 'ChlorV-4..389', 'ChlorV-4..392', 'ChlorV-2..404', 'ChlorV-2..408', 'ChlorV-3..388', 'ChlorV-3..390'),
locus8 = c('ChlorV-1..092', 'ChlorV-1..095', 'ChlorV-4..086', 'ChlorV-4..089', 'ChlorV-2..089', 'ChlorV-2..095', 'ChlorV-3..076', 'ChlorV-3..081'),
locus9 = c('ChlorV-1..077', 'ChlorV-1..080', 'ChlorV-4..077', 'ChlorV-4..078', 'ChlorV-2..080', 'ChlorV-2..081', 'ChlorV-3..067', 'ChlorV-3..068')
)

locus_locs = tibble('seq_id'=character(), 'feat_id'=character(), 'start'=numeric(), 'end'=numeric())
for (locus in names(loci)) {
    genes_p <- vgenes %>% filter(feat_id %in% loci[[locus]]) %>% filter(seq_id == 'ChlorV-1')
    locus_locs = bind_rows(locus_locs, tibble(seq_id='ChlorV-1', feat_id=locus, start=min(genes_p$start), end=max(genes_p$end)))
    genes_p <- vgenes %>% filter(feat_id %in% loci[[locus]]) %>% filter(seq_id == 'ChlorV-2')
    locus_locs = bind_rows(locus_locs, tibble(seq_id='ChlorV-2', feat_id=locus, start=min(genes_p$start), end=max(genes_p$end)))
    genes_p <- vgenes %>% filter(feat_id %in% loci[[locus]]) %>% filter(seq_id == 'ChlorV-3')
    locus_locs = bind_rows(locus_locs, tibble(seq_id='ChlorV-3', feat_id=locus, start=min(genes_p$start), end=max(genes_p$end)))
    genes_p <- vgenes %>% filter(feat_id %in% loci[[locus]]) %>% filter(seq_id == 'ChlorV-4')
    locus_locs = bind_rows(locus_locs, tibble(seq_id='ChlorV-4', feat_id=locus, start=min(genes_p$start), end=max(genes_p$end)))
}

locus_locs <- locus_locs %>% filter(feat_id %in% c('A', 'B', 'D', 'E', 'F', 'G'))

pbase = gggenomes(seqs=vseqs, genes=vgenes, feats=(locus_locs)) %>%
    pick('ChlorV-1', 'ChlorV-3', 'ChlorV-2', 'ChlorV-4') %>%
    add_sublinks(vblast) %>%
    sync() +
    geom_seq() +
    geom_gene(size=6, data = genes(capsid == "not_found"), aes(fill=match), position=position_nudge(y = 0.12)) +
    geom_gene(size=6, data = genes(capsid == "in_capsid"), aes(fill=match, stroke=1, linetype='dashed'), position=position_nudge(y = 0.12)) +
    geom_link(fill='#878787', position=position_nudge(y = 0.12), color = rgb(0,0,0, alpha=0)) +
    scale_fill_manual(values = c('#E69F00', '#56B4E9'))
    # scale_fill_viridis(discrete = TRUE, option='G') #+
# pannot = pbase + geom_bin_label(expand_left = .01, size=16, hjust=-0.45, vjust=0.3) +
#                  geom_gene_tag(aes(label=X2), nudge_y=0.2, check_overlap=TRUE)
# ggsave("ChlorV1-4_gggenome_filtlinks_annot.pdf", width=100, height=10, limitsize=FALSE)

# pid = pbase + geom_bin_label(expand_left = .01, size=16, hjust=-0.45, vjust=0.3) +
#                  geom_gene_tag(aes(label=feat_id), nudge_y=0.2, check_overlap=TRUE)
# ggsave("ChlorV1-4_gggenome_filtlinks_ids.pdf", width=100, height=10, limitsize=FALSE)


lplots = list()
for (locus in names(loci)) {
    p2 <- pbase %>% focus(feat_id %in% loci[[locus]], .track_id = genes, .overhang='drop', .locus_id = gsub("ChlorV-(\\d)", "\\1", seq_id), .expand = 0, .locus_bin='locus', .max_dist=100000)+ 
                    # geom_text_repel(aes(x=(x+xend)/2, y=y, label=X2), size=9, min.segment.length = unit(0, 'lines'), data = genes(), stat='identity', max.overlaps=40) + 
                    geom_gene_tag(aes(label=X2), nudge_y=0.2, check_overlap=FALSE, size=9, position=position_dodge(width = 0.5), angle=45) +
                    geom_bin_label(expand_left = 0, size=12, hjust=.4, vjust=0.1)
    lplots[[locus]] = p2
    ggsave(paste("ChlorV1-4_gggenome_filtlinks_selection_", locus, ".pdf", sep=''), plot=p2, width=10, height=5, limitsize=FALSE)
}


pall <- pbase + geom_bin_label(expand_left = 0, size=12, hjust=.4, vjust=0.1) + 
                geom_feat(color = "black", position=position_nudge(y=0.5)) + 
                geom_feat_text(aes(label = feat_id), position=position_nudge(y=0.5), vjust = -0.5, size=9)

parrange = ggarrange(
#   ggarrange(lplots[['locus10']], lplots[['locus3']], lplots[['locus7']], ncol = 3, labels = c("A", "B", "C"), legend='none'), 
  ggarrange(lplots[['A']], lplots[['B']], ncol = 2, labels = c("A", "B"), legend='none'), 
  pall,
  ggarrange(lplots[['D']], lplots[['E']], ncol = 2, labels = c("D", "E"), legend='none'), 
  ggarrange(lplots[['F']], lplots[['G']], ncol = 2, labels = c("F", "G"), legend='none'), 
#   ggarrange(lplots[['locus1']], lplots[['locus9']], lplots[['locus5']], ncol = 3, labels = c("E", "F"), legend='none'), 
  nrow = 4, 
  labels = c("A","C","D","F"),
  legend='none'
  )
ggsave("ChlorV1-4_gggenome_filtlinks_arranged.pdf", plot=parrange, width=30, height=25, limitsize=FALSE)