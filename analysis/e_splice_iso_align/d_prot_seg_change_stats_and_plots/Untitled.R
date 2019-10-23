library(ggplot2)
library(tidyverse)
theme_set(theme_bw())



#### plotting insertion stats
dfi = read.delim('../d_insertion_deletion_lengths_insertions.tsv', header=TRUE, na.strings=c("","NA"))
# remove identical segment changes within the same gene, nr=non-redundant (sometimes all iso. of a gene has same indel change)
dfinr = dfi %>% distinct(gene, seg1, seg2, is_ragged, seg_len, .keep_all=TRUE)

# plot insertion lengths distr - original dataset
ggplot(dfi, aes(x=seg_len)) + geom_histogram(bins=70) + xlab('Insertion length (AAs)') + ylab('Count') + ggtitle('Distr. of insertion lengths')
ggsave('insertion_seg_lengths_all.pdf', width=3, height=2.5)
# zoomed in
ggplot(dfi, aes(x=seg_len)) + geom_histogram(binwidth=4) + xlab('Insertion length (AAs)') + ylab('Count') + ggtitle('Distr. of insertion lengths') + xlim(0, 200)
ggsave('insertion_seg_lengths_all_zoomed.pdf', width=3, height=2.5)

# plot insertion lengths distr - non-redundant within gene dataset
ggplot(dfinr, aes(x=seg_len)) + geom_histogram(bins=70) + xlab('Insertion length (AAs)') + ylab('Count') + ggtitle('Distr. of insertion lengths')
ggsave('insertion_seg_lengths_nonredund.pdf', width=3, height=2.5)
# zoomed in
ggplot(dfinr, aes(x=seg_len)) + geom_histogram(binwidth=4) + xlab('Insertion length (AAs)') + ylab('Count') + ggtitle('Distr. of insertion lengths') + xlim(0, 200)
ggsave('insertion_seg_lengths_nonredund_zoomed.pdf', width=3, height=2.5)

# plot insertion lengths distr - nonredund, facet by is_ragged
dfinr$is_ragged = factor(dfinr$is_ragged, labels=c('Clean', 'Ragged edge'))
ggplot(dfinr, aes(x=seg_len)) + geom_histogram() + xlab('Insertion length (AAs)') + ylab('Count') + ggtitle('Distr. of insertion lengths') + xlim(0, 200) + facet_grid(.~is_ragged)
ggsave('insertion_seg_lengths_nonredund_zoomed_facet_ragged.pdf', width=6, height=2.5)

#### plotting deletion stats
dfd = read.delim('../d_insertion_deletion_lengths_deletions.tsv', header=TRUE, na.strings=c("","NA"))
# remove identical segment changes within the same gene, nr=non-redundant (sometimes all iso. of a gene has same indel change)
dfdnr = dfd %>% distinct(gene, seg1, seg2, is_ragged, seg_len, .keep_all=TRUE)

# plot deletion lengths distr - original dataset
ggplot(dfd, aes(x=seg_len)) + geom_histogram(bins=70) + xlab('Deletion length (AAs)') + ylab('Count') + ggtitle('Distr. of deletion lengths')
ggsave('deletion_seg_lengths_all.pdf', width=3, height=2.5)
# zoomed in
ggplot(dfd, aes(x=seg_len)) + geom_histogram(binwidth=4) + xlab('Deletion length (AAs)') + ylab('Count') + ggtitle('Distr. of deletion lengths') + xlim(0, 200)
ggsave('deletion_seg_lengths_all_zoomed.pdf', width=3, height=2.5)

# plot deletion lengths distr - non-redundant within gene dataset
ggplot(dfdnr, aes(x=seg_len)) + geom_histogram(bins=70) + xlab('Deletion length (AAs)') + ylab('Count') + ggtitle('Distr. of deletion lengths')
ggsave('deletion_seg_lengths_nonredund.pdf', width=3, height=2.5)
# zoomed in
ggplot(dfdnr, aes(x=seg_len)) + geom_histogram(binwidth=4) + xlab('Deletion length (AAs)') + ylab('Count') + ggtitle('Distr. of deletion lengths') + xlim(0, 200)
ggsave('deletion_seg_lengths_nonredund_zoomed.pdf', width=3, height=2.5)

# plot Deletion lengths distr - nonredund, facet by is_ragged
dfdnr$is_ragged = factor(dfdnr$is_ragged, labels=c('Clean', 'Ragged edge'))
ggplot(dfdnr, aes(x=seg_len)) + geom_histogram() + xlab('Deletion length (AAs)') + ylab('Count') + ggtitle('Distr. of deletion lengths') + xlim(0, 200) + facet_grid(.~is_ragged)
ggsave('deletion_seg_lengths_nonredund_zoomed_facet_ragged.pdf', width=6, height=2.5)



### plotting substitution stats
df = read.table('../d_substitution_stats.tsv', header=TRUE)
colnames(df)
ggplot(df, aes(x=len_diff)) + geom_histogram(binwidth=5)

ggplot(df, aes(x=len_diff)) + geom_histogram()
ggplot(df, aes(x=perc_change)) + geom_histogram(binwidth = 1) + xlim(0, 100)

# plot ref vs alt segment length in substitutions
ggplot(df, aes(x=seg2_len, y=seg1_len)) + geom_point(size=1, alpha=0.5) + scale_x_log10() + scale_y_log10() + 
  xlab('Protein segment length - alt. iso.') + ylab('Protein segment length - ref. iso.') + ggtitle('Len. diff. between substituted prot. segments')
ggsave('compare_protein_seg_lengths_between_isoforms.pdf', width=8, height=5)

# df with clustal scores of comparing subst. segments
dfc = read.table('../d_substitution_stats_w_clustal_scores.txt', header=TRUE)
colnames(dfc)
ggplot(dfc, aes(x=perc_frameshift)) + geom_histogram(binwidth=1) + xlab('Percentage of subst. region that is frameshift') + ylab('Count') + xlim(0.00000001, 100)
ggsave('subst_region_frameshift_perc.pdf', width=4, height=3)


# perc. identity of subst. segments h=homology analysis
dfh = read.delim('../d_substitution_stats_w_clustal_seq_compare.txt', header=TRUE)
dfh = dfh %>% mutate(min_seg_len=pmin(seg1_len, seg2_len))
colnames(dfh)
ggplot(dfh, aes(x=min_seg_len, y=perc_ident)) + geom_point(size=0.5, alpha=0.3) + xlab('Length of shortest polypep. segment') + ylab('Percent identity')
ggsave('substitution_clustal_homology_len_vs_percident.pdf', width=4, height=3)




dfc_w_facet = transform(dfc, facet = cut(clustal_score2, seq(1, 200, 10)))
dfc_w_facet = transform(dfc, facet = cut(clustal_score1, seq(0, 100, 10)))
ggplot(dfc_w_facet, aes(x=seg2_len, y=seg1_len)) + geom_point(size=1, alpha=0.5) + scale_x_log10() + scale_y_log10() + coord_fixed(ratio=1) +
  xlab('Protein segment length - alt. iso.') + ylab('Protein segment length - ref. iso.') + ggtitle('Len. diff. between substituted prot. segments') +
  facet_wrap(~facet)
ggsave('compare_protein_seg_lengths_between_isoforms.pdf', width=5, height=5)

