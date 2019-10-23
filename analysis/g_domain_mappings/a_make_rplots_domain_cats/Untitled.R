library(ggplot2)
library(tidyverse)
theme_set(theme_bw())


df = read.table('../a_domain_per_res_prot_block_change.txt', header=TRUE, na.strings='-')
colnames(df)

# all-sub. domain blocks
dfs = df %>% filter(sub_perc_diff != 'NA')
dfs = dfs %>% mutate(is_zf = case_when(grepl('zf', dom_name, ignore.case=TRUE) ~ 'is_ZF', TRUE ~ 'not_ZF'))
dfs$same_dom_type_in_alt_seg = factor(dfs$same_dom_type_in_alt_seg, levels=c(0, 1), labels=c('Same domain-type missing in alt. iso.', 'Same domain-type in alt. iso.'))

# distr. of subst. seg. perc. diff. in len, facet by if domain mapped in alt. iso.
ggplot(dfs, aes(x=sub_perc_diff)) + geom_histogram() + facet_grid(same_dom_type_in_alt_seg~., scales='free_y')
ggsave('distr_perc_diff_sub_seg_lens_facet_dom_map_in_alt.pdf', width=5, height=5)

# distr. of subst. seg. perc. diff. in len, colored by if is ZF
ggplot(dfs, aes(x=sub_perc_diff, fill=is_zf)) + geom_histogram(position='identity', alpha=0.5, binwidth=1) + facet_grid(same_dom_type_in_alt_seg~., scales='free_y') + xlab('Subst. domain segment, perc. diff b/t ref/alt')
ggsave('distr_perc_diff_sub_seg_lens_facet_dom_map_in_alt_fill_zf.pdf', width=5, height=5)

