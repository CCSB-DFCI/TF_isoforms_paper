
library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(upbm)
library(upbmAux)
theme_set(theme_bw())

tf <- "TBX5"

motif <- "AGGTGT"
rc_motif <- "ACACCT"

ref_condition <- "TBX5C05-REF"

datdir <- "../../../../data/internal/pbms/gpr_files"
sampdir <- "../../../../data/internal/pbms/samp_sheets"

tf_file <- paste0("../../../../data/internal/pbms/upbm_processed/", tf, "-alexa-preprocessed.rds")
tf_file

tf_ap <- readRDS(tf_file)
head(tf_ap)

alexa_pfit <- probeFit(tf_ap, stratify = "condition")
head(alexa_pfit)

pfit_dat <- broom::tidy(alexa_pfit, assay = c("beta", "sd"))
head(pfit_dat)

ggplot(pfit_dat, aes(x = beta, color = cname)) +
    geom_density() +
    scale_color_brewer("condition", palette = "Set1") +
    theme_bw() +
    ggtitle(paste(tf, "cross-replicate summarized probe intensities"))

alexa_kfit <- kmerFit(alexa_pfit, kmers = uniqueKmers(8L),
                      baseline = ref_condition)
alexa_kfit

alexa_pa <- kmerTestAffinity(alexa_kfit)
alexa_pa

pa_dat <- broom::tidy(alexa_pa, c("affinityEstimate", "affinityQ"))

pa_dat %>%
    dplyr::mutate(affinityQ = cut(affinityQ, c(1e0, 1e-1, 1e-2, 1e-4, 0),
                                  include.lowest = TRUE),
                  affinityQ = factor(affinityQ, levels = rev(levels(affinityQ)))) %>%
    ggplot(aes(x = affinityEstimate, fill = affinityQ)) +
    geom_histogram(binwidth = .1, boundary = 10, color = 'black', position = "stack", alpha = 1) +
    scale_fill_brewer("q-value", palette = "Spectral", direction = -1, drop = FALSE,
                      na.value = "black") +
    theme_bw() +
    facet_wrap( ~ cname, nrow = 1) +
    ggtitle(paste(tf,"8-mer preferential affinity"))

kfit_dat <- broom::tidy(alexa_kfit, c("affinityEstimate", "affinityVariance", "affinityQ",
                                      "contrastDifference", "contrastAverage",
                                      "contrastVariance"))

kfit_dat <- kfit_dat %>%
    dplyr::mutate(contains_any_motif = case_when(grepl(motif, seq) | grepl(rc_motif, seq) ~ paste("*", tf, "k-mer"),
                                       !grepl(motif, seq) & !grepl(rc_motif, seq) ~ "other k-mer"))


head(kfit_dat)

alexa_da <- kmerTestContrast(alexa_kfit)
alexa_da

da_dat <- broom::tidy(alexa_da, c("contrastAverage", "contrastDifference",
                                  "contrastQ"))

da_dat <- da_dat %>%
    dplyr::filter(cname != ref_condition) %>%
    dplyr::mutate(contrastQ_cut = cut(contrastQ, c(1e0, 1e-1, 1e-2, 1e-3, 0),
                                  include.lowest = TRUE),
                  contrastQ_cut = factor(contrastQ_cut, levels = rev(levels(contrastQ_cut))))

da_dat$contains_motif <- grepl(motif, da_dat$seq)
da_dat$contains_rc_motif <- grepl(rc_motif, da_dat$seq)
da_dat$contains_any_motif <- ifelse(da_dat$contains_motif | da_dat$contains_rc_motif, paste(tf, "k-mer"), "*other k-mer")

da_dat %>%
    ggplot(aes(x = contrastAverage, y = contrastDifference,
               color = contrastQ_cut, shape = contains_any_motif, 
               size = contains_any_motif, alpha = contains_any_motif)) +
    geom_point() +
    scale_color_brewer("q-value", palette = "Spectral", direction = -1, drop = FALSE,
                       na.value = "black") +
    scale_size_manual("k-mer type", values=c(1, 3)) +
    scale_shape_manual("k-mer type", values=c(16, 21)) +
    scale_alpha_manual("k-mer type", values=c(0.2, 1)) +
    geom_hline(color = 'gray', yintercept = 0) + 
    theme_bw() +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    facet_grid(. ~ cname) +
    ggtitle(paste(tf, "8-mer differential affinity"))

write.csv(kfit_dat, file=file.path("../../../../data/processed/pbms", paste0(tf, "kfit_dat.csv")), row.names=FALSE)
write.csv(da_dat, file=file.path("../../../../data/processed/pbms", paste0(tf, "da_dat.csv")), row.names=FALSE)

source("../../../pbm_utils.r")

n_to_plot <- 50

pa_dat_sort <- pa_dat[order(pa_dat$affinityQ, decreasing = FALSE), ]

pa_ref <- head(pa_dat_sort[pa_dat_sort$cname == "TBX5C05-REF",], n_to_plot)
print(nrow(pa_ref))

pa_a <- head(pa_dat_sort[pa_dat_sort$cname == "TBX5A05",], n_to_plot)
print(nrow(pa_a))

pa_b <- head(pa_dat_sort[pa_dat_sort$cname == "TBX5B05",], n_to_plot)
print(nrow(pa_b))

pa_ref_pwm <- pwm_from_ordered_list(pa_ref$seq, orientation_motif="CAC")
pa_a_pwm <- pwm_from_ordered_list(pa_a$seq, orientation_motif="CAC")
pa_b_pwm <- pwm_from_ordered_list(pa_b$seq, orientation_motif="CAC")

pwms <- list('TBX5-1'=pa_ref_pwm, 'TBX5-2'=pa_b_pwm, 'TBX5-3'=pa_a_pwm)
ggseqlogo(pwms, ncol=1)

options(repr.plot.width=1.5, repr.plot.height=1, repr.plot.pointsize=9)

ggseqlogo(pwms[[1]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../../figures/TBX5-ref_pbm_seqlogo.pdf", width=1.5, height=1)

ggseqlogo(pwms[[2]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../../figures/TBX5-2_pbm_seqlogo.pdf", width=1.5, height=1)

ggseqlogo(pwms[[3]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../../figures/TBX5-3_pbm_seqlogo.pdf", width=1.5, height=1)

options(repr.plot.width=2.5, repr.plot.height=2, repr.plot.pointsize=9)

ggseqlogo(pwms[[1]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_text(margin = margin(r=1)),
                                     axis.line.x = element_line(color="black", size=0.5),
                                     axis.line.y = element_line(color="black", size=0.5),
                                     axis.ticks.x = element_line(color="black"),
                                     axis.ticks.y = element_line(color="black"))
ggsave("../../../../figures/TBX5-ref_pbm_seqlogo_larger.pdf", width=2.5, height=2)

ggseqlogo(pwms[[2]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_text(margin = margin(r=1)),
                                     axis.line.x = element_line(color="black", size=0.5),
                                     axis.line.y = element_line(color="black", size=0.5),
                                     axis.ticks.x = element_line(color="black"),
                                     axis.ticks.y = element_line(color="black"))
ggsave("../../../../figures/TBX5-2_pbm_seqlogo_larger.pdf", width=2.5, height=2)

ggseqlogo(pwms[[3]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_text(margin = margin(r=1)),
                                     axis.line.x = element_line(color="black", size=0.5),
                                     axis.line.y = element_line(color="black", size=0.5),
                                     axis.ticks.x = element_line(color="black"),
                                     axis.ticks.y = element_line(color="black"))
ggsave("../../../../figures/TBX5-3_pbm_seqlogo_larger.pdf", width=2.5, height=2)

da_q_thresh <- 0.001

da_b <- da_dat[da_dat$contrastQ < da_q_thresh, ]
da_b <- da_b[da_b$cname=="TBX5B05", ]
print(nrow(da_b))

ref_highest <- pa_ref$seq[[1]]
ref_highest

da_b_pwm <- pwm_from_ordered_list(c(ref_highest, da_b$seq), remove_first = TRUE, orientation_motif = "CAC")

pwms <- list('TBX5-2 differential'=da_b_pwm)
ggseqlogo(pwms, ncol=1)

options(repr.plot.width=1.5, repr.plot.height=1, repr.plot.pointsize=9)

ggseqlogo(pwms[[1]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../../figures/TBX5-2_differential_pbm_seqlogo.pdf", width=1.5, height=1)

options(repr.plot.width=2.5, repr.plot.height=2, repr.plot.pointsize=9)

ggseqlogo(pwms[[1]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_text(margin = margin(r=1)),
                                     axis.line.x = element_line(color="black", size=0.5),
                                     axis.line.y = element_line(color="black", size=0.5),
                                     axis.ticks.x = element_line(color="black"),
                                     axis.ticks.y = element_line(color="black"))
ggsave("../../../../figures/TBX5-2_differential_pbm_seqlogo_larger.pdf", width=2.5, height=2)
