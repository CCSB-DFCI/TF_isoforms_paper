library(tidyverse)
library(cowplot)
library(ggrepel)
library(readxl)
library(upbm)
library(upbmAux)
theme_set(theme_bw())

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r * 10)
}

# writing the design files for GEO submission
write.table(pbm_8x60k_v1@design, file=file.path("../../../data/internal/pbms/gpr_geo/pbm_8x60k_v1.txt"), 
            row.names=FALSE,
            quote=FALSE, sep="\t")
write.table(broom::tidy(refcy3_8x60k_v1), file=file.path("../../../data/internal/pbms/gpr_geo/refcy3_8x60k_v1.txt"), 
            row.names=FALSE,
            quote=FALSE, sep="\t")

head(refcy3_8x60k_v1)

tf <- "CREB1"

motif <- "CGTCA"
rc_motif <- "TGACG"

ref_condition <- "CREB1-ref"

datdir <- "../../../data/internal/pbms/gpr_geo"
sampdir <- "../../../data/internal/pbms/gpr_geo"

alexa_samples <- read_csv(file.path(sampdir, paste0(tf, "-alexa.csv")), col_types = cols())
alexa_samples

cy3_samples <- read_csv(file.path(sampdir, paste0(tf, "-cy3.csv")), col_types = cols())
cy3_samples

# properly prefix file paths for data
alexa_samples <- dplyr::mutate(alexa_samples, gpr = file.path(datdir, gpr))
cy3_samples <- dplyr::mutate(cy3_samples, gpr = file.path(datdir, gpr))

tf_a <- gpr2PBMExperiment(alexa_samples, probes = pbm_8x60k_v1)
tf_c <- gpr2PBMExperiment(cy3_samples, probes = pbm_8x60k_v1)

tf_a_df <- broom::tidy(tf_a, long = TRUE) 
tf_c_df <- broom::tidy(tf_c, long = TRUE) 

tf_a_df %>%
    dplyr::group_by(id, condition, pmt) %>%
    dplyr::summarize(nsat = sum(log2(fore) > 15.5 |
                                log2(fore) < 4, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(pmt, nsat) %>%
    print(n = 50)

tf_a_df %>%
    dplyr::group_by(id, condition, pmt) %>%
    dplyr::summarize(nna = sum(is.na(fore))) %>%
    dplyr::ungroup() %>%
    tidyr::spread(pmt, nna) %>%
    print(n = 50)

chosen_pmt <- 500

tf_a_pmt <- tf_a[, colData(tf_a)$pmt == chosen_pmt]

tf_a2 <- broom::tidy(tf_a_pmt, "fore", long = TRUE)

ggplot(tf_a2, aes(x = fore, group = cname, color = condition)) +
    geom_density() +
    scale_x_continuous("probe intensity (log2-scale)", trans = "log2") + 
    facet_grid(id ~ ., labeller = label_both) + 
    theme_bw() +
    ggtitle("Distributions of probe intensities")

tf_a_pmt <- backgroundSubtract(tf_a_pmt, assay = "fore", assayb = "back")
tf_c <- backgroundSubtract(tf_c, assay = "fore", assayb = "back")

tf_a2 <- broom::tidy(tf_a_pmt, "fore", long = TRUE)

ggplot(tf_a2, aes(x = fore, group = cname, color = condition)) +
    geom_density() +
    scale_x_continuous("BSI (log2-scale)", trans = "log2") + 
    facet_grid(id ~ ., labeller = label_both) + 
    theme_bw() +
    ggtitle("Distributions of background subtracted intensities")

tf_c_e <- cy3FitEmpirical(tf_c, refcy3_8x60k_v1)

tf_a_pmt <- cy3Normalize(pe = tf_a_pmt, cy3pe = tf_c_e, match_by = "id_idx")

tf_a_pmt <- spatiallyAdjust(tf_a_pmt)

tf_a_pmt <- normalizeWithinReplicates(tf_a_pmt)

broom::tidy(tf_a_pmt, long = TRUE) %>%
    ggplot(aes(x = normalized, group = cname, color = condition)) +
    geom_density() +
    scale_x_continuous("normalized intensity (log2-scale)", trans = "log2") + 
    facet_grid(id ~ ., labeller = label_both) + 
    theme_bw() +
    ggtitle("Distributions of within replicate normalized intensities")

tf_a_pmt <- normalizeAcrossReplicates(tf_a_pmt)

broom::tidy(tf_a_pmt, long = TRUE) %>%
    ggplot(aes(x = normalized, group = cname, color = condition)) +
    geom_density() +
    scale_x_continuous("normalized intensity (log2-scale)", trans = "log2") + 
    facet_grid(id ~ ., labeller = label_both) + 
    theme_bw() +
    ggtitle("Distributions of across replicate normalized intensities")

tf_ap_df <- broom::tidy(tf_a_pmt, long = TRUE)

l <- length(unique(tf_ap_df$probeID))

tf_ap_df %>%
    dplyr::group_by(id, condition) %>%
    dplyr::summarize(nna = sum(is.na(normalized)), percent_filtered_out = 100*nna/l) %>%
    dplyr::filter(percent_filtered_out<20) %>%
    dplyr::ungroup() %>%
    print(n = 100)

tf_ap8 <- summarizeKmers(pe = tf_a_pmt,
                         metrics = "median")

tf_ap8_df <- broom::tidy(tf_ap8, long = TRUE)

tf_ap8_df %>%
    dplyr::filter(pmt == chosen_pmt) %>%
    dplyr::filter(condition == "CREB1-ref") %>%
    dplyr::select(kmer, condition, id_idx, medianIntensity) %>%
    dplyr::mutate(medianIntensity=log2(medianIntensity)) %>%
    tidyr::spread(id_idx, medianIntensity) %>%
    dplyr::select(-kmer, -condition) %>%
    as.matrix() %>%
    pairs(col = rgb(0, 0, 0, 1/10), cex = .4, upper.panel = panel.cor,
          main = "Cross-replicate comparison (8-mer median intensities): CREB1-ref")

tf_ap8_df %>%
    dplyr::filter(pmt == chosen_pmt) %>%
    dplyr::filter(condition == "CREB1-alt") %>%
    dplyr::select(kmer, condition, id_idx, medianIntensity) %>%
    dplyr::mutate(medianIntensity=log2(medianIntensity)) %>%
    tidyr::spread(id_idx, medianIntensity) %>%
    dplyr::select(-kmer, -condition) %>%
    as.matrix() %>%
    pairs(col = rgb(0, 0, 0, 1/10), cex = .4, upper.panel = panel.cor,
          main = "Cross-replicate comparison (8-mer median intensities): CREB1-alt")

alexa_pfit <- probeFit(tf_a_pmt, stratify = "condition")
head(alexa_pfit)

pfit_dat <- broom::tidy(alexa_pfit, assay = c("beta", "sd"))
head(pfit_dat)

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

kfit_dat$contains_motif <- grepl(motif, kfit_dat$seq)
kfit_dat$contains_rc_motif <- grepl(rc_motif, kfit_dat$seq)
kfit_dat$contains_any_motif <- ifelse(kfit_dat$contains_motif | kfit_dat$contains_rc_motif, paste(tf, "k-mer"), "*other k-mer")
kfit_dat[kfit_dat$contains_any_motif == paste(tf, "k-mer"), ]
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

source("../../pbm_utils.r")

n_to_plot <- 50

pa_dat_sort <- pa_dat[order(pa_dat$affinityQ, decreasing = FALSE), ]

pa_ref <- head(pa_dat_sort[pa_dat_sort$cname == "CREB1-ref",], n_to_plot)
print(nrow(pa_ref))

pa_a <- head(pa_dat_sort[pa_dat_sort$cname == "CREB1-alt",], n_to_plot)
print(nrow(pa_a))

pa_ref_pwm <- pwm_from_ordered_list(pa_ref$seq, orientation_motif="GTCA")
pa_a_pwm <- pwm_from_ordered_list(pa_a$seq, orientation_motif="GTCA")

pwms <- list('CREB1-ref'=pa_ref_pwm, 'CREB1-alt'=pa_a_pwm)

options(repr.plot.width=1.5, repr.plot.height=1, repr.plot.pointsize=9)

ggseqlogo(pwms[[1]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../figures/fig3/CREB1-ref_pbm_seqlogo.pdf", width=1.5, height=1)

ggseqlogo(pwms[[2]], ncol=1) + theme(axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_blank())
ggsave("../../../figures/fig3/CREB1-alt_pbm_seqlogo.pdf", width=1.5, height=1)

write.csv(kfit_dat, file=file.path("../../../data/processed/pbms", paste0(tf, "kfit_dat.csv")), row.names=FALSE)
write.csv(da_dat, file=file.path("../../../data/processed/pbms", paste0(tf, "da_dat.csv")), row.names=FALSE)

## for GEO, need affinityEstimate per protein
head(pa_dat)

geo <- pa_dat %>% reshape2::dcast(seq ~ cname, value.var = "affinityEstimate")
head(geo)

write.csv(geo, file=file.path("../../../data/processed/pbms", paste0(tf, "-geo.csv")), row.names=FALSE,
          quote=FALSE)

## for supp, want affinity estimates + da q-value
supp <- da_dat %>% 
    select(seq, cname, contrastAverage, contrastDifference, contrastQ)

supp <- full_join(supp, geo, by="seq") %>%
    select(seq, "CREB1-ref", "CREB1-alt", contrastAverage, contrastDifference, contrastQ)

write.table(supp, file="../../../supp/SuppTable_CREB1-PBM.txt", row.names=FALSE, quote=FALSE, sep="\t")
