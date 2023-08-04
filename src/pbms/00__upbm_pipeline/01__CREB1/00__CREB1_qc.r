
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

tf <- "CREB1"

motif <- "CGTCA"
rc_motif <- "TGACG"

ref_condition <- "CREB1-ref"

datdir <- "../../../../data/internal/pbms/gpr_files"
sampdir <- "../../../../data/internal/pbms/samp_sheets"

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

chosen_pmt <- 500

tf_a_df %>%
    dplyr::group_by(id, condition, pmt) %>%
    dplyr::summarize(nna = sum(is.na(fore))) %>%
    dplyr::ungroup() %>%
    tidyr::spread(pmt, nna) %>%
    print(n = 50)

options(repr.plot.width=10, repr.plot.height=8)

tf_a_df %>%
    dplyr::filter(pmt > 300) %>%
    dplyr::mutate(tf = paste0(tf, " (", id, ")")) %>%
    ggplot(aes(x = fore)) +
    geom_density(aes(fill = condition), alpha = 1/2) +
    facet_grid(pmt ~ tf) +
    scale_fill_brewer(palette = "Set1") + 
    scale_x_continuous("probe intensity (log2 scale)", trans = "log2") +
    ggtitle("Comparison of probe intensity distributions")

options(repr.plot.width=16, repr.plot.height=8)

tf_a_df %>%
    dplyr::filter(pmt > 300) %>%
    dplyr::mutate(condition = paste0("(array:", id_idx, ")\n", condition)) %>%
    ggplot(aes(x = Column, y = Row, fill = log2(fore))) +
    geom_tile() +
    scale_fill_distiller(palette = "Spectral") +
    scale_x_continuous(NULL, expand = c(0, 0)) +
    scale_y_continuous(NULL, expand = c(0, 0)) +
    facet_grid(pmt ~ condition) +
    ggtitle("Comparison of probe intensities") +
    theme(strip.text.x = element_text(size = 6))

tf_c_df %>%
    dplyr::mutate(fore = scales::squish(fore, range = c(2^9, 2^14))) %>%
    dplyr::mutate(condition = paste0("(array:", id_idx, ")")) %>%
    ggplot(aes(x = Column, y = Row, fill = log2(fore))) +
    geom_tile() +
    scale_fill_distiller("log2(fore)\n[9,14]", palette = "Spectral") +
    scale_x_continuous(NULL, expand = c(0, 0)) +
    scale_y_continuous(NULL, expand = c(0, 0)) +
    facet_wrap(~ condition, nrow = 2) +
    ggtitle("Comparison of probe intensities (Cy3) - truncated limits") +
    theme(strip.text.x = element_text(size = 8))

options(repr.plot.width=8, repr.plot.height=8)

tf_a_df %>%
    dplyr::filter(pmt == chosen_pmt) %>%
    dplyr::filter(condition == "CREB1-ref") %>%
    dplyr::select(Column, Row, condition, id_idx, fore) %>%
    dplyr::mutate(fore=log2(fore)) %>%
    tidyr::spread(id_idx, fore) %>%
    dplyr::select(-Column, -Row, -condition) %>%
    as.matrix() %>%
    pairs(col = rgb(0, 0, 0, 1/10), cex = .4, upper.panel = panel.cor,
          main = "Alexa488 probe intensities: CREB1-ref")

tf_a_df %>%
    dplyr::filter(pmt == chosen_pmt) %>%
    dplyr::filter(condition == "CREB1-alt") %>%
    dplyr::select(Column, Row, condition, id_idx, fore) %>%
    dplyr::mutate(fore=log2(fore)) %>%
    tidyr::spread(id_idx, fore) %>%
    dplyr::select(-Column, -Row, -condition) %>%
    as.matrix() %>%
    pairs(col = rgb(0, 0, 0, 1/10), cex = .4, upper.panel = panel.cor,
          main = "Alexa488 probe intensities: CREB1-alt")

tf_a_pmt <- tf_a[, colData(tf_a)$pmt == chosen_pmt]

# given that array 462_2 has generally low probe-level correlations, filter out
tf_a_pmt <- tf_a_pmt[, colData(tf_a_pmt)$id_idx != "462_2"]

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

tf_file <- paste0(tf, "-alexa-preprocessed.rds")
saveRDS(tf_a_pmt, file.path("../../../../data/internal/pbms/upbm_processed", tf_file))

tf_ap_df <- broom::tidy(tf_a_pmt, long = TRUE)

l <- length(unique(tf_ap_df$probeID))

tf_ap_df %>%
    dplyr::group_by(id, condition) %>%
    dplyr::summarize(nna = sum(is.na(normalized)), percent_filtered_out = 100*nna/l) %>%
    dplyr::filter(percent_filtered_out<20) %>%
    dplyr::ungroup() %>%
    print(n = 100)

tf_file <- paste0(tf, "-alexa-preprocessed8.rds")
tf_ap8 <- summarizeKmers(pe = tf_a_pmt,
                         metrics = "median")
saveRDS(tf_ap8, file.path(file.path("../../../../data/internal/pbms/upbm_processed", tf_file)))

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

tf_ap8_df %>%
    dplyr::mutate(tf = paste0(tf, " (", id, ")")) %>%
    ggplot(aes(x = medianIntensity)) +
    geom_density(aes(group = condition, color = condition), alpha = 1/2) +
    facet_wrap(~ tf, nrow = 2, dir = 'v') +
    scale_fill_brewer(palette = "Set1") + 
    scale_x_continuous("8-mer median intensity (log2 scale)", trans = "log2") +
    ggtitle("Comparison of 8-mer median intensity distributions") +
    coord_cartesian(xlim = 2^c(7, 12))
