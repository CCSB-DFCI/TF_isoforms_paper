
library(DECIPHER)
library(seqLogo) 
library(seqinr)
library(ggseqlogo)

fix_orientations <- function(oo, orientation_motif = orientation_motif) {
    oo_fixed <- list()
    for (x in 1:length(oo)) {
        flip <- grepl(orientation_motif, oo[x])
        if (flip) {
            seq <- reverseComplement(oo[x])
        } else {
            seq <- oo[x]
        }
        oo_fixed[[x]] <- as.character(seq)
    }

    oo_fixed <- DNAStringSet(unlist(oo_fixed))
    return(oo_fixed)
}

pwm_from_ordered_list <- function(seqs, rc = FALSE, remove_first = FALSE, orientation_motif = orientation_motif) {
    print(length(seqs))
    ss <- DNAStringSet(seqs)
    oo <- OrientNucleotides(ss, type="both", orientation="both")[[2]]
    
    if (orientation_motif != "none") {
        oo <- fix_orientations(oo, orientation_motif = orientation_motif)
    }
    
    aa <- AlignSeqs(oo, gapExtension=0, useStructures=FALSE)
    aj <- AdjustAlignment(aa)
    
    if ( remove_first == TRUE ) {
        n_seqs <- length(aj)
        aj <- aj[2:n_seqs]
    }
    
    pre_pwm <- consensusMatrix(aj)[c(1,2,3,4), ]
    pre_pwm <-  t(pre_pwm) + (length(seqs)-colSums(pre_pwm))/4
    pre_pwm <- t(pre_pwm/rowSums(pre_pwm))
    
    if ( rc == TRUE ) {
        reverse_pre_pwm <- pre_pwm
        reverse_pre_pwm[] <- pre_pwm[nrow(pre_pwm):1, ncol(pre_pwm):1]
        pwm <- reverse_pre_pwm
    } else {
        pwm <- pre_pwm
    }
    
    return(pwm)    
}

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
