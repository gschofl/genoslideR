## Die Funktion stellt zwei Parameter dar, Eta und Eta_e.
## Eta ist die Mutationen auf dem Teil Genomgebiet.
## Eta_e ist die Mutationen auf den extern Zweig im Phylogie-Baum.
Eta <- function (dna, outgroup = NULL) {
  if (!is.null(outgroup) && outgroup > length(dna)) {
    stop("Outgroup out of range")
  }
  codes <- c(A=1L, C=2L, G=4L, T=8L)
  cm <- .Call2("XStringSet_consensus_matrix", dna, shift = 0L, 
                width = NULL, baseOnly = TRUE, codes,
               PACKAGE = "Biostrings")
  gaps <- cm[5, ] > 0L
  if (all(gaps)) {
    return( c(Eta = NA_integer_, Eta_e = NA_integer_) )
  }
  cm <- cm[1:4, !gaps]
  if (is.null(outgroup)) {
    if (sum(!gaps) == 1L) {
      Eta <- sum(cm != 0) - 1
      if (any(cm == 1)) {
        Eta_e <- sum(cm == 1)
      } else {
        Eta_e <- 0
      }
    } else {
      Eta <- sum(cm != 0) - dim(cm)[2L]
      Eta_e <- sum(cm[, which(colSums(cm == 1) != 0)] == 1)
    }
  } else {
    cm_in <- .Call2("XStringSet_consensus_matrix", dna[-outgroup],
                    shift = 0L,  width = NULL, baseOnly = TRUE, codes,
                    PACKAGE = "Biostrings")
    cm_in <- cm_in[1:4, !gaps]
    cm_out <- .Call2("XStringSet_consensus_matrix", dna[outgroup],
                     shift = 0L,  width = NULL, baseOnly = TRUE, codes,
                     PACKAGE = "Biostrings")
    cm_out  <- cm_out[1:4, !gaps]
    if (sum(!gaps) == 1L) {
      Eta <- sum(cm_in != 0) - 1
      if (any(cm_in[cm_out != 1] == 1)) {
        Eta_e <- sum(cm_in[cm_out != 1] == 1)
      } else {
        Eta_e <- 0
      }
    } else {
      Eta <- sum(cm_in != 0) - dim(cm)[2L]
      if (any(cm_in[cm_out != 1] == 1)) {
        Eta_e <- sum(cm_in[cm_out != 1L] == 1L)
      } else {
        Eta_e <- 0
      }
    }
  }
  
  return( c(Eta = Eta, Eta_e = Eta_e) )
}
