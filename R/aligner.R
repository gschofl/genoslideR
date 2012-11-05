#' Sequence alignment with fsa
#' 
#' @param seqfile Sequence files in fasta format.
#' @param outfile Multifasta alignment file.
#' @param ... Additional arguments to fsa
#' 
#' @export
fsa <- function (seqfile, outfile = "fsa.mfa", ...) {
  
  if (missing(seqfile)) {
    system("fsa --help")
    return(invisible(NULL))
  }
  
  if (!all(file.exists(seqfile)))
    stop("Can not open input files")
  
  if (length(seqfile) > 1)
    infiles <- paste(infiles, collapse=" ")
  
  SysCall("fsa", ..., stdin = seqfile, stdout = outfile, redirection = FALSE,
          style = "gnu")
}

#' Sequence alignment with mavid
#' 
#' @param seqfile Sequence files in multifasta format.
#' @param treefile Phylogenetic tree file in Newick format.
#' @param outfile Multifasta alignment file.
#' @param ... Additional arguments to mavid
#' 
#' @export
mavid <- function (seqfile, treefile = "treefile", outfile = "mavid.mfa", ...) {
  
  if (missing(seqfile)) {
    system("mavid -h")
    return(invisible(NULL))
  }
  
  if (!file.exists(seqfile)) {
    stop("Can not open input file")
  }
  
  if (!file.exists(treefile)) {
    stop("Can not open tree file")
  }
  
  SysCall("mavid", ..., stdin = paste(treefile, seqfile), stdout = outfile,
          redirection = FALSE, style = "unix")
}


#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
align_mercator_segments <- function (segments_dir, aligner = "fsa",
                                     force = TRUE, mask = TRUE) {
  
  
  segments <- normalizePath(dir(segments_dir, "^\\d+$",
                                full.names=TRUE))
  if (!force &&
        all(vapply(file.path(segments, "mavid.mfa"), file.exists, logical(1)))) {
    return(invisible(segments_dir))
  }
  
  aligner <- match.arg(aligner, c("fsa", "mavid"))
  if (aligner == "fsa") {
    aln.opts <- list(mercator="cons", exonerate=TRUE, softmasked=mask)
  } else {
    aln.opts <- list(c="cons")
  }
  ALN <- match.fun(aligner)
  
  cwd <- getwd()
  ncores <- detectCores() - 1
  mclapply(segments, function (seg) {
    setwd(seg)
    ALN(seqfile="seqs.fasta", outfile="mavid.mfa", args=aln.opts)
    setwd(cwd)
  }, mc.cores=ncores)
  setwd(cwd)
  
  return(invisible(segments_dir))
}


#' Sequence alignment with MAUVE
#' 
#' not implemented
#'  
#' @export
mauve <- function () {
  NULL
}


