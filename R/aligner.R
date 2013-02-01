#' Sequence alignment with fsa
#' 
#' @param seqfile Sequence files in fasta format.
#' @param outfile Multifasta alignment file.
#' @param opts a named list of options for fsa.
#' @param ... Named values interpreted as options for fsa.
#' 
#' @export
fsa <- function (seqfile, outfile = "fsa.mfa", opts = list(), ...) {
  
  if (missing(seqfile)) {
    system("fsa --help")
    return(invisible(NULL))
  }
  
  if (!all(file.exists(seqfile)))
    stop("Can not open input files")
  
  if (length(seqfile) > 1)
    infiles <- paste(infiles, collapse=" ")
  
  args <- merge_list(opts, list(...))
  SysCall("fsa", args = args, stdin = seqfile, stdout = outfile,
          redirection = FALSE, style = "gnu")
}


#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
align_mercator_segments <- function (seg_dir, force, opts = list()) {
  
  segments <- normalizePath(dir(seg_dir, "^\\d+$", full.names=TRUE))
  if (!force &&
      all(vapply(file.path(segments, "mavid.mfa"), file.exists, logical(1)))) {
    return(invisible(seg_dir))
  }
  
  aln.opts <- merge_list(list(mercator="cons"), opts)
  cwd <- getwd()
  ncores <- detectCores() - 1
  mclapply(segments, function (seg) {
    setwd(seg)
    fsa(seqfile="seqs.fasta", outfile="mavid.mfa", opts=aln.opts)
    setwd(cwd)
  }, mc.cores=ncores)
  setwd(cwd)
  
  return(invisible(seg_dir))
}

