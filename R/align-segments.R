#' Align orthologous genome segments with \code{fsa} 
#' 
#' @param seg_dir Path to mercator segments.
#' @param fsa.opts a named list of options for \code{\link{fsa}}. The default
#' options are optimisiing fsa for long sequence alignment.
#' @param ... Named values interpreted as options for fsa.
#' @param skip.completed if \code{TRUE}, don't realign an existing alignment.
#' @param ncores Number of cores to utilised for parallelisation
#' 
#' @export
alignSegments <- function (seg_dir,
                           fsa.opts = list(anchored = TRUE,
                                           translated = TRUE,
                                           exonerate = TRUE,
                                           softmasked = TRUE),
                           ...,
                           skip.completed = TRUE,
                           ncores = detectCores() - 1) {
  
  if (is_segments_dir(seg_dir)) {
    seg_dir <- segments_dir(seg_dir)
  } else {
    stop("No valid path to mercator segments")
  }
  
  fsa.opts <- merge_list(fsa.opts, list(...))
  seg_dir <- fsaAlignSegmentDirs(initdir=seg_dir, seqfile="seqs.fasta",
                                 outfile="fsa.mfa", constraints="cons",
                                 skip.completed=skip.completed,
                                 fsa.opts=fsa.opts, ncores=ncores)

  message("Pass the path to the aligned mercator segments to the 'annotatedAlignment()' constructor to create an 'annotatedAlignment' object")
  return(invisible(NULL))
}


is_segments_dir <- function (dir) {
  if (length(dir) > 1 || !file.exists(dir) || !file.info(dir)$isdir)
    return(FALSE)
  
  if (split_path(dir) == "segments" && 
        all(grepl("genomes|map|\\d+|fsa\\.mfa", dir(dir))))
    return(TRUE)
  
  if (file.exists(file.path(dir, "segments")) && 
        all(grepl("genomes|map|\\d+|fsa.mfa", dir(file.path(dir, "segments")))))
    return(TRUE)
  
  else 
    return(FALSE)
}


segments_dir <- function (dir) {
  dir <- rmisc::trim(dir, paste0(.Platform$file.sep, "$"))
  dir <- if (split_path(dir) == "segments")
    dir
  else
    file.path(dir, "segments")
  
  return(dir)
}










