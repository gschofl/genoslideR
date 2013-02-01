#' Align orthologous genome segments with \code{fsa} 
#' 
#' @param seg_dir Path to mercator segments.
#' @param force Realign an existing alignment.
#' @param opts a named list of options for \code{\link{fsa}}. The default
#' options are optimizes for long sequences.
#' @param ... Named values interpreted as options for fsa.
#' 
#' @export
alignSegments <- function (seg_dir, force = FALSE,
                           opts = list(anchored = TRUE,
                                       translated = TRUE,
                                       exonerate = TRUE,
                                       softmasked = TRUE),
                           ...) {
  
  if (is_segments_dir(seg_dir)) {
    seg_dir <- segments_dir(seg_dir)
  } else {
    stop("No valid path to mercator segments")
  }
  
  opts <- merge_list(opts, list(...))
  seg_dir <- align_mercator_segments(seg_dir, force, opts)

  message("Pass the path to the aligned mercator segments to the 'annotatedAlignment()' constructor to create an 'annotatedAlignment' object")
  return(seg_dir)
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
  dir <- if (split_path(dir) == "segments")
    dir
  else
    file.path(dir, "segments")
  
  return(dir)
}










