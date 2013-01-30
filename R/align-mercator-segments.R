#' Align orthologous genome segments
#' 
#' @param seg_dir Path to mercator segments.
#' @param aligner Which alignment program to use; 'fsa' or 'mavid'.
#' @param ... Other params.
#' 
#' @export
alignMercatorSegments <- function (seg_dir, aligner = "fsa", force = FALSE, mask = TRUE) {
  
  aligner <- match.arg(aligner, c("fsa"))
  
  if (is_segments_dir(seg_dir)) {
    seg_dir <- segments_dir(seg_dir)
  } else {
    stop("No valid path to mercator segments")
  }
  
  seg_dir <- align_mercator_segments(seg_dir, aligner, force, mask)

  message("Pass the path to the aligned mercator segments to the 'annotatedAlignment()' constructor to create an 'annotatedAlignment' object")
  return(seg_dir)
}


is_segments_dir <- function (dir) {
  if (length(dir) > 1 || !file.exists(dir) || !file.info(dir)$isdir)
    return(FALSE)
  
  if (split_path(dir) == "segments" && all(grepl("genomes|map|\\d+", dir(dir))))
    return(TRUE)
  
  if (file.exists(file.path(dir, "segments")) && all(grepl("genomes|map|\\d+", dir(file.path(dir, "segments")))))
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










