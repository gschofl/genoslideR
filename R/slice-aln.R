#' @importFrom XVector compact
#' @importFrom GenomicRanges seqlevels
#' @importFrom Biostrings xscat reverse
NULL


#' [INTERNAL] Slice GRangesList from an alignment
#'
#' @param grl A \code{\linkS4class{GrangesList}} object.
#' @param dss A \code{\linkS4class{DNAStringSet}} object.
#' 
#' @keywords internal
sliceGRL <- function (grl, dss, compact = FALSE) {
  dss@metadata <- list()
  grl <- grl[unlist(lapply(grl, length), use.names=FALSE) != 0]
  listData <- lapply(grl, sliceGR, dss = dss, compact = compact)
  ans <- Biostrings:::XStringSetList("DNA", listData)
  metadata(ans) <- list(alignment_positions = grl)
  ans
}


#' [INTERNAL] Slice GRanges from an alignment
#'
#' @param dss A \code{\linkS4class{DNAStringSet}} object.
#' @param gr A \code{\linkS4class{Granges}} object.
#'
#' @keywords internal
sliceGR <- function(gr, dss, compact = FALSE) {
  start <- start(gr)
  width <- width(gr)
  reverse <- as.integer(strand(gr)) == 2L
  .slice(dss, start, width, reverse, compact=compact)
}


.slice <- function (x, start, width, reverse = FALSE, complement = FALSE, compact = FALSE) {
  if (length(start) == 1)
    return( narrowAlignment(x, start, width, reverse, complement, compact) )
  
  xl <- mapply(narrowAlignment, start = start, width = width, reverse = reverse,
               complement = complement, compact = compact, MoreArgs=list(x = x))
  xl <- setNames(.Call2("XStringSet_xscat", xl, PACKAGE = "Biostrings"), 
           nm = names(xl[[1L]]))
  xl
}


narrowAlignment <- function (x, start, width, reverse=FALSE, complement=FALSE,
                             compact=FALSE) {
  
  range <- narrowAlnrange(x@ranges, start, width)
  slot(x, "ranges", check = FALSE) <- range
  
  if (compact || complement || reverse) {
    
    if (complement)
      lkup <- Biostrings:::getDNAComplementLookup()
    else
      lkup <- NULL
    
    x <- xvcopy(x, lkup = lkup, reverse = reverse)
  }
  
  x
}


narrowAlnrange <- function (x, start, width) {
  SEW <- .Call2("solve_user_SEW", refwidths = x@width, start = start,
                end = NA_integer_, width = width,
                translate.negative.coord = TRUE, 
                allow.nonnarrowing = FALSE, PACKAGE = "IRanges")
  slot(x, "start", check=FALSE) <- x@start + SEW@start - 1L
  slot(x, "width", check=FALSE) <- SEW@width
  x
}




