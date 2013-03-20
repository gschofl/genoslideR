#' @importFrom IRanges compact
#' @importFrom GenomicRanges seqlevels
#' @importFrom Biostrings xscat
#' @importFrom Biostrings subseq
#' @importFrom Biostrings reverse
#' @importFrom rmisc %|null|%
NULL


#' INTERNAL: Slice genomic ranges mapped to an alignment
#'
#' @param ranges A \code{\linkS4class{Granges}} object.
#' @param aln A \code{\linkS4class{AnnotatedAlignment}} object.
#' 
#' @keywords internal
slice_aln <- function (ranges, aln, targetGenomes = NULL) {
  
  targetGenomes <- targetGenomes %|null|% seqlevels(aa)
  aln <- alignment(aa)[targetGenomes]
  
  mapply(.slice_aln, range=ranges, MoreArgs=list(aln = aln))
}


.slice_aln <- function(range, aln) {
  dss <- mapply(subseq2, start = start(range), end = end(range),
                strand = as.character(strand(range)), MoreArgs=list(x = aln))
  if (length(cut) == 1) {
    dss <- IRanges::compact(dss[[1]])
  } else {
    dss <- IRanges::compact(setNames(do.call(xscat, dss), nm=names(dss[[1]])))
  }
  metadata(dss) <- list(alignment_position = range)
  dss
}


subseq2 <- function(x, start, end, strand) {
  if (as.character(strand) == "-") {
    reverse(subseq(x, start, end))
  } else {
    subseq(x, start, end)
  }
} 
