#' @importFrom IRanges compact
#' @importFrom GenomicRanges seqlevels
#' @importFrom Biostrings xscat
#' @importFrom Biostrings subseq
#' @importFrom Biostrings reverse
#' @importFrom Biostrings XStringSetList
#' @importFrom rmisc %|null|%
NULL


#' INTERNAL: Slice genomic ranges mapped to an alignment
#'
#' @param ranges A \code{\linkS4class{Granges}} object.
#' @param aln A \code{\linkS4class{AnnotatedAlignment}} object.
#' 
#' @keywords internal
sliceAlnRanges <- function (ranges, aln) {
  metadata(aln) <- list()
  listData <- mapply(.slicer, r=ranges, MoreArgs=list(aln = aln))
  ans <- Biostrings:::XStringSetList("DNA", listData)
  metadata(ans) <- list(alignment_positions = ranges)
  ans
}


.slicer <- function(r, aln) {
  dss <- mapply(subseq2, start = start(r), end = end(r),
                strand = as.character(strand(r)), MoreArgs=list(x = aln))
  if (length(dss) == 1) {
    dss <- IRanges::compact(dss[[1]])
  } else {
    dss <- IRanges::compact(setNames(do.call(xscat, dss), nm=names(dss[[1]])))
  }
  dss
}


subseq2 <- function(x, start, end, strand) {
  if (as.character(strand) == "-") {
    reverse(subseq(x, start, end))
  } else {
    subseq(x, start, end)
  }
} 
