#' @importFrom IRanges compact
#' @importFrom GenomicRanges seqlevels
#' @importFrom Biostrings xscat
#' @importFrom Biostrings subseq
#' @importFrom Biostrings reverse
#' @importFrom rmisc %|null|%
NULL


#' [INTERNAL] Slice genomic ranges mapped to an alignment
#'
#' @param ranges A \code{\linkS4class{Granges}} object.
#' @param aln A \code{\linkS4class{AnnotatedAlignment}} object.
#' 
#' @keywords internal
sliceAlnranges <- function (alnranges, aln) {
  metadata(aln) <- list()
  ar <- alnranges[unlist(lapply(alnranges, length), use.names=FALSE) != 0]
  listData <- lapply(ar, .slicer, aln)
  ans <- Biostrings:::XStringSetList("DNA", listData)
  metadata(ans) <- list(alignment_positions = alnranges)
  ans
}


.slicer <- function(ar, aln) {
  dss <- mapply(subseq2, start = start(ar), end = end(ar),
                strand = as.character(strand(ar)), MoreArgs=list(x = aln))
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
