#' @importFrom IRanges IRanges
#' @importFrom IRanges mapply
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @importFrom IRanges shift
#' @importFrom GenomicRanges strand
NULL

update_alignment_position <- function (cr, ghr, ahr) {
  rev <- strand(ghr) == "-"
  if (any(rev)) {
    cr <- reverse_position(cr, hr=ahr, rev)
  }
  shift <- end(ghr) - end(ahr)
  shift(cr, shift)
}


update_genomic_position <- function(cr, ghr, ahr) {
  rev <- strand(ghr) == "-"
  if (any(rev)) {
    cr <- reverse_position(cr, ghr, rev)
  }
  shift <- start(ahr) - start(ghr)
  shift(cr, shift)
}


reverse_position <- function(cr, hr, rev) {
  rev <- which(as.logical(rev))
  cr_start <- start(cr)
  cr_end <- end(cr)
  shift <- start(hr)[rev] + end(hr)[rev] - cr_start[rev] - cr_end[rev]
  cr_start[rev] <- cr_start[rev] + shift
  cr_end[rev] <- cr_end[rev] +  shift
  rev_range <- IRanges(cr_start, cr_end, names=names(cr))
  ranges(cr) <- rev_range
  cr
}



