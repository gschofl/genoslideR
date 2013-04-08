#' @importFrom IRanges ranges
#' @importFrom IRanges split
#' @importFrom IRanges mapply
#' @importFrom IRanges findOverlaps 
#' @importFrom IRanges queryHits 
#' @importFrom IRanges subjectHits 
#' @importFrom IRanges IRangesList
#' @importFrom GenomicRanges seqlevels
#' @importFrom GenomicRanges GRangesList
#' @importFrom rmisc %|null|%
NULL


#' [INTERNAL] Map a alignment ranges to genomes
#'
#' @param ranges An IRangesList instance
#' @param gmap A GrangesList instance.
#' @param amap A GrangesList instance.
#' @param gaps A GrangesList instance.
#' @param normalize Normalize ovelapping query ranges.
#' @return A GrangesList instance. 
#' 
#' @keywords internal
aln2map <- function (ranges, gmap, amap, gaps, normalize = FALSE) {
  ar <- unlist(ranges)
  if (is.null(names(ar))) {
    names(ar) <- rep(seq_along(ranges), width(ranges@partitioning)) 
  } 
  mapping_ranges <- IRangesList(ungap_alignment_position(ar, gaps))
  if (normalize) {
    mapping_ranges <- lapply(mapping_ranges, as, "NormalIRanges")
    mapping_ranges <- IRangesList(lapply(mapping_ranges, as, "IRanges"))
  }
#   i <- 8
#   x <- .aln2map(mr=mapping_ranges[[i]], gm=gmap[[i]], am=amap[[i]])
#   x
  GRangesList(mapply(.aln2map, mr = mapping_ranges, gm = gmap, am = amap))
}




ungap_alignment_position <- function (cr, gaps) {
  gaprange <- ranges(gaps)
  start <- start(cr)
  end <- end(cr)
  nm <- names(cr) %|null|% rep("", length(cr))
  lapply(gaprange, make_ungapped_ranges, start = start, end = end, nm = nm)
}


## mr : mapping ranges
## gr : genomic map
## am : alignment map
.aln2map <- function(mr, gm, am) {
  if (length(mr) == 0 || all(width(mr) == 0))  
    return(GRanges(runValue(seqnames(gm)), IRanges(1, width=0, names=""), "*"))

  ar <- ranges(am)
  ovl <- findOverlaps(mr, ar, type="any")
  query_hits <- queryHits(ovl)
  subject_hits <- subjectHits(ovl)
  ovl_ranges <- ranges(ovl, mr, ar)
  names(ovl_ranges) <- names(mr)[query_hits]
  genomic_hits <- gm[subject_hits, ]
  gh_ranges <- ranges(genomic_hits)
  gh_strand <- as.integer(strand(genomic_hits))
  ah_ranges <- ar[subject_hits, ]
  update_alignment_position_cpp(ovl_ranges, gh_ranges, ah_ranges, gh_strand)
  ranges(genomic_hits) <- ovl_ranges
  genomic_hits
}

