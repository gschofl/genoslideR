#' @importFrom IRanges IRanges
#' @importFrom IRanges IRangesList
#' @importFrom IRanges subjectHits
#' @importFrom IRanges findOverlaps
#' @importFrom GenomicRanges unlist
NULL

ungap_range <- function (cr, gmap, amap, gaps) {
  # ungapped genomic ranges
  ugr <- IRangesList(lapply(gaps, make_ungapped_ranges,
                            start = start(cr), end = end(cr)))
  # remap ranges
  GRangesList(mapply(remap_ranges, ugr = ugr, gmap = gmap,
                     amap = amap, genome = as.list(names(ugr))))
}


remap_ranges <- function (ugr, gmap, amap, genome) {
  hit <- subjectHits(findOverlaps(ugr, amap, type="any"))
  hrg <- gmap[hit, ] # hit range genome
  hra <- amap[hit, ] # hit range alignment
  map <- IRanges(end(hrg) - (end(hra) - start(ugr)), end(hrg) - (end(hra) - end(ugr)))
  if (length(map) > 0) {
    GRanges(genome, map)
  } else {
    GRanges()
  }
}


gap_range <- function (start, end, gaprange, maprange) {
  
  gap_start <- GenomicRanges::start(gaprange)
  gap_width <- GenomicRanges::width(gaprange)
  map_start <- min(GenomicRanges::start(maprange))
  map_end <- max(GenomicRanges::end(maprange)) - map_start + 1L
  aln_len <- unique(GenomicRanges::seqlengths(gaprange))
  if (is.na(aln_len)) {
    stop("No seqlengths in gapranges", call.=TRUE)
  }
  
  start <- maximum(minimum(start - map_start + 1L, 1L), map_end + 1L)
  end <- maximum(minimum(end - map_start + 1L, 0L), map_end)
  end <- ifelse(start > end, start, end)
  make_gapped_ranges(start, end, gap_start, gap_width, aln_len)
}


minimum <- function(a, b) {
  ifelse(a < b, b, a)
}


maximum <- function(a, b) {
  ifelse(a > b, b, a)
}

