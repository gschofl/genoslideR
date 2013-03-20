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


aln2map <- function (ranges, gmap, amap, gaps) {
  mapping_ranges <- mapply(ungap_alignment_position, cr=ranges, MoreArgs=list(gaps = gaps))
  mapped_ranges <- vector("list", length(mapping_ranges))
  for (i in seq_along(mapping_ranges)) {
    mapped_ranges[[i]] <- 
      unlist(GRangesList(mapply(.aln2map, ranges=mapping_ranges[[i]],
                                gmap=gmap, amap=amap)))
  }
  GRangesList(mapped_ranges)
}

#.aln2map(ranges=mapping_ranges[[i]][[1]], gmap=gmap[[1]], amap=amap[[1]])


ungap_alignment_position <- function (cr, gaps) {
  gaps <- ranges(gaps)
  start <- start(cr)
  end <- end(cr)
  nm <- names(cr) %|null|% rep("", length(start))
  IRangesList(lapply(gaps, make_ungapped_ranges,
                     start = start, end = end, nm = nm))
}


.aln2map <- function(ranges, gmap, amap) {
  ovl <- findOverlaps(ranges, ranges(amap), type="any")
  ovl_ranges <- ranges(ovl, ranges, ranges(amap))
  names(ovl_ranges) <- rep(unique(names(ranges)), length(ovl_ranges))
  subject_hits <- subjectHits(ovl)

  if (length(ranges) == 1) {
    ahr <- amap[subject_hits, ]
    ghr <- gmap[subject_hits, ]
    cr <- ghr
    ranges(cr) <- ovl_ranges
    cuts <- update_alignment_position(cr, ghr, ahr)
    return( cuts )
  } else {
    query_hits <- queryHits(ovl)
    cuts <- vector("list", length(ranges))
    for (i in unique(query_hits)) {
      query <- which(query_hits == i)
      subject <- subject_hits[query]
      o <- order(start(amap)[subject])
      query_order <- query[o]
      subject_order <- subject[o]
      ahr <- amap[subject_order, ]
      ghr <- gmap[subject_order, ]
      cr <- ghr
      ranges(cr) <- ovl_ranges[query_order, ]
      cuts[i] <- update_alignment_position(cr, ghr, ahr)
    } 
    return( unlist(GRangesList(cuts)) )
  }
}

