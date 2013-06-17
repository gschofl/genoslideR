#' @importFrom IRanges runValue
#' @importFrom IRanges ranges
#' @importFrom IRanges which
#' @importFrom IRanges order
#' @importFrom IRanges split
#' @importFrom IRanges findOverlaps
#' @importFrom IRanges subjectHits
#' @importFrom IRanges queryHits
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges GRanges
#' @importFrom XVector subseq
NULL

#' [INTERNAL] Map a genomic range to the alignment
#'
#' @param ranges An IRanges instance
#' @param gmap A Granges instance.
#' @param amap A Granges instance.
#' @param gaps A Granges instance.
#' @return A GrangesList instance. 
#' 
#' @keywords internal
map2aln <- function (ranges, gmap, amap, gaps) {
  map_min <- map_min(gmap)
  map_max <- map_max(gmap, map_min)
  ovl <- findOverlaps(ranges, ranges(gmap), type="any")
  ovl_ranges <- ranges(ovl, ranges, ranges(gmap))
  subject_hits <- subjectHits(ovl)
  query_hits <- queryHits(ovl)
  
  order <- lapply(split(start(gmap)[subject_hits], query_hits), order)
  query_split <- split(seq_along(query_hits), query_hits)
  subject_split <- split(subject_hits, query_hits)
  query_order <- unlist(mapply(`[`, x = query_split, i = order), use.names=FALSE)
  subject_order <- unlist(mapply(`[`, x = subject_split, i = order), use.names=FALSE)
  
  genomic_hits <- gmap[subject_order, ]
  gh_ranges <- ranges(genomic_hits)
  gh_strand <- as.integer(strand(genomic_hits))
  ah_ranges <- ranges(amap)[subject_order, ]
  
  # cut_ranges
  cr <- ovl_ranges[query_order, ]
  # use side effect to directly update start positions in cr
  update_genomic_position_cpp(cr, gh_ranges, ah_ranges, gh_strand)                                 
#   start <- maximum(minimum(start(cr) - map_min + 1L, 1L), map_max + 1L)
#   end <- maximum(minimum(end(cr) - map_min + 1L, 0L), map_max)
  start <- maximum(minimum(start(cr), 1L), map_max + 1L)
  end <- maximum(minimum(end(cr), 0L), map_max)
  end <- ifelse(start > end, start, end)
  gapped_genomic_ranges <- make_gapped_ranges(start, end, gaps)
  seqlengths(genomic_hits) <- NA
  ranges(genomic_hits) <- gapped_genomic_ranges
  
  gsplit <- factor(c(query_hits, setdiff(seq_along(ranges), query_hits)))
  genomic_ranges_list <- split(genomic_hits, gsplit)
  names(genomic_ranges_list) <- names(ranges)
  genomic_ranges_list
}


# map2aln <- function (ranges, gmap, amap, gaps) {
#   map_min <- map_min(gmap)
#   map_max <- map_max(gmap, map_min)
#   ovl <- findOverlaps(ranges, ranges(gmap), type="any")
#   ovl_ranges <- ranges(ovl, ranges, ranges(gmap))
#   subject_hits <- subjectHits(ovl)
#   query_hits <- queryHits(ovl)
#   
#   cuts <- vector("list", length(ranges))
#   for (i in seq_along(ranges)) {
#     query <- which(query_hits == i)
#     subject <- subject_hits[query]
#     o <- order(start(gmap)[subject])
#     query_order <- query[o]
#     subject_order <- subject[o]
#     ahr <- amap[subject_order, ]
#     ghr <- gmap[subject_order, ]
#     cr <- ghr
#     ranges(cr) <- ovl_ranges[query_order, ]
#     cuts[i] <- .map2aln(cr, ghr, ahr, gaps, map_min, map_max)
#   }
#   
#   cuts[-unique(query_hits)] <- GRanges()
#   cuts <- GRangesList(cuts)
#   cuts
# }


#' INTERNAL: Map a genomic range to the alignment
#'
#' @param cr cut range (the range within a genomic hit range which we want
#' to slice from an alignment)
#' @param map_min Minumum position in the genomic map.
#' @param map_max Maximum position in the genomic map
#' @param ghr genomic hit range (range of an orthologous genomic segment within
#' a genome)
#' @param ahr alignment hit range (degapped range of an orthologous genomic
#' segment within a genome)
#' @param gaps gapranges
#' 
#' @keywords internal
# .map2aln <- function (cr, ghr, ahr, gaps, map_min, map_max) {
#   # update cut ranges to the degapped alignment positions (dapos)
#   dapos <- update_genomic_position(cr, ghr, ahr)
#   # update degapped alignment postions to the gapped alignment positions (gapos)
#   gapos <- gap_alignment_position(dapos, map_min, map_max, gaps)
#   ranges(dapos) <- gapos
#   seqlengths(dapos) <- NA
#   dapos
# }
# 
# 
# gap_alignment_position <- function (dapos, map_min, map_max, gaps) {
#   start <- maximum(minimum(start(dapos) - map_min + 1L, 1L), map_max + 1L)
#   end <- maximum(minimum(end(dapos) - map_min + 1L, 0L), map_max)
#   end <- ifelse(start > end, start, end)
#   make_gapped_ranges(start, end, gaps)
# }


minimum <- function(a, b) {
  ifelse(a < b, b, a)
}


maximum <- function(a, b) {
  ifelse(a > b, b, a)
}


map_min <- function(map) {
  min(start(map)[width(map) > 0])
}


map_max <- function(map, map_min = 1) {
  max(end(map)[width(map) > 0]) - map_min + 1L
}

