#' @importFrom IRanges which
#' @importFrom IRanges order
#' @importFrom IRanges compact
#' @importFrom IRanges subseq
#' @importFrom IRanges mendoapply
#' @importFrom IRanges ranges
#' @importFrom IRanges split
#' @importFrom IRanges runValue
#' @importFrom GenomicRanges seqnames
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings xscat
NULL

#' INTERNAL: Map a genomic range to the alignment
#'
#' @param range A \code{\linkS4class{Granges}} object.
#' @param aln A \code{\linkS4class{AnnotatedAlignment}} object.
#' 
#' @keywords internal
map2aln <- function (ranges, aln) {

  mapping_ranges <- ranges(ranges)
  genome <- unique(as.character(runValue(seqnames(ranges))))
  if (length(genome) > 1) {
    stop("Provide only one genome for mapping", call.=FALSE)
  }
  
  if (is(mapping_ranges, "IRangesList")) {
    mapping_ranges <- unlist(mapping_ranges)
  }
  
  gmap <- gMap(aln)[[genome]]
  amap <- aMap(aln)[[genome]]
  gaps <- genoslideR::gaps(aln)[[genome]]
  map_min <- min(start(gmap))
  map_max <- max(end(gmap)) - map_min + 1L

  ovl <- findOverlaps(mapping_ranges, ranges(gmap), type="any")
  ovl_ranges <- ranges(ovl, mapping_ranges, ranges(gmap))
  subject_hits <- subjectHits(ovl)
  query_hits <- queryHits(ovl)
  cuts <- vector("list", length(mapping_ranges))
  for (i in unique(query_hits)) {
    query <- which(query_hits == i)
    subject <- subject_hits[query]
    o <- order(start(gmap)[subject])
    query_order <- query[o]
    subject_order <- subject[o]
    ahr <- amap[subject_order, ]
    ghr <- gmap[subject_order, ]
    cr <- ghr
    ranges(cr) <- ovl_ranges[query_order, ]
    cuts[i] <- .map2aln(cr, ghr, ahr, gaps, map_min, map_max)
  }

  cuts[-unique(query_hits)] <- GRanges()
  cuts <- GRangesList(cuts)
  cuts
}


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
.map2aln <- function (cr, ghr, ahr, gaps, map_min, map_max) {
  # update cut ranges to the degapped alignment positions (dapos)
  dapos <- update_genomic_position(cr, ghr, ahr)
  # update degapped alignment postions to the gapped alignment positions (gapos)
  gapos <- gap_alignment_position(dapos, map_min, map_max, gaps)
  ranges(dapos) <- gapos
  seqlengths(dapos) <- NA
  dapos
}


gap_alignment_position <- function (dapos, map_min, map_max, gaps) {
  start <- maximum(minimum(start(dapos) - map_min + 1L, 1L), map_max + 1L)
  end <- maximum(minimum(end(dapos) - map_min + 1L, 0L), map_max)
  end <- ifelse(start > end, start, end)
  make_gapped_ranges(start, end, gaps)
}


minimum <- function(a, b) {
  ifelse(a < b, b, a)
}


maximum <- function(a, b) {
  ifelse(a > b, b, a)
}

