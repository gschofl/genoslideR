#' @importFrom Biostrings reverseComplement
#' @importFrom IRanges unlist
#' @importFrom IRanges levels
#' @importFrom IRanges runValue
#' @importFrom IRanges runLength
#' @importFrom IRanges PartitioningByWidth
#' @importFrom IRanges relist
NULL


#' Map genomic positions to alignment positions.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} or
#' \code{\linkS4class{GRangesList}} containing genomic
#' positions to map to the alignment.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param simplify if \code{TRUE}, return \code{GRanges}, otherwise
#' return a \code{GRangesList}.
#' @return A \code{\linkS4class{GRangesList}} object.
#' object.
#' @seealso \code{\link{alignment2Genome}}, \code{\link{genome2Genome}},
#' \code{\link{sliceAlignment}}, \code{\link{findAnnotation}}.
#' @export
genome2Alignment <- function(ranges, aln, simplify = FALSE) {
  
  if (missing(ranges) || !(is(ranges, "GRanges") || is(ranges, "GRangesList"))) {
    stop("Provide 'GRanges' or 'GRangesList' to map genomic positions to an alignment.",
         call.=FALSE)
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide 'annotatedAlignment' to map genomic positions to an alignment")
  }
  
  mapping_ranges <- ranges(ranges)
  genome <- unique(as.character(runValue(seqnames(ranges))))
  if (length(genome) > 1) {
    stop("Provide only one genome for mapping", call.=FALSE)
  }
  
  if (is(mapping_ranges, "IRangesList")) {
    mapping_ranges <- unlist(mapping_ranges)
  }
  
  gmap <- gMap(aln, compact = TRUE)[[genome]]
  amap <- aMap(aln, compact = TRUE)[[genome]]
  gaps <- genoslideR::gaps(aln)[[genome]]
  alnrange <- map2aln(ranges=mapping_ranges, gmap, amap, gaps)
  if (simplify)
    alnrange <- unlist(alnrange)
  alnrange 
}


#' Map alignment positions to genomic positions.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} or
#' \code{\linkS4class{GRangesList}} containing alignment
#' positions to map to the \code{targetGenomes}.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param targetGenomes a character vector of genome names.
#' @return A \code{\linkS4class{GRangesList}} object.
#' @seealso \code{\link{genome2Alignment}}, \code{\link{genome2Genome}},
#' \code{\link{sliceAlignment}}, \code{\link{findAnnotation}}.
#' @export
alignment2Genome <- function (ranges, aln, targetGenomes = NULL) {
  
  if (missing(ranges)) {
    stop("Provide 'IRanges' or 'GRanges' to map alignment positions to the genomes." )
  }
  
  if (missing(aln) || !is_mapped_alignment(alignment(aln))) {
    stop("Provide an 'annotatedAlignment' to map alignment positions to the genomes.")
  }
  
  if (is(ranges, "GRanges") || is(ranges, "GRangesList")) {
    ranges <- ranges(ranges)  
  }
  
  if (!is(ranges, "IRangesList")) {
    ranges <- unname(split(ranges, seq_along(ranges)))
  }
  
  genome <- targetGenomes %|null|% seqlevels(aln)
  if (!all(genome %in% seqlevels(aln))) {
    stop("Invalid target genomes provided")
  }
  
  gaps <- genoslideR::gaps(aln)[genome]
  gmap <- gMap(aln, compact = TRUE)[genome]
  amap <- aMap(aln, compact = TRUE)[genome]

  maprange <- aln2map(ranges, gmap, amap, gaps, normalize=FALSE)

  maprange
}


#' Map genomic positions in one genome to genomic positions in other genomes.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} or
#' \code{\linkS4class{GRangesList}} containing genomic
#' positions to map to the \code{targetGenomes}.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param targetGenomes a character vector of genome names.
#' @return A \code{\linkS4class{GRangesList}} object.
#' @seealso  \code{\link{alignment2Genome}}, \code{\link{genome2Alignment}},
#' \code{\link{sliceAlignment}}, \code{\link{findAnnotation}}.
#' @export
genome2Genome <- function (ranges, aln, targetGenomes = NULL) {
  
  if (missing(ranges) || !(is(ranges, "GRanges") || is(ranges, "GRangesList"))) {
    stop("Provide 'GRanges' or 'GRangesList' to map genomic positions.",
         call.=FALSE)
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide 'annotatedAlignment' to map genomic positions")
  }
  
  mapping_ranges <- ranges(ranges)
  sourceGenome <- unique(as.character(runValue(seqnames(ranges))))
  if (length(sourceGenome) > 1) {
    stop("Provide only one genome for mapping", call.=FALSE)
  }
  
  if (is(mapping_ranges, "GRangesList")) {
    mapping_ranges <- unlist(mapping_ranges)
  }
  
  gmap <- gMap(aln, compact = TRUE)
  amap <- aMap(aln, compact = TRUE)
  gaps <- genoslideR::gaps(aln)
  alnrange <- map2aln(ranges=mapping_ranges, gmap[[sourceGenome]],
                      amap[[sourceGenome]], gaps[[sourceGenome]])
  
  targetGenomes <- targetGenomes %|null|% seqlevels(aln)
  if (!all(targetGenomes %in% seqlevels(aln))) {
    stop("Invalid target genomes provided")
  }
  
  aln2map(ranges=ranges(alnrange), gmap=gmap[targetGenomes],
          amap=amap[targetGenomes], gaps=gaps[targetGenomes])
}


#' Slice genomic ranges from an alignment.
#' 
#' The ranges passed to \code{sliceAlignment} are converted
#' from genomic ranges to alignment ranges. The alignment ranges
#' are returned as metadata with the resulting
#' \code{\linkS4Class{DNAStringSetList}} object.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} or
#' \code{\linkS4class{GRangesList}} containing genomic positions.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param targetGenomes a character vector of genome names.
#' @return A \code{\linkS4class{DNAStringSetList}} object.
#' @seealso  \code{\link{alignment2Genome}}, \code{\link{genome2Alignment}},
#' \code{\link{genome2Genome}}, \code{\link{findAnnotation}}.
#' @export
sliceAlignment <- function (ranges, aln, targetGenomes = NULL) {
  
  if (missing(ranges) || !(is(ranges, "GRanges") || is(ranges, "GRangesList"))) {
    stop("Provide 'GRanges' or 'GRangesList' to slice an alignment",
         call.=FALSE)
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide an 'annotatedAlignment' to map genomic positions to the alignment")
  }
  
  genome <- unique(as.character(runValue(seqnames(ranges))))
  if (length(genome) > 1) {
    stop("Provide only one genome for mapping", call.=FALSE)
  }
  
  mapping_ranges <- ranges(ranges)
  if (is(mapping_ranges, "IRangesList")) {
    mapping_ranges <- unlist(mapping_ranges)
  }
  
  targetGenomes <- targetGenomes %|null|% seqlevels(aln)
  gmap <- gMap(aln, compact = TRUE)[[genome]]
  amap <- aMap(aln, compact = TRUE)[[genome]]
  gaps <- genoslideR::gaps(aln)[[genome]]
  aln <- alignment(aln)[targetGenomes]
  alnranges <- map2aln(mapping_ranges, gmap, amap, gaps)
  sliceAlnranges(alnranges, aln)
}


#' Find annotation associated with alignment positions.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} or
#' \code{\linkS4class{GRangesList}} containing alignment positions.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param targetGenomes a character vector of genome names.
#' @return A list of \code{\linkS4class{AnnotationList}} objects.
#' @seealso  \code{\link{alignment2Genome}}, \code{\link{genome2Alignment}},
#' \code{\link{genome2Genome}}, \code{\link{sliceAlignment}}.
#' @export
findAnnotation <- function (ranges, aln, targetGenomes = NULL) {
  
  if (missing(ranges)) {
    stop("Provide 'IRanges(List)' or 'GRanges(List)'" )
  }
  
  if (missing(aln) || !is_mapped_alignment(alignment(aln))) {
    stop("Provide 'annotatedAlignment' to map alignment positions to the genomes.")
  }
  
  if (is(ranges, "list") && names(ranges) == "alignment_positions") {
    ranges <- ranges$alignment_positions
  }
  
  if (is(ranges, "GRanges") || is(ranges, "GRangesList")) {
    ranges <- ranges(ranges)  
  }
  
  if (!is(ranges, "IRangesList")) {
    ranges <- unname(split(ranges, seq_along(ranges)))
  }
  
  targetGenomes <- targetGenomes %|null|% seqlevels(aln)
  if (!all(targetGenomes %in% seqlevels(aln))) {
    stop("Invalid target genomes provided")
  }
  
  gmap <- gMap(aln, compact = TRUE)[targetGenomes]
  amap <- aMap(aln, compact = TRUE)[targetGenomes]
  gaps <- genoslideR::gaps(aln)[targetGenomes]
  anno <- annotation(aln)[targetGenomes]
  
  nm <- names(ranges)[unlist(lapply(ranges, length) != 0)]
  if (is.null(nm)) {
    nm <- rep(seq_along(ranges), width(ranges@partitioning)) 
  } 

  mr <- aln2map(ranges, gmap, amap, gaps)
  if (length(nm) != unlist(unique(lapply(mr, length)))) {
    stop("Unequal number of mapped ranges")
  }
  
  hl <- findOverlaps(ranges(mr), ranges(anno))
  splitter <- factor(unlist(lapply(hl, queryHits), use.names=FALSE))
  hits <- split(anno@unlistData[subjectHits(hl), ], splitter)
  
  res <- lapply(hits, .newAnnotationList)
  names(res) <- nm
  res
}


.newAnnotationList <- function (range) {
  if (!is(range, "Granges") && !all(names(mcols(range)) %in% 
       c("type","gene","synonym","geneID","proteinID","product"))) {
    stop("Invalid range")
  }
  sqn <- seqnames(range)
  width <- rep(0L, nlevels(sqn))
  width[IRanges::levels(sqn) %in% IRanges::runValue(sqn)] <- IRanges::runLength(sqn)
  partitioning <- PartitioningByWidth(width, names=seqlevels(range))
  new("annotationList", relist(range, partitioning))
}


setAs("DNAStringSet", "DNAMultipleAlignment",
      function (from) {
        rmask <- as(IRanges(), "NormalIRanges")
        cmask <- as(IRanges(), "NormalIRanges")
        DNAMultipleAlignment(from, rowmask=rmask, colmask=cmask)
      })

