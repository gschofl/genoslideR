#' @importFrom Biostrings reverseComplement
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges as.list
NULL

#' Map genomic positions to alignment positions.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} object containing genomic
#' positions to be mapped to the alignment.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param simplify if \code{TRUE}, return a \code{GRanges} instance.
#' @return A \code{\linkS4class{GRanges}} or \code{\linkS4class{GRangesList}}
#' object.
#' @export
map2Alignment <- function(ranges, aln, simplify = TRUE) {
  
  if (missing(ranges) || !(is(ranges, "GRanges") || is(ranges, "GRangesList"))) {
    stop("Provide a 'GRanges' or 'GRangesList' object to map genomic positions to the alignment.",
         call.=FALSE)
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide an 'annotatedAlignment' object to map genomic positions to the alignment")
  }
  
  names <- names(ranges)
  alnrange <- map2aln(ranges, aln)
  names(alnrange) <- names
  if (simplify)
    alnrange <- unlist(alnrange)
  alnrange 
}


#' Map alignment positions to genomic positions.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} object containing alignment
#' positions to be mapped to the genomes.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param targetGenomes
#' @return A \code{\linkS4class{GRangesList}} object.
#' @export
map2Genome <- function(ranges, aln, targetGenomes = NULL) {
  
  if (missing(ranges)) {
    stop("Provide 'IRanges' or 'GRanges' to map alignment positions to the genomes." )
  }
  
  if (missing(aln) || !is_mapped_alignment(alignment(aln))) {
    stop("Provide an 'annotatedAlignment' to map alignment positions to the genomes.")
  }
  
  names <- names(ranges)
  maprange <- aln2map(ranges, aln, genome=targetGenomes)
  names(maprange) <- names
  maprange
}


#' Extract genomic ranges from an alignment.
#' 
#' @param ranges a \code{\linkS4class{GRanges}} object containing genomic
#' positions to be mapped to the alignment.
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param addAnnotations
#' @param targetGenomes
#' @return A \code{\linkS4class{GRangesList}} object.
#' @export
sliceAlignment <- function(ranges, aln, addAnnotations = TRUE, targetGenomes = NULL) {
  
  if (missing(ranges) || !(is(ranges, "GRanges") || is(ranges, "GRangesList"))) {
    stop("Provide a 'GRanges' or 'GRangesList' to map genomic positions to the alignment",
         call.=FALSE)
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide an 'annotatedAlignment' to map genomic positions to the alignment")
  }
  
  aln_ranges <- map2aln(ranges, aln)
  aln_slices <- slice_aln(ranges=aln_ranges, aln, targetGenomes = targetGenomes)  
  for (i in seq_along(aln_slices)) {
    rev <- which(strand(metadata(aln_slices[[i]])[["alignment_position"]]) == "-")
    aln_slices[rev] <- lapply(aln_slices[rev], function (slcd) {
      reverseComplement(as(slcd, "DNAStringSet"))
    })
  }

  if (addAnnotations) {
    aln_slices <- add_annotation(aln_slices, aln, type="any",
                                 targetGenomes = targetGenomes)
  }
    
  return(aln_slices)
}

# slcd=aln_slices[[1]]
# aln=aa

add_annotation <- function (aln_slices, aln, targetGenome, type = "any") {
  gaps <- genoslideR::gaps(aln)[targetGenomes]
  gmap <- gMap(aln)[targetGenomes]
  amap <- aMap(aln)[targetGenomes]
  
  aa <- lapply(aln_slices, function (slcd) {
    apos <- ranges(metadata(slcd)[["alignment_position"]])
    mapping_ranges <- ungap_alignment_position(apos, gaps)
    mapped_ranges <- vector("list", length(mapping_ranges))
    for (i in seq_along(mapping_ranges)) {
      mapped_ranges[[i]] <- GRangesList(mapply(.aln2map,
                                               
                                               ranges = mapping_ranges[i],
                                               gmap = gmap, amap = amap))
    }
    GRangesList(mapped_ranges)
    
    
    
    gpos <- aln2map(ranges=apos, aln, genome=names(slcd))
    gpos <- split(gpos, seqnames(gpos))
    hitlist <- findOverlaps(ranges(gpos), ranges(aln), type=type)
    shits <- lapply(hitlist, unique%.%subjectHits)
    an <- new("annotationList", GRangesList(Map(function (anno, hit) {
      anno[hit,] }, anno = as.list(annotation(aln)), hit = shits)))
    new("annotatedAlignment", annotation = an, alignment = slcd)
  })
  aa
}
