#' @importFrom Biostrings reverseComplement
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges as.list
NULL

#' Slice an Annotated Alignment
#' 
#' @param aln an \code{\linkS4class{annotatedAlignment}} instance.
#' @param slices a \code{\linkS4class{GRanges}} instance.
#' @param addAnnotation if \code{TRUE} an \code{\linkS4class{annotatedAlignment}}
#' is returned, otherwise a \code{\linkS4class{DNAStringSet}} object.
#' @export
sliceAnnotatedAlignment <- function(aln, slices, addAnnotation = TRUE) {
  
  if (!is(slices, "GRanges")) {
    stop("Provide a 'GRanges' object to slice the alignment")
  }
  
  if (!is_mapped_alignment(alignment(aln))) {
    stop("Provide a mapped 'annotatedAlignment' object to slice the alignment")
  }
  
  sliced <- get_slices(aln, slices)  
  rev <- which(strand(slices) == "-")
  sliced[rev] <- lapply(sliced[rev], function (slcd) {
    reverseComplement(as(slcd, "DNAStringSet"))
  })
  
  if (addAnnotation)
    sliced <- add_annotation(aln, sliced, type="any")
  
  sliced
}


add_annotation <- function (aln, sliced, type = "any") {
  aa <- lapply(sliced, function (slcd) {
    gpos <- metadata(slcd)[["genomic_position"]]
    hitlist <- findOverlaps(ranges(gpos), ranges(aln), type=type)
    shits <- lapply(hitlist, unique%.%subjectHits)
    an <- new("annotationList", GRangesList(Map(function (anno, hit) {
      anno[hit,] }, anno = as.list(annotation(aln)), hit = shits)))
    new("annotatedAlignment", annotation = an, alignment = slcd)
  })
  aa
}