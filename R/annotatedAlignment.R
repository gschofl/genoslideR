#' @include annotationList.R
#' @importClassesFrom Biostrings XStringSet
#' @importFrom BiocGenerics annotation
#' @importFrom biofiles getAccession
#' @importFrom biofiles getDefinition
NULL

#' Construct an \sQuote{annotatedAlignment}.
#' 
#' \dQuote{annotatedAlignment} is an S4 class that provides a container for
#' multiple alignments stored as an \dQuote{\linkS4class{XStringSet}} object
#' and genomic annotations stored as an \dQuote{\linkS4class{annotationList}}
#' object.
#' @export
setClass("annotatedAlignment",
         representation(annotation="annotationList",
                        alignment="XStringSet"))


#' @param alnpath Path to a directory containing mercator segments,
#' a single file in mfa (multi fasta) or maf (multiple alignment file)
#' format containing the genome aligment.
#' @inheritParams importAnnotation
#' @seealso \code{\link{mercator}}, \code{\link{alignSegments}}.
#' @export
annotatedAlignment <- function (alnpath, annopath, type, ...) {
  
  if (missing(alnpath)) {
    stop("No alignment provided")
  }

  if (missing(annopath)) {
    stop("No annotation provided")
  }
  
  if (missing(type)) {
    stop("No annotation type provided")
  }
  
  if(!is_mapped_alignment(alnpath)) {
    alignment <- importAlignment(alnpath)
  }
  
  genomes <- names(alignment)
  alignment <- alignment[order(genomes)]
  genome_annotations <- strip_ext(basename(annopath))
  
  idx <- genomes %in% genome_annotations
  if (all(idx == FALSE)) {
    stop("The annotation file names don't match with the genome names in the alignment")
  }
  
  annopath <- annopath[idx][order(genome_annotations[idx])]
  annotation <- annotationList(annopath, type, ...)
  
  new("annotatedAlignment", annotation = annotation, alignment = alignment)
}


setValidity("annotatedAlignment", function (object) {
  
  n <- names(object)
  if (!all(n$annotation %in% n$alignment)) {
    return("Name mismatch between annotations and alignments")
  }

  met <- metadata(alignment(object))
  if (!all(idx <- c("gMap", "aMap", "gaps") %in% names(met))) {
    return(sprintf("%s missing from alignment metadate", c("map", "gaps")[idx]))
  }

  if (!is(met$gaps, "GRangesList")) {
    return("Gaps must be a 'GRangesList'")
  }
  
  if (!is(met$gMap, "GRangesList")) {
    return("Genomic map must be a 'GRangesList'")
  }
  
  if (!is(met$aMap, "GRangesList")) {
    return("Alignment map must be a 'GRangesList'")
  }
  
  if (!all(seqlevels(met$gaps) %in% n$alignment)) {
    return("Name mismatch between alignments and gap metadata")
  }
  
  if (!all(seqlevels(met$gMap) %in% n$alignment)) {
    return("Name mismatch between alignments and genomic map metadata")
  }
  
  if (!all(seqlevels(met$aMap) %in% n$alignment)) {
    return("Name mismatch between alignments and alignment map metadata")
  }
  
  if (!all(unlist(lapply(met$aMap, length)) == unlist(lapply(met$gMap, length)))) {
    return("Length mismatch between genomic map and alignment map")
  }
  
  TRUE
})


#' @export
setGeneric("annotation")
setMethod("annotation", "annotatedAlignment",
          function (object) object@annotation)


#' @export
setGeneric("alignment", signature="x",
           function (x, ...) {
             standardGeneric("alignment")
           })

setMethod("alignment", "annotatedAlignment",
          function (x) x@alignment)


setMethod("show", "annotatedAlignment",
          function (object) {
            cat(sprintf("An %s instance of length %s\n",
                        sQuote(class(object)), length(alignment(object))))
            show(alignment(object))
            cat("\n")
            show(annotation(object))  
          })


setMethod("getAccession", "annotatedAlignment",
          function (x) {
            getAccession(annotation(x))
          })


setMethod("getDefinition", "annotatedAlignment",
          function (x) {
            getDefinition(annotation(x))
          })


setMethod("seqlengths", "annotatedAlignment",
          function (x) {
            seqlengths(annotation(x))
          })


setMethod("length", "annotatedAlignment",
          function (x) {
            c(alignment = length(x@alignment),
              annotation =  length(x@annotation))
          })


setMethod("seqnames", "annotatedAlignment",
          function (x) {
            seqnames(annotation(x))
          })


setMethod("names", "annotatedAlignment",
          function (x) {
            list(alignment = names(x@alignment),
                 annotation = names(x@annotation))
          })


setMethod("start", "annotatedAlignment",
          function (x) {
            start(annotation(x))
          })


setMethod("end", "annotatedAlignment",
          function (x) {
            end(annotation(x))
          })


setMethod("ranges", "annotatedAlignment",
          function (x) {
            ranges(annotation(x))
          })


setMethod("strand", "annotatedAlignment",
          function (x, ...) {
            strand(annotation(x))
          })


setMethod("seqlevels", "annotatedAlignment",
          function (x) {
            seqlevels(annotation(x))
          })


setMethod("seqinfo", "annotatedAlignment",
          function (x) {
            seqinfo(annotation(x))
          })


setMethod("genome", "annotatedAlignment",
          function (x) {
            genome(annotation(x))
          })


setMethod("[", "annotatedAlignment",
          function (x, i, j, ..., drop = TRUE) {
            if (!missing(j)) {
              stop("invalid subsetting")
            }
            aln <- alignment(x)
            ann <- annotation(x)
            if (is.numeric(i)) {
              i <- names(aln)[i]
            }
            aln <- aln[i]
            
            gm <- which(seqlevels(gMap(x)) %in% i)
            gMap <- gMap(x)[gm]
            seqinfo(gMap, new2old = gm) <- seqinfo(gMap)[i]
              
            am <- which(seqlevels(aMap(x)) %in% i)
            aMap <- aMap(x)[am]
            seqinfo(aMap, new2old = am) <- seqinfo(aMap)[i]
            
            gps <- which(seqlevels(genoslideR::gaps(x)) %in% i)
            gaps <- genoslideR::gaps(x)[gps]
            seqinfo(gaps, new2old = gps) <- seqinfo(gaps)[i]
            
            metadata(aln) <- list(gMap =gMap, aMap = aMap, gaps = gaps)
            
            i <- names(ann)[names(ann) %in% i]
            if (all_empty(i)) {
              ann <- new("annotationList")
            } else {
              ann <- ann[i]  
            }
            new("annotatedAlignment", annotation = ann, alignment = aln)
          })


setMethod("product", "annotatedAlignment",
          function(x) {
            product(annotation(x))
          })


setMethod("locusTag", "annotatedAlignment",
          function(x) {
            locusTag(annotation(x))
          })


setMethod("geneID", "annotatedAlignment",
          function(x) {
            geneID(annotation(x))
          })


setMethod("proteinID", "annotatedAlignment",
          function(x) {
            proteinID(annotation(x))
          })


setMethod("type", "annotatedAlignment",
          function (x) {
            type(annotation(x))
          })


#' @export
setGeneric("gMap", function (x, ...) standardGeneric("gMap"))
setMethod("gMap", "annotatedAlignment",
          function (x, compact =  FALSE) {
            gmap <- metadata(alignment(x))[["gMap"]]
            if (compact) {
              gmap <- gmap[width(gmap) != 0, ]
            }
            gmap
          })


#' @export
setGeneric("aMap", function (x, ...) standardGeneric("aMap"))
setMethod("aMap", "annotatedAlignment",
          function (x, compact =  FALSE) {
            amap <- metadata(alignment(x))[["aMap"]]
            if (compact) {
              amap <- amap[width(amap) != 0, ]
            }
            amap
          })


#' @export
setGeneric("gaps", function (x, ...) standardGeneric("gaps"))
setMethod("gaps", "annotatedAlignment",
          function (x) {
            metadata(alignment(x))[["gaps"]]
          })

