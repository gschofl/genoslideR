#' @importClassesFrom GenomicRanges GenomicRangesORGRangesList GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GenomicRangesList GRangesList Seqinfo
#' @importClassesFrom IRanges CompressedList Annotated List Vector
#' @importFrom GenomicRanges GRanges GRangesList Seqinfo seqinfo "seqinfo<-" seqnames
#' @importFrom GenomicRanges genome seqlengths seqlevels ranges strand mcols
#' @importFrom IRanges CharacterList IntervalTree metadata "metadata<-" elementMetadata
#' @importFrom IRanges "elementMetadata<-" Rle start end width elementLengths lapply gaps
#' @importFrom Biostrings type
#' @importFrom biofiles getAccession getDefinition qualif
NULL

#' Construct an \dQuote{annotationList}.
#' 
#' \dQuote{annotationList} is an S4 class that provides a container for
#' genomic annotations stored in \dQuote{\linkS4class{GRanges}}.
#'
#' @export
setClass("annotationList", contains="GRangesList")


#' @inheritParams importAnnotation
#' @seealso \code{\link{mercator}}, \code{\link{alignSegments}}.
#' @export
annotationList <- function(annopath, type = "gff", ...) {
  new("annotationList", importAnnotation(annopath, type, ...))
}


.show_AnnotationList <- function(x) {
  len <- length(x)
  cat(sQuote(class(x)), " of length ", len, "\n", sep = "")
  if (len == 0L)
    return(invisible(NULL))
  nm <- seqlevels(x)
  len <- elementLengths(x)
  slen <- seqlengths(x)
  cat(labeledLine("Id", els=nm, count=TRUE))
  cat(labeledLine("Annotations", els=len, count=FALSE))
  cat(labeledLine("Seqlengths", els=slen, count=FALSE))
}


setMethod("show", "annotationList",
          function(object) {
            .show_AnnotationList(object)
          })


#' @export
setGeneric("getAccession")
setMethod("getAccession", "annotationList",
          function (x) seqnames(seqinfo(x)))


#' @export
setGeneric("getDefinition")
setMethod("getDefinition", "annotationList",
          function (x) genome(seqinfo(x)))


#### GRanges and GRangesList accesssors ####

#' @export
setGeneric("seqlengths")
setMethod("seqlengths", "annotationList",
          function (x) callNextMethod())


setMethod("length", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("seqnames")
setMethod("seqnames", "annotationList",
          function (x) callNextMethod())


setMethod("names", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("start")
setMethod("start", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("end")
setMethod("end", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("ranges")
setMethod("ranges", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("strand")
setMethod("strand", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("seqinfo")
setMethod("seqinfo", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("seqinfo<-")
setMethod("seqinfo<-", "annotationList",
          function (x, new2old=NULL, force=FALSE, value) callNextMethod())


#' @export
setGeneric("seqlevels")
setMethod("seqlevels", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("genome")
setMethod("genome", "annotationList",
          function (x) callNextMethod())


#' @export
setGeneric("mcols")
setMethod("mcols", "annotationList",
          function (x) callNextMethod())


# setMethod("[", "annotationList",
#           function (x, i, j, ..., drop = TRUE) {
#             if (!missing(i)) {
#               if (is.numeric(i)) {
#                 i <- names(x)[i]
#               }
#               x <- callNextMethod(x = x, i = i, ...)
#               seql <- which(seqlevels(x) %in% i)
#               seqinfo(x, new2old = seql) <- seqinfo(x)[i]
#             } else {
#               x <- callNextMethod(x = x, ...)
#             }
#             
#             x
#           })
# 
# 
# setMethod("[[", "annotationList",
#           function (x, i, j, ...) {
#             if (!missing(i)) {
#               if (is.numeric(i)) {
#                 i <- names(x)[i]
#               }
#               x <- callNextMethod(x = x, i = i, ...)
#               seql <- which(seqlevels(x) %in% i)
#               seqinfo(x, new2old = seql) <- seqinfo(x)[i]
#             } else {
#               x <- callNextMethod(x = x, ...)
#             }
#             
#             x
#           })


#' @export
setGeneric("product", function (x, ...) standardGeneric("product"))
setMethod("product", "annotationList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "product"))
          })

setMethod("product", "GRanges",
          function (x) {
            pt <- mcols(x)[["product"]]
            if (is.null(pt))
              message("No 'product' metadata available")
            pt
          })


#' @export
setGeneric("locusTag", function (x, ...) standardGeneric("locusTag"))
setMethod("locusTag", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "synonym"))
          })

setMethod("locusTag", "GRanges",
          function (x) {
            lt <- mcols(x)[["synonym"]]
            if (is.null(lt))
              message("No 'locus_tag' metadata available")
            lt
          })


#' @export
setGeneric("geneID", function (x, ...) standardGeneric("geneID"))
setMethod("geneID", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "geneID"))
          })

setMethod("geneID", "GRanges",
          function (x) {
            mcols(x)[["geneID"]]
          })


#' @export
setGeneric("proteinID", function (x, ...) standardGeneric("proteinID"))
setMethod("proteinID", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "proteinID"))
          })

setMethod("proteinID", "GRanges",
          function (x) {
            mcols(x)[["proteinID"]]
          })


#' @export
setGeneric("type")
setMethod("type", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "type"))
          })

setMethod("type", "GRanges",
          function (x) {
            mcols(x)[["type"]]
          })
