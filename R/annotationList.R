#' @importClassesFrom GenomicRanges GenomicRangesORGRangesList
#' @importClassesFrom GenomicRanges GenomicRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GenomicRangesList
#' @importClassesFrom GenomicRanges GRangesList
#' @importClassesFrom GenomicRanges Seqinfo
#' @importClassesFrom IRanges CompressedList
#' @importClassesFrom IRanges Annotated
#' @importClassesFrom IRanges List
#' @importClassesFrom IRanges Vector
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges Seqinfo
#' @importFrom GenomicRanges seqinfo
#' @importFrom GenomicRanges "seqinfo<-"
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges genome
#' @importFrom GenomicRanges seqlengths
#' @importFrom GenomicRanges seqlevels
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges mcols
#' @importFrom IRanges CharacterList
#' @importFrom IRanges IntervalTree
#' @importFrom IRanges metadata
#' @importFrom IRanges metadata<-
#' @importFrom IRanges elementMetadata
#' @importFrom IRanges elementMetadata<-
#' @importFrom IRanges Rle
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @importFrom IRanges width
#' @importFrom IRanges elementLengths
#' @importFrom IRanges lapply
#' @importFrom IRanges gaps
#' @importFrom Biostrings type
#' @importFrom biofiles accession
#' @importFrom biofiles definition
#' @importFrom biofiles qualif
NULL

#' annotationList
#' 
#' \dQuote{annotationList} is an S4 class that provides a container for
#' genomic annotations stored in \dQuote{\linkS4class{gbRange}}s.
#'
#' @export
setClass("annotationList", contains="GRangesList")


annotationList <- function(files, type = "gff", ...) {
  
  type <- match.arg(type, c("gff", "ptt", "ftb", "table"))
  importer <- match.fun(paste0("import_annotation_from_", type))
  
  if (type %in% c("gff", "ftable")) {
    # import_annotation_from_gff
    features <- list(...)[["features"]] %|null|% c("CDS", "RNA")
    anno <- GRangesList(lapply(files, importer, features=features))
    names(anno) <- seqnames(seqinfo(anno))
  }
  
  if (type == "ptt") {
    if (is.null(seqid <- list(...)[["seqid"]])) {
      anno <- GRangesList(lapply(files, importer))
    } else {
      if (length(seqid) != length(files))
        stop("Provide as many sequence identifiers (accession numbers) as annotation files")
      anno <- GRangesList(mapply(importer, file=files, seqid=seqid,
                                 SIMPLIFY=FALSE, USE.NAMES=FALSE))
    }
    names(anno) <- seqnames(seqinfo(anno))
  }
  
  if (type == "table") {
    if (is.null(seqid <- list(...)[["seqid"]])) {
      stop("Provide sequence identifiers (accession numbers) for the annotation files.")
    } else {
      if (length(seqid) != length(files))
        stop("Provide as many sequence identifiers as annotation files")
      sep <- list(...)[["sep"]] %||% "\t"
      anno <- GRangesList(mapply(importer, file=files, seqid=seqid,
                                 sep=sep, SIMPLIFY=FALSE, USE.NAMES=FALSE))
    }
    names(anno) <- seqnames(seqinfo(anno))
  } 
  
  new("annotationList", anno)
}


.show_AnnotationList <- function(x) {
  len <- length(x)
  cat(sQuote(class(x)), " of length ", len, "\n", sep = "")
  if (len == 0L)
    return(invisible(NULL))
  nm <- names(x)
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
setGeneric("accession")
setMethod("accession", "annotationList",
          function (x) seqnames(seqinfo(x)))


#' @export
setGeneric("definition")
setMethod("definition", "annotationList",
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


setMethod("[", "annotationList",
          function (x, i, j, ..., drop = TRUE) {
            if (!missing(i)) {
              if (is.numeric(i)) {
                i <- names(x)[i]
              }
              x <- callNextMethod(x = x, i = i, ...)
              seql <- which(seqlevels(x) %in% i)
              seqinfo(x, new2old = seql) <- seqinfo(x)[i]
            } else {
              x <- callNextMethod(x = x, ...)
            }
            
            x
          })


setMethod("[[", "annotationList",
          function (x, i, j, ...) {
            if (!missing(i)) {
              if (is.numeric(i)) {
                i <- names(x)[i]
              }
              x <- callNextMethod(x = x, i = i, ...)
              seql <- which(seqlevels(x) %in% i)
              seqinfo(x, new2old = seql) <- seqinfo(x)[i]
            } else {
              x <- callNextMethod(x = x, ...)
            }
            
            x
          })


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
