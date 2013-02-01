#' annotationList
#' 
#' \dQuote{annotationList} is an S4 class that provides a container for
#' genomic annotations stored in \dQuote{\linkS4class{gbRange}}s.
#'
#' @export
setClass("annotationList", contains="GRangesList")


annotationList <- function(files = anno, type = "genbank", ...) {
  
  type <- match.arg(type, c("genbank", "gff", "ptt", "ftable", "table"))
  
  if (type == "genbank") {
    stop("Not implemented yet")
  }
  
  if (type == "gff") {
    anno <- GRangesList(lapply(files, import_annotation_from_gff))
    names(anno) <- seqnames(seqinfo(anno))
  }
  
  if (type == "ptt") {
    if (is.null(seqid <- list(...)[["seqid"]])) {
      anno <- GRangesList(lapply(files, import_annotation_from_ptt))
      
    } else {
      if (length(seqid) != length(files))
        stop("Provide as many sequence identifiers (accession numbers) as annotation files")
      anno <- GRangesList(mapply(import_annotation_from_ptt, file=files,
                                 seqid=seqid, SIMPLIFY=FALSE, USE.NAMES=FALSE))
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
      anno <- GRangesList(mapply(import_annotation_from_tbl, file=files, 
                                 seqid=seqid, sep=sep, SIMPLIFY=FALSE,
                                 USE.NAMES=FALSE))
    }
    names(anno) <- seqnames(seqinfo(anno))
  } 
  
  if (type == "ftable") {
    anno <- GRangesList(lapply(files, import_annotation_from_ftb))
    names(anno) <- seqnames(seqinfo(anno))
  } 
  
  new("annotationList", anno)
}


setMethod("accession", "annotationList",
          function (x) seqnames(seqinfo(x)))

setMethod("definition", "annotationList",
          function (x) unname(genome(seqinfo(x))))

setMethod("[[", "annotationList",
          function (x, i, j, ...) {
            callNextMethod()
          })

setMethod("[", "annotationList",
          function (x, i, j, ..., drop = TRUE) {
            callNextMethod()
          })

setMethod("length", "annotationList",
          function (x) callNextMethod())

setMethod("names", "annotationList",
          function (x) callNextMethod())

setMethod("strand", "annotationList",
          function (x) callNextMethod())

setMethod("mcols", "annotationList",
          function (x) callNextMethod())


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


setGeneric("product", function (x, ...) standardGeneric("product"))
setGeneric("locusTag", function (x, ...) standardGeneric("locusTag"))
setGeneric("geneID", function (x, ...) standardGeneric("geneID"))
setGeneric("proteinID", function (x, ...) standardGeneric("proteinID"))


setMethod("product", "GRangesList",
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

setMethod("geneID", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "geneID"))
          })

setMethod("geneID", "GRanges",
          function (x) {
            mcols(x)[["geneID"]]
          })

setMethod("proteinID", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "proteinID"))
          })

setMethod("proteinID", "GRanges",
          function (x) {
            mcols(x)[["proteinID"]]
          })

setMethod("type", "GRangesList",
          function(x) {
            elm <- elementMetadata(x, level="within")
            CharacterList(lapply(elm, `[[`, "type"))
          })

setMethod("type", "GRanges",
          function (x) {
            mcols(x)[["type"]]
          })


#' annotatedAlignment
#' 
#' \dQuote{annotatedAlignment} is an S4 class that provides a container for
#' multiple alignments stored as an \dQuote{\linkS4class{XStringSet}} and
#' genomic annotations stored as an \dQuote{\linkS4class{annotationList}}.
#'
#' @export
setClass("annotatedAlignment",
         representation(annotation="annotationList",
                        alignment="XStringSet"))


#' Construct an annotated alignment
#' 
#' @param aln Path to a directory containing mercator segments,
#' a single file in mfa (multi fasta) or maf (multiple alignment file)
#' format containing the genome aligment.
#' @param anno Paths to annotation files.
#' @param type Type of annotation files.
#' @param ... seqid (accession numbers)
#' 
#' @export
annotatedAlignment <- function (aln, anno, type, ...) {
  
  
  if (missing(aln)) {
    stop("No alignment provided")
  }

  if (missing(anno)) {
    stop("No annotation provided")
  }
  
  if (missing(type)) {
    stop("No annotation type provided")
  }
  
  alignment <- if(!is_mapped_alignment(aln)) {
    importAlignment(aln)
  }
    
  genomes <- names(alignment)
  alignment <- alignment[order(genomes)]
  genome_annotations <- strip_ext(basename(anno))
  
  idx <- genomes %in% genome_annotations
  if (all(idx == FALSE)) {
    stop("The annotation file names don't match with the genome names in the alignment")
  }
  
  anno <- sort(anno[idx])
  annotation <- annotationList(anno, type, ...)
  
  res <- new("annotatedAlignment", annotation = annotation, alignment = alignment)
  res
  
}


setValidity("annotatedAlignment", function (object) {
  
  n <- names(object)
  if (!all(n$annotation %in% n$alignment)) {
    return("Name mismatch between annotations and alignments")
  }

  met <- metadata(alignment(object))
  if (!all(idx <- c("map", "gaps") %in% names(met))) {
    return(sprintf("%s missing from alignment metadate", c("map", "gaps")[idx]))
  }

  if (!all(names(met$gaps) %in% n$alignment)) {
    return("Name mismatch between alignments and gap metadata")
  }
  
  if (!all(names(met$map) %in% n$alignment)) {
    return("Name mismatch between alignments and map metadata")
  }
  
  if (any(vapply(met$gaps, class, character(1), USE.NAMES=FALSE) != "IRanges")) {
    return("Gaps must be of class 'IRanges'")
  }
   
  if (any(vapply(met$map, class, character(1), USE.NAMES=FALSE) != "data.frame")) {
    return("Map must be of class 'data.frame'")
  }

  TRUE
})




setMethod("annotation", "annotatedAlignment",
          function (object) object@annotation)


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

# setMethod("[", "annotatedAlignment",
#           function (x, i, j, ..., drop = TRUE) {
#             # ignore j
#             if (!missing(j)) {
#               stop("invalid subsetting")
#             }
#             aln <- alignment(x)
#             met <- metadata(aln)
#             ann <- annotation(x)
#             
#             if (is.numeric(i)) {
#               i <- names(aln)[i]
#             }
#             
#             aln <- callNextMethod(x = aln, i = i, ... = ..., drop = drop)
#             metadata(aln) <- list(map = met[["map"]][i], gaps = met[["gaps"]][i])
#             
#             i <- names(ann)[names(ann) %in% i]
#             ann <- if (is_empty(i)) {
#               new("annotationList")
#             } else {
#               callNextMethod(x = ann, i = i, ... = ..., drop = drop)
#             }
#             
#             new("annotatedAlignment", annotation = ann,
#                 alignment = aln)
#           })

setMethod("length", "annotatedAlignment",
          function (x) {
            c(alignment = length(x@alignment),
              annotation =  length(x@annotation))
          })

setMethod("names", "annotatedAlignment",
          function (x) {
            list(alignment = names(x@alignment),
                 annotation =  names(x@annotation))
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

setMethod("seqnames", "annotatedAlignment",
          function (x) {
            seqnames(annotation(x))
          })

setMethod("strand", "annotatedAlignment",
          function (x, ...) {
            GenomicRanges::strand(annotation(x))
          })

setMethod("seqlevels", "annotatedAlignment",
          function (x) {
            seqlevels(annotation(x))
          })

setMethod("seqinfo", "annotatedAlignment",
          function (x) {
            seqinfo(annotation(x))
          })

setMethod("seqlengths", "annotatedAlignment",
          function (x) {
            seqlengths(annotation(x))
          })

setMethod("genome", "annotatedAlignment",
          function (x) {
            genome(annotation(x))
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

setGeneric("getMap", function (x, ...) standardGeneric("getMap"))
setMethod("getMap", "annotatedAlignment",
          function (x, i = NULL) {
            ans <- metadata(alignment(x))[["map"]]
            if (!is.null(i)) {
              ans <- ans[i]
            }
            return(ans)
          })

setGeneric("getGaps", function (x, ...) standardGeneric("getGaps"))
setMethod("getGaps", "annotatedAlignment",
          function (x, i = NULL) {
            ans <- metadata(alignment(x))[["gaps"]]
            if (!is.null(i)) {
              ans <- ans[i]
            }
            return(ans)
          })

