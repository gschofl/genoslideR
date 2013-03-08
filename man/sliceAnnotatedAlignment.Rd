\name{sliceAnnotatedAlignment}
\alias{sliceAnnotatedAlignment}
\title{Slice an Annotated Alignment}
\usage{
  sliceAnnotatedAlignment(ranges, aln,
    getAlignmentPositions = TRUE, getDNAStringSet = FALSE,
    getGenomicPositions = TRUE, getAnnotations = FALSE,
    targetGenomes = NULL)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} instance
  containing genomic positions to be mapped to the
  alignment.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{getAlignmentPositions}{if \code{TRUE}, return a
  \code{\linkS4class{GRangesList}} with genomic positions
  mapped to the alignment.}

  \item{getDNAStringSet}{if \code{TRUE}, return a
  \code{\linkS4class{DNAStringSetList}}.}

  \item{getGenomicPositions}{if \code{TRUE}}

  \item{getAnnotations}{if \code{TRUE}, return the sliced
  alignments as \code{\linkS4class{annotatedAlignment}}s.}

  \item{targetGenomes}{}
}
\description{
  Slice an Annotated Alignment
}

