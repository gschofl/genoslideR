\name{sliceAnnotatedAlignment}
\alias{sliceAnnotatedAlignment}
\title{Slice an Annotated Alignment}
\usage{
  sliceAnnotatedAlignment(aln, slices,
    addAnnotation = TRUE)
}
\arguments{
  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{slices}{a \code{\linkS4class{GRanges}} instance.}

  \item{addAnnotation}{if \code{TRUE} an
  \code{\linkS4class{annotatedAlignment}} is returned,
  otherwise a \code{\linkS4class{DNAStringSet}} object.}
}
\description{
  Slice an Annotated Alignment
}

