\name{sliceAlignment}
\alias{sliceAlignment}
\title{Extract genomic ranges from an alignment.}
\usage{
  sliceAlignment(ranges, aln, addAnnotations = TRUE,
    targetGenomes = NULL)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} object
  containing genomic positions to be mapped to the
  alignment.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{addAnnotations}{}

  \item{targetGenomes}{}
}
\value{
  A \code{\linkS4class{GRangesList}} object.
}
\description{
  Extract genomic ranges from an alignment.
}

