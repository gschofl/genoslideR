\name{sliceAlignment}
\alias{sliceAlignment}
\title{Slice genomic ranges from an alignment.}
\usage{
  sliceAlignment(ranges, aln, targetGenomes = NULL)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} containing genomic
  positions to be mapped to the alignment.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{targetGenomes}{}
}
\value{
  A \code{\linkS4class{DNAStringSetList}} object.
}
\description{
  The ranges passed to \code{sliceAlignment} are first
  converted from genomic ranges to alignment ranges. The
  alignment ranges are returned as metadata with the
  resulting \code{\linkS4Class{DNAStringSetList}} object.
}

