\name{alignment2Annotation}
\alias{alignment2Annotation}
\title{Map alignment positions to annotations.}
\usage{
  alignment2Annotation(ranges, aln, targetGenome,
    type = "any")
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} object
  containing alignment positions.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{targetGenomes}{}
}
\value{
  A \code{\linkS4class{GRangesList}} object.
}
\description{
  Map alignment positions to annotations.
}

