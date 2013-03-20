\name{alignment2Genome}
\alias{alignment2Genome}
\title{Map alignment positions to genomic positions.}
\usage{
  alignment2Genome(ranges, aln, targetGenomes = NULL)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} object
  containing alignment positions to be mapped to the
  genomes.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{targetGenomes}{}
}
\value{
  A \code{\linkS4class{GRangesList}} object.
}
\description{
  Map alignment positions to genomic positions.
}

