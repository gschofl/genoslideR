\name{genome2Genome}
\alias{genome2Genome}
\title{Map genomic positions in one genome to genomic positions in other genomes.}
\usage{
  genome2Genome(ranges, aln, targetGenomes = NULL)
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
  Map genomic positions in one genome to genomic positions
  in other genomes.
}

