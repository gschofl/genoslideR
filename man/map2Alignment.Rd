\name{map2Alignment}
\alias{map2Alignment}
\title{Map genomic positions to alignment positions.}
\usage{
  map2Alignment(ranges, aln, simplify = TRUE)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} object
  containing genomic positions to be mapped to the
  alignment.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{simplify}{if \code{TRUE}, return a \code{GRanges}
  instance.}
}
\value{
  A \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} object.
}
\description{
  Map genomic positions to alignment positions.
}

