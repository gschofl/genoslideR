\name{genome2Alignment}
\alias{genome2Alignment}
\title{Map genomic positions to alignment positions.}
\usage{
  genome2Alignment(ranges, aln, simplify = FALSE)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} containing genomic
  positions to map to the alignment.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{simplify}{if \code{TRUE}, return \code{GRanges},
  otherwise return a \code{GRangesList}.}
}
\value{
  A \code{\linkS4class{GRangesList}} object. object.
}
\description{
  Map genomic positions to alignment positions.
}
\seealso{
  \code{\link{alignment2Genome}},
  \code{\link{genome2Genome}},
  \code{\link{sliceAlignment}},
  \code{\link{findAnnotation}}.
}

