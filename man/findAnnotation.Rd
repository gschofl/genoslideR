\name{findAnnotation}
\alias{findAnnotation}
\title{Find annotation associated with alignment positions.}
\usage{
  findAnnotation(ranges, aln, targetGenomes = NULL)
}
\arguments{
  \item{ranges}{a \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} containing alignment
  positions.}

  \item{aln}{an \code{\linkS4class{annotatedAlignment}}
  instance.}

  \item{targetGenomes}{a character vector of genome names.}
}
\value{
  A list of \code{\linkS4class{AnnotationList}} objects.
}
\description{
  Find annotation associated with alignment positions.
}
\seealso{
  \code{\link{alignment2Genome}},
  \code{\link{genome2Alignment}},
  \code{\link{genome2Genome}},
  \code{\link{sliceAlignment}}.
}

