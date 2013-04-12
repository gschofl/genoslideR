\name{aggregateAlignment}
\alias{aggregateAlignment}
\title{Generate summaries on sliding windows over an alignment.}
\usage{
  aggregateAlignment(aln, ranges = NULL, window = NULL,
    step = NULL, f = "Eta", ...)
}
\arguments{
  \item{aln}{An \code{\linkS4class{annotatedAlignment}} or
  \code{\linkS4class{DNAStringSet}} object.}

  \item{ranges}{A \code{\linkS4class{GRanges}} object
  holding genomic windows on which summaries are
  generated.}

  \item{window}{Sliding window size. If no \code{ranges}
  are specified sliding windows are constructed over the
  complete alignment, otherwise sliding windows are
  constructed within each genomic window specified by
  \code{ranges}.}

  \item{step}{Step size of sliding windows.}

  \item{f}{Function applied to each window of \code{aln}.}

  \item{...}{Further arguments for \code{f}.}
}
\value{
  A \code{\linkS4class{RangedDataList}}.
}
\description{
  Generate summaries on sliding windows over an alignment.
}
\details{
  The function \code{f} applied to each window must accept
  a \code{\linkS4class{DNAStringSet}} as input and generate
  a named vector as output.
}

