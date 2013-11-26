\name{mergeOrthologyMatrix}
\alias{mergeOrthologyMatrix}
\title{Merge a list of pairwise orthology matrices into an OrthoFrame}
\usage{
  mergeOrthologyMatrix(omat)
}
\arguments{
  \item{omat}{An \code{orthoMatrix} object}
}
\value{
  An \code{OrthoFrame} object.
}
\description{
  Merges all pairwise comparisons of genomes as provided by
  \code{\link{orthologyMatrix}} into a single overall
  orthology Matrix. (an object of class
  \sQuote{\bold{OrthoFrame}}). It checks for conflicts and
  provides those as an attribute \sQuote{\emph{conflicts}}.
  The attribute \sQuote{\emph{labels}} holds the species
  names.
}

