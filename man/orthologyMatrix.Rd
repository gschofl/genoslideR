\name{orthologyMatrix}
\alias{orthologyMatrix}
\title{Construct orthology matrices}
\usage{
  orthologyMatrix(base_names, blat_dir, cutoff = 0.001,
    name_sep = "-", log = TRUE)
}
\arguments{
  \item{base_names}{Species for which an orthology matrix
  should be constructed. Must correspond to the relevant
  filenames minus extension.}

  \item{blat_dir}{Directory where the anchor and hit files
  produced by \code{link{mercator}} live.}

  \item{cutoff}{Evalue cutoff.}

  \item{name_sep}{String that separates two compared
  species in filenames.}

  \item{log}{Write logfiles}
}
\value{
  A list of class \sQuote{orthoMatrix}, which is
  essentially a list of \code{\link[Matrix]{dgCMatrix}}es
  with attributes \sQuote{labels}, \sQuote{nGenes}, and
  \sQuote{combinations} attached.
}
\description{
  Constructs pairwise orthology matrices for all
  combinations of species as provided by the
  \code{base_names} argument based on reciprocal best BLAT
  hits. Use \code{\link{mergeOrthologyMatrix}} to merge
  these pairwise comparison into a single data frame.
}

