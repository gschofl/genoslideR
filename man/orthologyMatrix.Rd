% Generated by roxygen2 (4.0.0): do not edit by hand
\name{orthologyMatrix}
\alias{orthologyMatrix}
\title{Construct orthology matrices}
\usage{
orthologyMatrix(blat_dir, base_names = NULL, pident_threshold = 60,
  paln_threshold = 40, name_sep = "-", log = TRUE)
}
\arguments{
\item{blat_dir}{Directory where the anchor and hit files produced by
\code{link{mercator}} live.}

\item{base_names}{Species for which an orthology matrix should be
constructed. Must correspond to the relevant filenames minus extension.
If \code{NULL} all are selected}

\item{pident_threshold}{Threshold for (cumulative) percentage of sequence
identity}

\item{paln_threshold}{Threshhold for (cumulative) percentage of alignment
length (query coverage).}

\item{name_sep}{String that separates two compared species in filenames.}

\item{log}{Write logfiles}
}
\value{
A list of class \sQuote{orthoMatrix}, which is essentially a
list of \code{\link[Matrix]{dgCMatrix}}es with attributes
\sQuote{labels}, \sQuote{nGenes}, and \sQuote{combinations} attached.
}
\description{
Constructs pairwise orthology matrices for all combinations of species
as provided by the \code{base_names} argument based on reciprocal best
BLAT hits. Use \code{\link{mergeOrthologyMatrix}} to merge these
pairwise comparison into a single data frame.
}

