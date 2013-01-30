\name{mercator}
\alias{mercator}
\title{Run mercator to build a homology map and create orthologous segments
that can be aligned using FSA}
\usage{
  mercator(seq_files, anno_files = NULL,
    anno_type = "glimmer3", wd = NULL, mask = TRUE)
}
\arguments{
  \item{seq_files}{Path to sequence files ('fasta' or
  'genbank').}

  \item{anno_files}{(Optional) Path to annotation files. If
  no annotation files are provided an \emph{ab initio}
  annotation using \code{\link{glimmer3}} is performed.}

  \item{anno_type}{Type of annotation ('glimmer3',
  'genbank', 'gff', 'ptt', or 'ftable').}

  \item{wd}{Working directory. Defaults to the parent
  directory of the provided sequence files.}

  \item{mask}{Softmask sequence before aligning using
  \code{\link{maskSequence}}.}
}
\description{
  Run mercator to build a homology map and create
  orthologous segments that can be aligned using FSA
}
\details{
  \code{mercator} will create a hidden directory
  \sQuote{.mercator} in the parent directory of the
  submitted sequence files and place all results of the
  intermediate computational steps in this directory.

  The results of a \code{mercator} run are stored in a
  directory \sQuote{segments} in the parent directory of
  the submitted sequence files. The \sQuote{segments}
  directory should contain numbered subdirectories
  containing orthologous genomic segments in multi fasta
  format, a file named \sQuote{map}, and a file named
  \sQuote{genomes}.

  \code{mercator} segements can be aligned using
  \code{\link{alignMercatorSegments}}, merged into a single
  multi-fasta alignment file using
  \code{\link{importAlignment}} or imported to R as an
  \code{\linkS4class{annotatedAlignment}} using
  \code{\link{annotatedAlignment}}.
}

