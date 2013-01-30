\name{alignMercatorSegments}
\alias{alignMercatorSegments}
\title{Align orthologous genome segments produced by \code{mercator}}
\usage{
  alignMercatorSegments(seg_dir, aligner = "fsa",
    force = FALSE, mask = TRUE, ...)
}
\arguments{
  \item{seg_dir}{Path to mercator segments.}

  \item{aligner}{Which alignment program to use; 'fsa' or
  'mavid'.}

  \item{force}{Realign existing alignments.}

  \item{mask}{set \code{TRUE} if input sequences are
  softmasked.}

  \item{...}{Aligner options}
}
\description{
  Align orthologous genome segments produced by
  \code{mercator}
}

