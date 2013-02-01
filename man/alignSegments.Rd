\name{alignSegments}
\alias{alignSegments}
\title{Align orthologous genome segments with \code{fsa}}
\usage{
  alignSegments(seg_dir, force = FALSE,
    opts = list(anchored = TRUE, translated = TRUE, exonerate = TRUE, softmasked = TRUE),
    ...)
}
\arguments{
  \item{seg_dir}{Path to mercator segments.}

  \item{force}{Realign an existing alignment.}

  \item{opts}{a named list of options for
  \code{\link{fsa}}. The default options are optimizes for
  long sequences.}

  \item{...}{Named values interpreted as options for fsa.}
}
\description{
  Align orthologous genome segments with \code{fsa}
}

