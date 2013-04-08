\name{alignSegments}
\alias{alignSegments}
\title{Align orthologous genome segments with \code{fsa}}
\usage{
  alignSegments(seg_dir,
    fsa.opts = list(anchored = TRUE, translated = TRUE, exonerate = TRUE, softmasked = TRUE),
    ..., skip.completed = TRUE, ncores = detectCores() - 1)
}
\arguments{
  \item{seg_dir}{Path to mercator segments.}

  \item{fsa.opts}{a named list of options for
  \code{\link{fsa}}. The default options are optimisiing
  fsa for long sequence alignment.}

  \item{...}{Named values interpreted as options for fsa.}

  \item{skip.completed}{if \code{TRUE}, don't realign an
  existing alignment.}

  \item{ncores}{Number of cores to utilised for
  parallelisation}
}
\description{
  Align orthologous genome segments with \code{fsa}
}
\seealso{
  \code{\link{mercator}}, \code{\link{annotatedAlignment}}.
}

