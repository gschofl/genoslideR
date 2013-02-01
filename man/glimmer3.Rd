\name{glimmer3}
\alias{glimmer3}
\title{Ab initio genome annotation using glimmer3}
\usage{
  glimmer3(genome = seq_files[1],
    opts = list(o = 50, g = 110, t = 30), ...,
    cleanup = TRUE)
}
\arguments{
  \item{genome}{Path to a genome file in fasta format}

  \item{opts}{a named list of options for glimmer3}

  \item{...}{named values interpreted as options for
  glimmer3}

  \item{cleanup}{Clean up intermediate files.}
}
\description{
  Ab initio genome annotation using glimmer3
}

