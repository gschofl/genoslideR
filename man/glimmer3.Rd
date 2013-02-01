\name{glimmer3}
\alias{glimmer3}
\title{Ab initio genome annotation using glimmer3}
\usage{
  glimmer3(fasta, opts = list(o = 50, g = 110, t = 30),
    ..., cleanup = TRUE)
}
\arguments{
  \item{fasta}{Path to a genome file(s) in fasta format}

  \item{opts}{a named list of options for glimmer3}

  \item{...}{named values interpreted as options for
  glimmer3}

  \item{cleanup}{Clean up intermediate files.}
}
\value{
  Character vector. Path to glimmer3 file(s).
}
\description{
  Ab initio genome annotation using glimmer3
}

