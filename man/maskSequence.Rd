\name{maskSequence}
\alias{maskSequence}
\title{Softmask repeats and low-complexity regions using
RepeatScout and RepeatMasker}
\usage{
  maskSequence(fasta)
}
\arguments{
  \item{fasta}{Path to fasta file(s).}
}
\value{
  Character vector. Path to softmasked file(s).
}
\description{
  \code{maskSequence} will create a hidden directory
  \sQuote{.masked} in the parent directory of the submitted
  fasta files and place all all intermediate files as well
  as the final softmasked sequence file into this
  directory.
}
\details{
  If multiple fasta files share the same parent directory,
  the softmasked files are placed in individual
  subdirectories within \sQuote{.masked}.
}

