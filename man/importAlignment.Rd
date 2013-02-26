\name{importAlignment}
\alias{importAlignment}
\title{Import a genome alignment}
\usage{
  importAlignment(aln, writeToFile = FALSE, outfile = NULL)
}
\arguments{
  \item{aln}{Path to a directory containing mercator
  segments, a single file in mfa (mult-fasta) or maf
  (multiple alignment file) format containing the genome
  aligment.}

  \item{writeToFile}{}

  \item{outfile}{}
}
\value{
  A \code{\linkS4class[Biostrings]{DNAStringSet}} instance.
}
\description{
  Import a genome alignment
}
\details{
  If the alignment is provided as a multi fasta file the
  each header must contain a map of the orthologous
  segments in the following format:

  GenomeID GenomeLength ChromosomeID:start:end:strand ...
  ChromosomeID:start:end:strand

  where strand is encodes as + or - and GenomeLength is the
  length of the genomes in the alignment which can be
  shorter than the length of the original genome.
}

