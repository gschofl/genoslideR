\name{annotatedAlignment}
\alias{annotatedAlignment}
\alias{annotatedAlignment-class}
\title{annotatedAlignment}
\usage{
  annotatedAlignment(aln, anno, type, ...)
}
\arguments{
  \item{aln}{Path to a directory containing mercator
  segments, a single file in mfa (multi fasta) or maf
  (multiple alignment file) format containing the genome
  aligment.}

  \item{anno}{Paths to annotation files.}

  \item{type}{Type of annotation files.}

  \item{...}{seqid (accession numbers)}
}
\description{
  \dQuote{annotatedAlignment} is an S4 class that provides
  a container for multiple alignments stored as an
  \dQuote{\linkS4class{XStringSet}} and genomic annotations
  stored as an \dQuote{\linkS4class{annotationList}}.

  Construct an annotated alignment
}

