\name{annotatedAlignment}
\alias{annotatedAlignment}
\title{annotatedAlignment}
\usage{
  annotatedAlignment(aln, anno, type, ...)
}
\arguments{
  \item{aln}{Path to a multifasta alignment file or a
  mercator directory containing aligned orthologous
  segments.}

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

