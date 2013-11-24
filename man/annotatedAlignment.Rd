\name{annotatedAlignment}
\alias{annotatedAlignment}
\title{Construct an \sQuote{annotatedAlignment}.}
\usage{
  annotatedAlignment(alnpath, annopath, type, ...)
}
\arguments{
  \item{alnpath}{Path to a directory containing mercator
  segments, a single file in mfa (multi fasta) or maf
  (multiple alignment file) format containing the genome
  aligment.}

  \item{annopath}{Path to annotation files.}

  \item{type}{Type of annotation. Supported file types are
  \emph{gff}, \emph{ptt}, \emph{ftb} (NCBI feature table),
  \emph{genbank}/\emph{gbk}, and \emph{table} (custom
  tab-delimted annotation files, minimally containing
  \emph{start} and \emph{end} values).}

  \item{...}{Optionally \code{features} (default: CDS,
  RNA), \code{seqid} (sequence identifier, e.g. NCBI
  accession numbers), and \code{sep}.}
}
\description{
  \dQuote{annotatedAlignment} is an S4 class that provides
  a container for multiple alignments stored as an
  \dQuote{\linkS4class{XStringSet}} object and genomic
  annotations stored as an
  \dQuote{\linkS4class{annotationList}} object.
}
\seealso{
  \code{\link{mercator}}, \code{\link{alignSegments}}.
}

