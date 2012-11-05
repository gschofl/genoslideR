\name{mercator}
\alias{mercator}
\title{Run mercator to build a homology map and create orthologous segments
that can be aligned using FSA}
\usage{
  mercator(seq_files, seq_type = "fasta",
    anno_files = NULL, anno_type = "glimmer3", mask = TRUE)
}
\arguments{
  \item{seq_files}{Path to sequence files.}

  \item{seq_type}{Type of sequence file ('fasta' or
  'genbank').}

  \item{anno_files}{(Optional) Path to annotations.}

  \item{anno_type}{Type of annotation ('glimmer3',
  'genbank', 'gff', 'ptt', or 'feature_table').}

  \item{mask}{Softmask sequence before aligning.}
}
\description{
  Run mercator to build a homology map and create
  orthologous segments that can be aligned using FSA
}

