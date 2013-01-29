\name{ncbi_bacteria}
\alias{ncbi_bacteria}
\title{Retrieve bacterial genomes from NCBI FTP site}
\usage{
  ncbi_bacteria(which, what = "gbk|gff|fna",
    where = "~/Bacteria", ignore.case = TRUE)
}
\arguments{
  \item{which}{A regular expression to match the
  directories of completed bacterial genomes.}

  \item{what}{A regular expression to match file
  extensions. Possible values are \emph{gbk} (Genbank),
  \emph{gff} (GFF3 file), \emph{fna} (FASTA), \emph{ptt}
  (protein table).}

  \item{where}{Local target directory for downloading}

  \item{ignore.case}{Case sensitive pattern matching if
  \code{FALSE}}
}
\description{
  Get genomes from
  \url{ftp://ftp.ncbi.nih.gov/genomes/Bacteria/}
}

