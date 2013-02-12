\name{get_intergenic_ranges}
\alias{get_intergenic_ranges}
\title{Generate intergenic regions}
\usage{
  get_intergenic_ranges(gff = chps_gff)
}
\arguments{
  \item{gff}{A gff annotation file to calculate the
  intergenic regions.}
}
\value{
  An IRanges object including the start and end plus the
  intergenic name of each intergenic region.
}
\description{
  Generate intergenic regions
}
\details{
  Zero length ranges and unknown locus_tags are excluded.
}

