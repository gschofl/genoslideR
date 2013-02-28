#' @importFrom Biostrings readDNAStringSet
#' @importFrom parallel mclapply
#' @importFrom rmisc %.%
#' @importFrom IRanges runValue
#' @importFrom IRanges Map
#' @importFrom GenomicRanges seqlengths<-
NULL

#' Import a genome alignment
#' 
#' @param aln Path to a directory containing mercator segments,
#' a single file in mfa (mult-fasta) or maf (multiple alignment file)
#' format containing the genome aligment.
#' @param writeToFile
#' @param outfile
#'
#' @details
#' If the alignment is provided as a multi fasta file the each header
#' must contain a map of the orthologous segments in the following 
#' format: 
#'  
#' GenomeID GenomeLength ChromosomeID:start:end:strand ... ChromosomeID:start:end:strand
#'
#' where strand is encodes as + or - and GenomeLength is the length of the
#' genomes in the alignment which can be shorter than the length of the
#' original genome.
#'  
#' @return A \code{\linkS4class[Biostrings]{DNAStringSet}} instance.
#' @export
importAlignment <- function (aln, writeToFile=FALSE, outfile=NULL) {
  
  if (!file.exists(aln)) {
    stop("No alignment file provided")
  }
  
  if (is_segments_dir(aln)) {
    seq <- merge_mercator_segments(seg_dir=segments_dir(aln), writeToFile, outfile)
  } else if (is_maf(aln)) {
    seq <- readMAF(aln)
  } else if (is_mfa(aln)) {
    seq <- readDNAStringSet(aln)
  } else {
    stop("No valid alignment file provided")
  }
  
  map <- header2map(seq)
  metadata(seq) <- map
  names(seq) <- strsplitN(names(seq), ' ', 1)
  seq
}


header2map <- function (seq) {
  headers <- names(seq)
  genomes <- strsplitN(headers, ' ', 1)
  coverage <- strsplitN(headers, ' ', 2)
  
  gapranges <- get_alignment_gaps(seq, genomes)
  seqlengths(gapranges) <- width(seq)
  
  map <- strsplit(strsplitN(headers, ' ', -c(1:2)), ' ')
  map <- lapply(map, function (m) do.call(rbind, strsplit(m, ":"))) 
  map <- lapply(map, as.data.frame.matrix, stringsAsFactors=FALSE)
  map <- lapply(map, `colnames<-`, value=c("chr", "genomic_start",
                                           "genomic_end", "strand"))
  for (i in seq_along(map)) {
    map[[i]][,c("genomic_start")] <- as.integer(map[[i]][,c("genomic_start")])
    map[[i]][,c("genomic_end")] <- as.integer(map[[i]][,c("genomic_end")])
  }
  
  genomicMap <- GRangesList(Map(function(g, cov, m) {
    seqinfo <- Seqinfo(g, as.integer(cov))
    GRanges(seqnames=g, IRanges(start=m[["genomic_start"]],
                                end=m[["genomic_end"]],
                                names=m[["chr"]]),
            strand=m[["strand"]], seqinfo=seqinfo)
  }, g = as.list(genomes), cov = as.list(coverage), m = map))
  names(genomicMap) <- IRanges::Map(as.character%.%runValue, seqnames(genomicMap))
  
  alignmentMap <- GRangesList(Map(update_position, m = map, g = as.list(genomes)))
  names(alignmentMap) <- IRanges::Map(as.character%.%runValue, seqnames(alignmentMap))
  
  list(gMap=genomicMap, aMap=alignmentMap, gaps=gapranges) 
}


update_position <- function (m, g) {
  s <- 1
  w <- m[["genomic_end"]] - m[["genomic_start"]] + 1L
  mat <- matrix(rep(0, nrow(m)*2), ncol=2)
  for (i in seq_len(nrow(mat))) {
    mat[i,1] <- s
    mat[i,2] <- w[i] + s - 1L
    s <- mat[i,2] + 1L
  }
  GRanges(seqnames=g, IRanges(mat[,1], mat[,2], names=m[["chr"]]),
          strand="*")
}


readMAF <- function (maf = "~/local/workspace/Chlamydia/pomago.maf") {
  tempfile <- tempfile()
  # exec <- system.file("src", "maf2mfa.pl", package="genoslideR")
  exec <- "~/R/Projects/Devel/genoslideR/inst/src/maf2mfa.pl"
  res <- pipe(paste(exec, maf))
  write(scan(res, sep="\n", quiet=TRUE, what=character()),
        file = tempfile)
  close(res)
  seq <- readDNAStringSet(tempfile)
  metadata(seq) <- list(path=tempfile)
  seq
}


get_alignment_gaps <- function(x, genomes) {
  seqs <- as.list(as(x, "character"))
  GRangesList(mapply(function(genome, seq) {
    GRanges(genome, find_gaps(seq), strand="*")
  }, genomes, seqs))
}
