#' @importFrom Biostrings readDNAStringSet
#' @importFrom parallel mclapply
#' @importFrom rmisc %.%
#' @importFrom IRanges runValue
#' @importFrom IRanges Map
#' @importFrom GenomicRanges seqlengths<-
NULL

#' Import a genome alignment
#' 
#' @param alnpath Path to the directory containing mercator segments,
#' a single file in \sQuote{mfa} (multi-fasta) or \sQuote{maf} (multiple
#' alignment file) format containing the genome aligment.
#' @param outfile if not \code{NULL} the alignment is writen to file in
#' multi-fasta format.
#'
#' @details
#' If the alignment is provided as a multi-fasta file the each header
#' must contain a map of the orthologous segments in the following 
#' format: 
#'  
#' GenomeID GenomeLength ChromosomeID:start:end:strand ... ChromosomeID:start:end:strand
#'
#' where strand is encodes as + or - and GenomeLength is the length of the
#' genomes in the alignment which can be shorter than the length of the
#' original genome.
#'  
#' @return A \code{\linkS4class{DNAStringSet}} instance.
#' @seealso \code{\link{mercator}}.
#' @export
importAlignment <- function (alnpath, outfile=NULL) {
  
  if (!file.exists(alnpath)) {
    stop("No alignment file provided")
  }
  
  if (is_segments_dir(alnpath)) {
    seq <- merge_mercator_segments(segments_dir(alnpath), outfile)
  } else if (is_maf(alnpath)) {
    seq <- readMAF(alnpath)
  } else if (is_mfa(alnpath)) {
    seq <- readDNAStringSet(alnpath)
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
  
  gapranges <- get_alignment_gaps(x=seq, genomes)
  seqlengths(gapranges) <- width(seq)
  
  map <- strsplit(strsplitN(headers, ' ', -c(1:2)), ' ')
  map <- lapply(map, function (m) do.call(rbind, strsplit(m, ":"))) 
  map <- lapply(map, as.data.frame.matrix, stringsAsFactors=FALSE)
  map <- lapply(map, `colnames<-`, value=c("chr", "genomic_start",
                                           "genomic_end", "strand"))
  for (i in seq_along(map)) {
    map[[i]][["genomic_start"]] <-
      suppressWarnings(as.integer(map[[i]][["genomic_start"]]))
    map[[i]][["genomic_end"]] <-
      suppressWarnings(as.integer(map[[i]][["genomic_end"]]))
    map[[i]][["strand"]] <-
      ifelse(map[[i]][["strand"]] == 'NA', '*', map[[i]][["strand"]])
  }
  
  genomicMap <- GRangesList(Map(function(g, cov, m) {
    seqinfo <- Seqinfo(g, as.integer(cov))
    GRanges(seqnames=g, map2range(name=m[["chr"]],
                                  start=m[["genomic_start"]],
                                  end=m[["genomic_end"]]),
            strand=m[["strand"]], seqinfo=seqinfo)
  }, g = as.list(genomes), cov = as.list(coverage), m = map))
  names(genomicMap) <- IRanges::Map(as.character%.%runValue, seqnames(genomicMap))
  
  alignmentMap <- GRangesList(Map(update_position, m = map, g = as.list(genomes)))
  names(alignmentMap) <- IRanges::Map(as.character%.%runValue, seqnames(alignmentMap))
  
  list(gMap=genomicMap, aMap=alignmentMap, gaps=gapranges) 
}


map2range <- function(name, start, end) {
  width <- end - start + 1L
  start <- ifelse(is.na(start), 1L, start)
  width <- ifelse(is.na(width), 0L, width)
  new("IRanges", start=start, width=width, NAMES=name)
}


update_position <- function (m, g) {
  s <- 1L
  w <- m[["genomic_end"]] - m[["genomic_start"]] + 1L
  w <- ifelse(is.na(w), 0L, w)
  mat <- matrix(rep(0L, nrow(m)*2), ncol=2)
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
  }, genome=genomes, seq=seqs))
}
