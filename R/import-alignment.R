#' @importFrom Biostrings readDNAStringSet
#' @importFrom parallel mclapply
#' @importFrom rmisc %.%
#' @importFrom IRanges runValue
#' @importFrom IRanges Map
NULL

#' Import a genome alignment
#' 
#' @param aln Path to a directory containing mercator segments,
#' a single file in mfa (mult-fasta) or maf (multiple alignment file)
#' format containing the genome aligment.
#'
#' @details
#' If the alignment is provided as a multi fasta file the each header
#' must contain a map of the orthologous segments in the following 
#' format: 
#'  
#' GenomeID ChromosomeID:start:end:strand ... ChromosomeID:start:end:strand
#'
#' where strand is encodes as + or -
#'  
#' @return A \code{\linkS4class[Biostrings]{DNAStringSet}} instance.
#' @export
importAlignment <- function (aln) {
  
  if (!file.exists(aln)) {
    stop("No alignment file provided")
  }
  
  if (is_segments_dir(aln)) {
    seq <- merge_mercator_segments(seg_dir=segments_dir(aln))
    aln_path <- normalizePath(metadata(seq)[["path"]])
  } else if (is_maf(aln)) {
    seq <- readMAF(aln)
    aln_path <- normalizePath(metadata(seq)[["path"]])
    on.exit(unlink(aln_path))
  } else if (is_mfa(aln)) {
    seq <- readDNAStringSet(aln)
    aln_path <- aln
  } else {
    stop("No valid alignment file provided")
  }
  
  headers <- names(seq)
  metadata(seq) <- header2map(headers, aln_path)
  names(seq) <- strsplitN(headers, ' ', 1)
  seq
}


header2map <- function (headers, aln_path) {
  if (is.list(headers)) {
    headers <- unlist(unname(headers))
  }
  genomes <- strsplitN(headers, ' ', 1)
  gapranges <- get_alignment_gaps(aln_path, genomes)
  
  map <- strsplit(strsplitN(headers, ' ', -1), ' ')
  map <- lapply(map, function (m) do.call(rbind, strsplit(m, ":"))) 
  map <- lapply(map, as.data.frame.matrix, stringsAsFactors=FALSE)
  map <- lapply(map, `colnames<-`, value=c("chr", "genomic_start",
                                           "genomic_end", "strand"))
  for (i in seq_along(map)) {
    map[[i]][,c("genomic_start")] <- as.integer(map[[i]][,c("genomic_start")])
    map[[i]][,c("genomic_end")] <- as.integer(map[[i]][,c("genomic_end")])
  }
  
  genomicMap <- GRangesList(Map(function(g, m) {
    GRanges(seqnames=g, IRanges(start=m[["genomic_start"]],
                                end=m[["genomic_end"]],
                                names=m[["chr"]]),
            strand=m[["strand"]])
  }, g = as.list(genomes), m = map))
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


get_alignment_gaps <- function (aln_path, genomes) {
  #aln_path <- "~/daten/alignment/aln.mfa"
  exec <- system.file("src", "runlvec.pl", package="genoslideR")
  #exec <- "~/R/Projects/Devel/genoslideR/inst/src/runlvec.pl"
  cores <- detectCores()
  gap_ranges <- GRangesList(mclapply(seq_along(genomes), function (genome_seq) {
    pipe_desc <- pipe(paste(exec, aln_path, genome_seq))
    l <- scan(pipe_desc, sep="\t", quiet=TRUE,
              what=list(pos = integer(), len = integer()))
    close(pipe_desc)
    GRanges(seqnames=genomes[genome_seq],
            ranges=IRanges(start=l[["pos"]], width=l[["len"]]),
            strand="*")
  }, mc.cores = cores - 1))
  names(gap_ranges) <- genomes
  gap_ranges
}

