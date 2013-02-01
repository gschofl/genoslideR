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
    aln <- normalizePath(metadata(seq)[["path"]])
  } else if (is_maf(aln)) {
    seq <- readMAF(aln)
    aln <- normalizePath(metadata(seq)[["path"]])
    on.exit(unlink(aln))
  } else if (is_mfa(aln)) {
    seq <- readDNAStringSet(aln)
  } else {
    stop("No valid alignment file provided")
  }
  
  map <- header2map(headers=names(seq), aln_path=aln)
  names(seq) <- names(map[["map"]])
  metadata(seq) <- map
  seq
}


header2map <- function (headers, aln_path = aln) {
  xs <- str_split_fixed(headers, ' ', 2)
  genomes <- xs[,1]
  gapranges <- get_alignment_gaps(aln_path, genomes)
  
  map <- strsplit(xs[,2], ' ')
  map <- Map(function (m) do.call(rbind, strsplit(m, ":")), m=map)
  map <- lapply(map, as.data.frame.matrix, stringsAsFactors=FALSE)
  map <- lapply(map, `colnames<-`, value=c("chr", "genomic_start",
                                           "genomic_end", "strand"))
  for (i in seq_along(map)) {
    map[[i]][,c("genomic_start")] <- as.integer(map[[i]][,c("genomic_start")])
    map[[i]][,c("genomic_end")] <- as.integer(map[[i]][,c("genomic_end")])
  }
  
  width <- lapply(map, function (m) {
    m[["genomic_end"]] - m[["genomic_start"]] + 1L
  })
  map <- Map(update_position, map, width, gapranges)
  map <- setNames(map, genomes)
  
  list(map = map, gaps = gapranges) 
}


update_position <- function (m, w, r) {
  s <- 1
  mat <- matrix(rep(0, nrow(m)*2), ncol=2,
                dimnames=list(NULL, c("aln_start", "aln_end")))
  for (i in seq_len(nrow(m))) {
    mat[i,1] <- s
    mat[i,2] <- w[i] + s - 1L
    s <- mat[i,2] + 1L
  }
#   mat <- t(apply(mat, 1, function(x) {
#     make_gapped_range(x[1], x[2], r)
#   }))
  
  cbind(m, mat)
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
  exec <- system.file("src", "runlvec.pl", package="genoslideR")
  #exec <- "~/R/Projects/Devel/genoslideR/inst/src/runlvec.pl"
  gap_ranges <- list()
  for (i in seq_along(genomes)) {
    pipe_desc <- pipe(paste(exec, aln_path, i))
    l <- scan(pipe_desc, sep="\t", quiet=TRUE,
              what=list(pos = integer(), len = integer()))
    close(pipe_desc)
    gap_ranges[[i]] <- IRanges(start=l[["pos"]], width=l[["len"]])
  }
  names(gap_ranges) <- genomes
  gap_ranges  
}

