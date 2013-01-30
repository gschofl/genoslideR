#' Import aligned mercator segments or a multifasta alignment
#' 
#' @param alnPath Path to mercator segments or a multifasta segment
#' @param writeToMFA if \code{TRUE} merge mercator segments into a
#' single multi fasta file and write to file "aln.mfa"  
#'
#' @export
importAlignment <- function (alnPath, writeToMFA = TRUE) {
  
  if (is_segments_dir(alnPath)) {
    seq <- merge_mercator_segments(segments_dir(alnPath), write=writeToMFA)
    alnPath <- metadata(seq)[["path"]]
  } else {
    seq <- readDNAStringSet(alnPath)
  }
  
  map <- header2map(names(seq), alnPath)
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


get_alignment_gaps <- function (aln_path, genomes) {
  exec <- system.file("src", "runlvec.pl", package="genoslideR")
  #exec <- "./src/runlvec.pl"
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

