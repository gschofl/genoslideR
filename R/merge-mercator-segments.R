#' @importFrom Biostrings writeXStringSet readDNAStringSet DNAStringSet xscat
NULL

merge_mercator_segments <- function (seg_dir, outfile =  NULL) {
  parent_dir <- normalizePath(strsplitN(seg_dir, "segments", 1))
  segments <- dir(seg_dir, "^\\d+$", full.names=TRUE)
  segments <- segments[order(as.numeric(split_path(segments)))]
  map <- read.table(file.path(seg_dir, "map"), as.is=TRUE, sep="\t", row.names=1)
  genomes <- scan(file.path(seg_dir, "genomes"), what="character", quiet=TRUE)
  
  # parse map into fasta headers
  headers <- map2header(map, genomes)
  
  # merge segments into one multi fasta file
  aln <- merge_segments(segments, headers)
  
  if (!is.null(outfile)) {
    aln_path <- outfile %||% file.path(seg_dir, "fsa.mfa")
    writeXStringSet(aln, aln_path)
  }
  
  return(invisible(aln))
}


map2header <- function (map, genomes) {
  len <- length(genomes) 
  by <- gl(len, nrow(map))
  idx <- Map(base::seq, seq(from=1, by=4, length.out=len),
             seq(from=4, by=4, length.out=len))
  m <- do.call("rbind", lapply(idx, function (j) {
    m <- map[,j]
    colnames(m) <- c("chr", "start", "end", "strand")
    m
  }))
  m[["start"]] <- m[["start"]] + 1L
  m <- setNames(split(m, by), nm=genomes)
  coverage <- Map(function(x) max(x[["end"]], na.rm=TRUE), m)
  Map(collapse_slices, df=m, name=genomes, cov=coverage)
} 


collapse_slices <- function (df, name, cov) {
  slices <- apply(df, 1, function (x) paste0(rmisc::trim(x), collapse=":"))
  paste(name, cov, paste(slices, collapse=" "))
}


merge_segments <- function(segments, headers) {
  genomes <- names(headers)
  paths <- file.path(segments, "fsa.mfa")
  seqlist <- lapply(paths, readDNAStringSet)
  
  empty_seqment <- vapply(seqlist, length, numeric(1)) == 0L
  if (any(empty_seqment)) {
    stop("No alignment in segment: ", paths[empty_seqment], call.=FALSE)
  }
  
  mat <- t(vapply(Map(names, seqlist), function (n) genomes%in%n,
                  logical(length(genomes))))
  colnames(mat) <- genomes
  incomplete <- which(!apply(mat, 1, all))
  for (i in incomplete) {
    present <- which(mat[i,] == TRUE)
    absent <- which(mat[i,] == FALSE)
    w <- unique(Biostrings::width(seqlist[[i]]))
    gaps <- dup('-', w)
    missing_seqs <- setNames(DNAStringSet(rep(gaps, length(absent))),
                             nm=names(absent))
    seqlist[[i]] <- 
      c(seqlist[[i]], missing_seqs)[order(c(present, absent))]
  }
  setNames(do.call(xscat, seqlist), unlist(headers))
}


