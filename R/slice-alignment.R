slice_alignment <- function (segment_dir, write = TRUE) {

  parent_dir <- normalizePath(strsplitN(segment_dir, "segments", 1))
  segments <- dir(segment_dir, "^\\d+$", full.names=TRUE)
  map <- read.table(file.path(segment_dir, "map"),
                    as.is=TRUE, sep="\t", row.names=1)
  genomes <- scan(file.path(segment_dir, "genomes"), what="character",
                  quiet=TRUE)
  
  # order segments
  segments <- segments[order(as.numeric(split_path(segments)))]
  
  # parse map into fasta headers
  headers <- map2header(map, genomes)
  
  # merge segments into one multifasta
  merged_aln <- merge_segments(segments, headers)
  
  if (write) {
    out <- file.path(parent_dir, "aln.mfa")
    writeXStringSet(merged_aln, filepath=out, format="fasta")
    metadata(merged_aln) <- list(path=out)
  }

  return(invisible(merged_aln))
}


map2header <- function (map, genomes) {
  idx <- Map(base::seq, seq(1, by=4, length.out=length(genomes)),
             seq(4, by=4, length.out=length(genomes)))
  m <- do.call("rbind", lapply(idx, function (j) {
    m <- map[,j]
    colnames(m) <- c("chr", "start", "end", "strand")
    m
  }))
  m <- m[complete.cases(m), ]
  m[["start"]] <- m[["start"]] + 1
  m <- split(m, as.factor(m$chr))
  
  Map(collapse_slices, df=m, name=genomes)
} 


collapse_slices <- function (df, name) {
  slices <- gsub("\\s+", "", (apply(df, 1, paste, collapse=":")))
  paste(name, paste(slices, collapse=" "))
}


merge_segments <- function(segments, headers) {
  genomes <- names(headers)
  paths <- file.path(segments, "mavid.mfa")
  seqlist <- lapply(paths, readDNAStringSet)
  mat <- t(vapply(Map(base::names, seqlist), function (n) genomes%in%n,
                  logical(length(genomes))))
  colnames(mat) <- genomes
  for (i in seq_along(seqlist)) {
    if (!all(mat[i,])) {
      present <- which(mat[i,])
      absent <- which(mat[i,] == FALSE)
      w <- unique(Biostrings::width(seqlist[[i]]))
      gaps <- dup('-', w)
      missing_seqs <- setNames(DNAStringSet(rep(gaps, length(absent))),
                               nm=names(absent))
      seqlist[[i]] <- 
        c(seqlist[[i]], missing_seqs)[order(c(present, absent))]
    }
  }
  setNames(do.call(xscat, seqlist), unlist(headers))
}


