interval <- function(x, y) {
  if (!is.numeric(c(x, y)))
    stop("Argument is not numeric")
  if (length(x) != length(y))
    stop("Arguments of unequal lengths")
  if (any(swap <- !(x <= y))) {
    tmp_x <- y[swap]
    y[swap] <- x[swap]
    x[swap] <- tmp_x
  }
  structure(matrix(c(x, y, y - x + 1), nrow=length(x), ncol=3,
                   dimnames=list(NULL, c("start", "end", "width"))),
            class=c("interval", "matrix"))
}

intervals_overlap <- function(int1, int2=NULL) {
  if (!is(int1, "interval") || missing(int1))
    stop("'int1' is not of class 'interval'")
  
  if (is(int2, "interval") &&
        (nrow(int1) == nrow(int2) || nrow(int1) == 1L))
    return(int1[,"start"] <= int2[,"end"] & int1[,"end"] >= int2[,"start"])
  
  if (is.null(int2)) {
    if (nrow(int1) < 2L)
      return(FALSE)
    else {
      combi <- combn(nrow(int1), 2)
      ol <- int1[combi[1,],"start"] <= int1[combi[2,],"end"] &
        int1[combi[1,],"end"] >= int1[combi[2,],"start"]
      fr <- data.frame(t(combi), ol)
      colnames(fr) <- c("int1", "int2", "overlap")
      fr
    }
  }    
}

combine_hits <- function(hit) {
  q_int <- interval(hit$qstart, hit$qend)
  q_overlap <- intervals_overlap(q_int)
  s_int <- interval(hit$sstart, hit$send)
  s_overlap <- intervals_overlap(s_int)
  if (any(q_overlap$overlap) || any(s_overlap$overlap)) {
    qol <- q_overlap$overlap
    sol <- s_overlap$overlap
    # check if any of the query intervals or subject intervals overlap
    # and get the combinations of intervals that do overlap (=ol_int)
    ol_int <- unique(c(unlist(s_overlap[sol, c(1,2)]), unlist(q_overlap[qol, c(1,2)])))
    # delete those intervals that overlap except for the largest
    ol_int_del <- ol_int[-which.max(q_int[ol_int,"width"])]
    hit <- hit[-ol_int_del,]
  }
  # combine hits
  qid <- unique(hit$qid)
  sid <- unique(hit$sid)
  pident <- (sum(hit$length) - sum(hit$mismatch))/sum(hit$length)
  length <- sum(hit$length)
  mismatch <- sum(hit$mismatch)
  gapopen <- sum(hit$gapopen)
  qstart <- min(hit$qstart)
  qend <- max(hit$qend)
  sstart <- min(hit$sstart)
  send <- max(hit$send)
  bitscore <- sum(hit$bitscore)
  evalue <- prod(hit$evalue)
  data.frame(
    stringsAsFactors=FALSE,
    qid=qid, sid=sid, pident=pident, length=length, mismatch=mismatch,
    gapopen=gapopen, qstart=qstart, qend=qend, sstart=sstart, send=send,
    evalue=evalue, bit_score=bitscore
  )
}

blat2hits <- function(blat_output) {
  stopifnot(require(blastr))
  blat <- blastr::blastTable(blat_output)[[1]]
  groups <- as.factor(paste0(blat@table$qid, "_", blat@table$sid))
  split_blat <- split.data.frame(blat@table, groups)
  if (any(mult_hits <- vapply(split_blat, nrow, 0L, USE.NAMES=FALSE) > 1L)) {     
    combined_hits <- lapply(split_blat[mult_hits], combine_hits)
    split_blat[mult_hits] <- combined_hits
  }
  blat <- do.call("rbind", split_blat)
  hits <- data.frame(
    qid=blat$qid, sid=blat$sid, pident=blat$pident,
    bit_score=blat$bit_score, evalue=ifelse(blat$evalue < 1e-50, 0e+00, blat$evalue)
  ) 
  hitfile <- replace_ext(blat_output, replacement="hits", level=1)
  write.table(x=hits, file=hitfile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE,
              col.names=FALSE)
  return(invisible(hits))
}
