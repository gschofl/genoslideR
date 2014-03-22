#' @importFrom plyr ddply summarise
NULL

interval <- function(x, y) {
  assert_that(is.numeric(c(x, y)))
  assert_that(length(x) == length(y))
  o <- order(x)
  x <- x[o]
  y <- y[o]
  if (any(swap <- !(x <= y))) {
    tmp_x <- y[swap]
    y[swap] <- x[swap]
    x[swap] <- tmp_x
  }
  structure(matrix(c(x, y, y - x + 1), nrow=length(x), ncol=3,
                   dimnames=list(NULL, c("start", "end", "width"))),
            class=c("interval", "matrix"))
}

intervals_overlap <- function(int1, int2 = NULL) {
  if (!is(int1, "interval") || missing(int1)) {
    stop("'int1' is not of class 'interval'")
  }
  if (is(int2, "interval") &&
        (nrow(int1) == nrow(int2) || nrow(int1) == 1L)) {
    return(int1[,"start"] <= int2[,"end"] & int1[,"end"] >= int2[,"start"])
  }
  if (is.null(int2)) {
    if (nrow(int1) < 2L) {
      return(FALSE)
    } else {
      combi <- combn(nrow(int1), 2)
      ol <- int1[combi[1,], "start"] <= int1[combi[2,], "end"] &
        int1[combi[1,], "end"] >= int1[combi[2,], "start"]
      fr <- data.frame(t(combi), ol)
      colnames(fr) <- c("int1", "int2", "overlap")
      fr
    }
  }    
}

merge_hsps <- function(hit) {
  q_int <- interval(x = hit$qstart, y = hit$qend)
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
    ol_int_del <- ol_int[-which.max(q_int[ol_int, "width"])]
    hit <- hit[-ol_int_del,]
  }
  # combine hits
  qid <- unique(hit$qid)                    ## query id
  sid <- unique(hit$sid)                    ## subject id
  length <- sum(hit$length)                 ## aligned length: sum of all hsp lengths
  mismatch <- sum(hit$mismatch)             ## cumulative mismatches
  pident <- round((length - mismatch)/length*100, 2)  ## cumulative percentage of identity
  qlen <- unique(hit$qlen)                  ## query length
  paln <- round(length/qlen*100, 2)         ## cumulative percentage of alignment length
  gapopen <- sum(hit$gapopen)
  qstart <- min(hit$qstart)
  qend <- max(hit$qend)
  sstart <- min(hit$sstart)
  send <- max(hit$send)
  bit_score <- sum(hit$bit_score)
  evalue <- prod(hit$evalue)
  data.frame(
    stringsAsFactors = FALSE,
    qid = qid, sid = sid, pident = pident, length = length,
    mismatch = mismatch, gapopen = gapopen, qstart = qstart, qend = qend,
    sstart = sstart, send = send, evalue = evalue, bit_score = bit_score,
    qlen = qlen, paln = paln
  )
}

blat_add_query_length <- function(blat_output) {
  blat <- blastr::blastTable(blat_output)[[1]]
  a <- file.path(dirname(blat_output),
                 paste0(usplit(strip_ext(basename(blat_output)), "-")[1], ".anchors"))
  qid <- blat@table$qid
  qa <- read.table(a)
  qa <- ddply(qa[qa$V1 %in% qid, ], .(V1), summarise, qid = V1, qlen = abs(V5 - V4)/3)[, -1]
  qrl <- rle(as.numeric(qid))
  qlen <- numeric()
  for (i in seq_along(qrl$lengths)) {
    qlen <- c(qlen, rep(qa[qa$qid %in% qrl$values[i], "qlen"], qrl$lengths[i]))
  }
  blat@table$qlen <- qlen
  blat@table$paln <- round((blat@table$qend - blat@table$qstart)/qlen*100, 2)
  blat
}

blat2hits <- function(blat_output) {
  stopifnot(require(blastr))
  blat <- blat_add_query_length(blat_output)
  groups <- as.factor(paste0(blat@table$qid, "_", blat@table$sid))
  split_blat <- split.data.frame(blat@table, groups)
  if (any(mult_hits <- vapply(split_blat, nrow, 0L, USE.NAMES=FALSE) > 1L)) {     
    combined_hits <- lapply(split_blat[mult_hits], merge_hsps)
    split_blat[mult_hits] <- combined_hits
  }
  blat <- do.call("rbind", split_blat)
  hits <- data.frame(
    qid = blat$qid, sid = blat$sid, pident = blat$pident, paln = blat$paln,
    bit_score = blat$bit_score, evalue = ifelse(blat$evalue < 1e-50, 0e+00, blat$evalue)
  ) 
  hitfile <- replace_ext(blat_output, replacement = "hits", level = 1)
  write.table(x = hits, file = hitfile, append = FALSE, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  return(invisible(hits))
}
