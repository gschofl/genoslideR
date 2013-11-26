count_lines <- function(file) {
  count <- vapply(file, function(f) {
    wc <- SysCall("wc", l=TRUE, stdin=f, redirection=FALSE, intern=TRUE)
    as.integer(strsplit(wc, ' ')[[1L]][1L])
  }, FUN.VALUE=integer(1))
  names(count) <- NULL
  count
}

#' Replace CDS indices in an oframe by the corresponding gbFeatureList indeces.
#' 
#' @param oframe
#' @param flist A (list of) \code{\link[biofiles]{gbFeatureList}}s
#' @importFrom biofiles select index start getAccession
#' @export
substituteCDSNumbers <- function(oframe, flist) {
  if (!is(oframe, "OrthoFrame")) {
    stop("oframe needs to be of class OrthoFrame")
  }
  if (!is(flist, "gbRecordList")) {
    stop("flists need to be of class gbRecordList")
  }
  if (!all(biofiles::getAccession(flist) %in% labels(oframe))) {
    stop("labels of oframe and accession numbers in flist do not match")
  }
  accn <- labels(oframe)
  names(oframe) <- accn
  df <- oframe
  order_names <- names(oframe)
  anch <- anchors(oframe)
  for (i in seq_len(ncol(oframe))) {
    f <- flist[[ which(biofiles::getAccession(flist)%in%accn[i]) ]]
    accession <- biofiles::getAccession(f)
    cds <- select(f, key="CDS") 
    idx <- index(cds)    
    s1 <- start(cds) 
    s2 <- anch[[i]][, 4] + 1
    idx <- idx[which(s1 %in% s2)]
    idx <- data.frame(cbind(seq_along(idx), idx))
    names(idx) <- c(accession, paste0(accession, ".idx"))
    df <- merge(x=df, y=idx, all.x=TRUE)
  }
  
  df <- df[, order(names(df))]
  df <- df[, grepl("idx$", names(df))]
  df <- df[ do.call(order, df), ]
  names(df) <- strip_ext(names(df))
  df <- df[, order_names]
  class(df) <- class(oframe)
  labels(df) <- labels(oframe)
  attr(df, "anchors") <- attr(oframe, "anchors")
  rownames(df) <- NULL
  df
}

