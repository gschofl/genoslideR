#' Generate summaries on sliding windows over an alignment.
#' 
#' @details
#' The function \code{f} applied to each window must accept a 
#' \code{\linkS4class{DNAStringSet}} as input and generate a named
#' vector as output.
#'
#' @param aln An \code{\linkS4class{annotatedAlignment}} or
#' \code{\linkS4class{DNAStringSet}} object.
#' @param ranges A \code{\linkS4class{GRanges}} object holding genomic windows
#' on which summaries are generated.
#' @param window Sliding window size. If no \code{ranges} are specified sliding
#' windows are constructed over the complete alignment, otherwise sliding
#' windows are constructed within each genomic window specified by \code{ranges}.
#' @param step Step size of sliding windows.
#' @param f Function applied to each window of \code{aln}.
#' @param ... Further arguments for \code{f}.
#' @return A \code{\linkS4class{RangedDataList}}.
#' @export
aggregateAlignment <- function(aln, ranges = NULL, window = NULL, step = NULL,
                               f = "Eta", ...) {
  
  if (is(aln, "annotatedAlignment")) {
    dss <- alignment(aln)
  } else if (is(aln, "DNAStringSet")) {
    dss <- aln
  } else {
    stop("Only DNAStringSet supported")
  }
  
  f <- match.fun(f)
  if (is.null(ranges)) {
    ranges <- slidingWindows(1L, unique(width(dss)), window, step)
    ans <- aggregateRanges(dss, ranges, f)
  } else {
    grl <- genome2Alignment(ranges, aln)
    grl <- mapply(slidingWindows, start = start(grl), width = width(grl),
                  MoreArgs=list(window = window, step = step))
    ans <- lapply(grl, aggregateRanges, dss = dss, f = f, ... = ...)
  }
  RangedDataList(ans)
}


aggregateRanges <- function(dss, ranges, f, ...) {
  
  space <- NULL
  if (is(ranges, "IRangesList")) {
    p <- ranges@partitioning
    space <- factor(unlist(Map(rep, seq_along(p), times=width(p))))
    ranges <- ranges@unlistData
  }
  
  ans <- lapply(seq_along(ranges), function(i) {
    f(narrowAlignment(dss, start=ranges@start[i], width=ranges@width[i]), ...)
  })
  
  nmr <- names(ranges)
  if (any(duplicated(nmr))) {
    dups <- duplicated(nmr)
    x <- rle(as.integer(dups))
    suff <- unlist(Map(seq_len, x$lengths[x$values == 1]))
    names(ranges)[dups] <- paste0(nmr[dups], '.', suff)
  }
  
  rd <- RangedData(ranges, DataFrame(do.call(rbind, ans)), space=space)
  rd
}


slidingWindows <- function(start = 1L, width, window = NULL, step = NULL) {
  if (length(start) == 1L)
    return( .sw(start, width, window, step) )
  
  IRangesList(mapply(.sw, start = start, width = width,
         MoreArgs=list(window = window, step = step)))
}


.sw <- function(start, width, window = NULL, step = NULL) {
  if (is.null(window) && is.null(step))
    r <- IRanges::new2("IRanges", start = start, width = width)
  else {
    nSteps <- (width - window) %/% step + 1
    stepWidth <- step * nSteps
    lastWindowWidth <- width - stepWidth
    startvec <- seq.int(start, by=step, length.out=nSteps + 1)
    widthvec <- as.integer(c(rep(window, nSteps), lastWindowWidth))
    r <- IRanges::new2("IRanges", start = startvec, width = widthvec)
  }
  r
}

