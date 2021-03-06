# taken from IRanges
labeledLine <- function(label, els, count = TRUE, labelSep = ":", sep = " ", ellipsis = "...") {
    if (count) {
      label <- paste0(label, "(", length(els), ")")
    }
    label <- paste0(label, labelSep, sep)
    width <- getOption("width") - nchar(label)
    line <- ellipsize(els, width, sep, ellipsis)
    paste0(label, line, "\n")
  }

ellipsize <- function(obj, width = getOption("width"), sep = " ", ellipsis = "...") {
    if (length(obj) > 2 * width) {
      obj <- c(head(obj, width), tail(obj, width))
    }
    str <- encodeString(obj)
    ## get order selectSome() would print
    half <- seq_len(ceiling(length(obj) / 2))
    ind <- as.vector(rbind(half, length(obj) - half + 1))
    nc <- cumsum(nchar(str[ind]) + nchar(sep)) - nchar(sep)
    last <- findInterval(width, nc)
    if (length(obj) > last) {
      ## make sure ellipsis fits
      while (last &&
               (nc[last] + nchar(sep)*2^(last>1) + nchar(ellipsis)) > width)
        last <- last - 1L
      if (last == 0) ## have to truncate the first element
        str <-
        paste0(substring(str[1L], 1, width - nchar(ellipsis)), ellipsis)
      else if (last == 1) ## can only show the first
        str <- c(str[1L], "...")
      else
        str <- selectSome(str, last + 1L)
    }
    paste(str, collapse = sep)
  }

## taken directly from Biobase
selectSome <- function (obj, maxToShow = 5) {
  len <- length(obj)
  if (maxToShow < 3) 
    maxToShow <- 3
  if (len > maxToShow) {
    maxToShow <- maxToShow - 1
    bot <- ceiling(maxToShow/2)
    top <- len - (maxToShow - bot - 1)
    nms <- obj[c(1:bot, top:len)]
    c(as.character(nms[1:bot]), "...", as.character(nms[-c(1:bot)]))
  } else {
    obj
  }
}

capitalize <- function(x) {
  substring(x, 1, 1) <- toupper(substring(x, 1, 1))
  x
}

is_mapped_alignment <- function (aln) {
  is(aln, "DNAStringSet") &&
    all(c("aMap", "gMap", "gaps") %in% names(metadata(aln)))
}

is_sliced_alignment <- function (aln) {
  is(aln, "DNAStringSet") &&
    all(c("genomic_position", "alignment_position") %in% names(metadata(aln)))
}

is_fasta <- function(f) {
  l <- lapply(f, readLines, n=2)
  vapply(l, function(l) {
    grepl("^>", l[1]) && 
      all(unique(strsplit(l[2], "")[[1]]) %in% c("A","T","G","C","a","t","g","c"))
  }, logical(1))
}

is_gbk <- function(f) {
  l <- lapply(f, readLines, n=2)
  vapply(l, function(l) {
    grepl("^LOCUS", l[1]) && grepl("^DEFINITION", l[2])
  }, logical(1))
}

is_mfa <- function (f) {
  l <- lapply(f, readLines, n=1)
  vapply(l, function(l) {
    grepl("^>[[:graph:]]+ [[:graph:]]+:\\d+:\\d+:[+-]{1}", l)
  }, logical(1))
}

is_maf <- function (f) {
  l <- lapply(f, function(f) {
    i <- 1
    while (grepl("^#", readLines(f, i))) {
      i <- i + 1
    }
    readLines(f, i + 1)
  })
  vapply(l, function(l) {
    grepl("^a", l[i]) && grepl("^s ", l[i+1])
  }, logical(1))
}

#' @importFrom assertthat see_if
has_dependencies <- function(cmd) {
  lcmd <- vapply(cmd, function (e) see_if(has_command(e)), logical(1))
  if (any(!lcmd)) {
    stop("The following external program(s) must be installed: ", 
         paste(cmd[!lcmd], collapse=", "),
         "\n\nTry running 'install_genoslider_dependencies()''",
         call. = FALSE)
  }
  invisible()
}



