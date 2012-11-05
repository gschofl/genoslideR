make_ungapped_range <- function (start, end, gaprange) {
  ungapped_start <- make_ungapped_position(start, gaprange, right = TRUE)
  ungapped_end <- make_ungapped_position(end, gaprange, right =  FALSE)
  if (ungapped_start - 1 > ungapped_end)
    ungapped_start <- ungapped_end <- integer()
  
  c(ungapped_start, ungapped_end)
}


make_ungapped_position <- function(pos, gaprange, right = FALSE) {
  gapstart <- IRanges::start(gaprange)
  gapend <- IRanges::end(gaprange)
  gapwidth <- IRanges::width(gaprange)
  preceding <- which(gapend < pos)
  if (length(overlapping <- setdiff(which(gapstart <= pos), preceding)) > 0) {
    if (right) {
      pos - sum(gapwidth[preceding]) - (pos - gapstart[overlapping])
    } else {
      pos - sum(gapwidth[preceding]) - (pos - gapstart[overlapping] + 1)
    }  
  } else {
    pos - sum(gapwidth[preceding]) 
  }
}


make_gapped_range <- function (start, end, gaprange) {
  c(make_gapped_position(start, gaprange), make_gapped_position(end, gaprange))
}


make_gapped_position <- function(pos, gaprange) {
  i <- 1
  gapstart <- IRanges::start(gaprange)
  gapwidth <- IRanges::width(gaprange)
  while (pos >= gapstart[i] && i <= length(gapstart)) {
    pos <- pos + gapwidth[i]
    i <- i + 1
  }
  pos
}


make_ungapped_genomic_pos <- function (cr, map, gaps) {
  
  mar <- IRangesList(lapply(map, function (m) {
    IRanges(m[["aln_start"]], m[["aln_end"]])
  }))
  
  mgr <- IRangesList(lapply(map, function (m) {
    IRanges(m[["genomic_start"]], m[["genomic_end"]])
  }))
  
  ugr <- lapply(gaps, function (gap) {
    Map(function(s, e, g) {
      r <- make_ungapped_range(s, e, g)
      IRanges(r[1], r[2])
    }, s = start(cr), e = end(cr), g = list(gap))
  })
  
  hits <- Map(function (ug, genomes) {
    m <- Map(function (u, g) {
      hit <- subjectHits(findOverlaps(u, mar[[g]], type="any"))
      hrg = mgr[[g]][hit, ] # hit range genome
      hra = mar[[g]][hit, ] # hit range alignment
      IRanges(end(hrg) - (end(hra) - start(u)),
              end(hrg) - (end(hra) - end(u)))
    }, u = ug, g = genomes)
    
    if (!all(vapply(m, length, numeric(1)) == 0)) {
      GRanges(genomes, do.call(c, m))
    } else {
      GRanges()
    }
  }, ug = ugr, genomes = names(ugr))
  
  GRangesList(hits)
}

