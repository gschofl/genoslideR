#' @importFrom IRanges IRanges
#' @importFrom IRanges IRangesList
#' @importFrom IRanges subjectHits
#' @importFrom IRanges findOverlaps
#' @importFrom GenomicRanges unlist
NULL

make_ungapped_range <- function (start, end, gaprange) {
  ungapped_start <- make_ungapped_position(start, gaprange, right = TRUE)
  ungapped_end <- make_ungapped_position(end, gaprange, right = FALSE)
  if (ungapped_start - 1 > ungapped_end)
    ungapped_start <- ungapped_end <- integer()
  
  IRanges(ungapped_start, ungapped_end)
}


make_ungapped_position <- function(pos, gaprange, right = FALSE) {
  gapstart <- GenomicRanges::start(gaprange)
  gapend <- GenomicRanges::end(gaprange)
  gapwidth <- GenomicRanges::width(gaprange)
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
  gapstart <- GenomicRanges::start(gaprange)
  gapwidth <- GenomicRanges::width(gaprange)
  IRanges(make_gapped_position(start, gapstart, gapwidth),
          make_gapped_position(end, gapstart, gapwidth))
}


cppFunction('
int make_gapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth) {
  int i = 0;
  int n = gapstart.size();
  while (pos >= gapstart[i] && i <= n) {
    pos += gapwidth[i];
    i++;
  }
  return pos;
}')


make_ungapped_genomic_pos <- function (cr, aln) {
  amap <- ranges(aMap(aln))
  gmap <- ranges(gMap(aln))
  gaps <- ranges(genoslideR::gaps(aln))
  
  # degap the alignment range cr - ungapped ranges
  ugr <- lapply(gaps, function (gap) {
    IRangesList(Map(make_ungapped_range, 
                    start = start(cr),
                    end = end(cr),
                    gaprange = list(gap)))
  })
  
  GRangesList(Map(function (ug, g) {
    m <- local({u <- GenomicRanges::unlist(ug)
                hit <- subjectHits(findOverlaps(u, amap[[g]], type="any"))
                hrg <- gmap[[g]][hit, ] # hit range genome
                hra <- amap[[g]][hit, ] # hit range alignment
                IRanges(end(hrg) - (end(hra) - start(u)),
                        end(hrg) - (end(hra) - end(u)))
    })
    if (length(m) > 0) {
      GRanges(g, m)
    } else {
      GRanges()
    }
  }, ug = ugr, g = as.list(names(ugr))))
}

