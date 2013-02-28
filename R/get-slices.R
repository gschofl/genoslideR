#' @importFrom IRanges which
#' @importFrom IRanges order
#' @importFrom IRanges compact
#' @importFrom IRanges subseq
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings xscat
NULL

get_slices <- function (aln, slices) {
  ranges <- ranges(slices)
  genome <- as.character(runValue(seqnames(slices)))
  ## make sure we extract slice with respect to only one genome
  if (length(genome) > 1) {
    stop("More than one genome provided for slicing", call.=FALSE)
  }
  
  # Map of alignment ranges
  g_amap <- aMap(aln)[[genome]]
  # Map of genomic ranges
  g_gmap <- gMap(aln)[[genome]]
  # Map of gaps
  g_gaps <- genoslideR::gaps(aln)[[genome]]

  # pb <- txtProgressBar(max=length(ranges))
  aln_cuts <- vector("list", length(ranges))
  names(aln_cuts) <- names(ranges)
  for (i in seq_along(ranges)) {
    hit <- subjectHits(findOverlaps(ranges[i], ranges(g_gmap), type="any"))
    if (length(hit) == 0) {
      aln_cuts[i] <- DNAStringSet()
    } else {
      hit_range <- g_gmap[hit,]
      hit_order <- IRanges::order(hit_range)
      hit <- hit[hit_order]
      hit_range <- hit_range[hit_order,]
      
      ### generate a cut_range out of a hit_range
      cut_range <- hit_range
      if (start(ranges[i])[1] > start(cut_range)[1]) {
        start(cut_range)[1] <- start(ranges[i])[1]
      }
      if (end(ranges[i])[1] < end(cut_range)[length(cut_range)]) {
        end(cut_range)[length(cut_range)] <- end(ranges[i])[1]
      }
      
      cut_range <- update_reverse_pos(hit_range, cut_range)        
      updated_g_range <- update_genome_pos(cut_range, hit_range, g_amap[hit,])
    
      # cut range with gaps 
      cr <- gap_range(start(updated_g_range), end(updated_g_range), g_gaps, g_gmap)
      aln_cuts[i] <- cut_alignment(cr, updated_g_range, aln)
    }
    # print(i)
    # setTxtProgressBar(pb, i)
  }
  # close(pb)
  aln_cuts
}


update_genome_pos <- function(cut_range, hit_range, pwog_hit_range) {
  start <- start(cut_range) - start(hit_range) + start(pwog_hit_range)
  end <- end(cut_range) - start(hit_range) + start(pwog_hit_range)
  GRanges(seqnames(cut_range), IRanges(start = start, end = end),
          strand = GenomicRanges::strand(cut_range))
}


update_reverse_pos <- function(hit_range, cut_range) {
  idx <- IRanges::which(strand(cut_range) == "-")
  if (all_empty(idx)) {
    return(cut_range)
  }
  
  new_start <- start(cut_range)
  new_end <- end(cut_range)
  new_start[idx] <- end(hit_range)[idx] - (end(cut_range)[idx] - start(hit_range)[idx])
  new_end[idx] <- end(hit_range)[idx] - (start(cut_range)[idx] - start(hit_range)[idx])
  GRanges(seqnames(cut_range), IRanges(new_start, new_end))
}


subseq2 <- function(x, start, end, strand) {
  if (strand == "-") {
    reverse(subseq(x, start, end))
  } else {
    subseq(x, start, end)
  }
} 


cut_alignment <- function (cr, updated_g_range, aln) {
  gmap <- ranges(gMap(aln))
  amap <- ranges(aMap(aln))
  gaps <- ranges(genoslideR::gaps(aln))
  strand <- as.character(strand(updated_g_range))
  seqname <- as.character(runValue(seqnames(updated_g_range)))
  cuts <- mapply(subseq2, start = start(cr), end = end(cr),
                 strand = strand, MoreArgs=list(x = alignment(aln)))
  if (length(cuts) == 1) {
    cuts <- IRanges::compact(cuts[[1]])
  } else {
    cuts <- IRanges::compact(setNames(do.call(xscat, cuts), nm=names(cuts[[1]])))
  }
  
  alignment_position <- GRanges(seqname, cr, strand)
  genomic_position <- ungap_range(cr, amap, gmap, gaps)
  metadata(cuts) <- list(alignment_position = alignment_position,
                         genomic_position = genomic_position)
  cuts
}

