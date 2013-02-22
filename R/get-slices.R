#' @importFrom IRanges which
#' @importFrom IRanges order
#' @importFrom IRanges compact
#' @importFrom Biostrings DNAStringSet
NULL

get_slices <- function (aln, slices) {
  ranges <- ranges(slices)
  genome <- as.character(runValue(seqnames(slices)))
  
  ## make sure we extract slice with respect to only one genome
  if (length(genome) > 1) {
    stop("More than one genome provided for slicing", call.=FALSE)
  }
  
  # Map of gaps
  mgps <- genoslideR::gaps(aln)[[genome]]
  # Map of genomic ranges
  mgr <- gMap(aln)[[genome]]
  # Map of alignment ranges
  mar <- aMap(aln)[[genome]]
  
  pb <- txtProgressBar(max=length(ranges))
  aln_cuts <- vector("list", length(ranges))
  names(aln_cuts) <- names(ranges)
  for (i in seq_along(ranges)) {
    hit <- subjectHits(findOverlaps(ranges[i], ranges(mgr), type="any"))
    if (length(hit) == 0) {
      aln_cuts[i] <- DNAStringSet()
    } else {
      hit_range <- mgr[hit,]
      hit_order <- IRanges::order(ranges(hit_range))
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
      updated_g_range <- update_genome_pos(cut_range, hit_range,
                                           pwog_hit_range=mar[hit,])
      
      # cut range with gaps
      cr <- unlist(IRangesList(Map(make_gapped_range,
                                   start = start(updated_g_range),
                                   end = end(updated_g_range),
                                   gaprange = list(mgps))))
      aln_cuts[i] <- cut_alignment(aln, cr, updated_g_range)
    }
    print(i)
    #setTxtProgressBar(pb, i)
  }
  close(pb)
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


cut_alignment <- function (aln, cr, updated_g_range) {
  
  strand <- as.character(strand(updated_g_range))
  seqname <- as.character(runValue(seqnames(updated_g_range)))
  cuts <- vector("list", length(strand))
  
  for (i in seq_along(cuts)) {
    if (strand[i] == "-") {
      cuts[[i]] <- reverse(subseq(alignment(aln), start(cr)[i], end(cr)[i]))
    } else {
      cuts[[i]] <- subseq(alignment(aln), start(cr)[i], end(cr)[i])
    }
  }
  
  if (length(cuts) == 1) {
    cuts <- IRanges::compact(cuts[[1]])
  } else {
    cuts <- IRanges::compact(setNames(do.call(xscat, cuts), nm=names(cuts[[1]])))
  }
  
  alignment_position <- GRanges(seqname, cr, strand)
  genomic_position <- make_ungapped_genomic_pos(cr, aln)
  metadata(cuts) <- list(alignment_position = alignment_position,
                         genomic_position = genomic_position)
  cuts
}

