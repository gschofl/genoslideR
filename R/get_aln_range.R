# ncRNAtest <- list(IRanges(start=c(2, 1), end=c(8, 93),
#                           name=c("ncRNA1", "ncRNA2")),
#                   IRanges(start=c(11, 1), end=c(18, 97), 
#                           name=c("ncRNA1", "ncRNA2")),
#                   IRanges(start=c(40, 1), end=c(80, 112), 
#                           name=c("ncRNA1", "ncRNA2")))
# ncRNAtest <- setNames(ncRNAtest, c("NC_002179", "NC_003361", "NC_015408"))

get_aln_range <- function (aln, cut_range) {
  
  if (!is(cut_range, "GRanges")) {
    stop("Provide a 'GRanges' object")
  }
  
  ranges <- ranges(cut_range)
  genome <- seqlevels(cut_range)
  r <- getGaps(aln)[[genome]]
  m <- getMap(aln)[[genome]]
  
  # Map of genomic ranges
  mgr <- GRanges(m[["chr"]], IRanges(m[["genomic_start"]],
                                     m[["genomic_end"]]),
                 strand=m[["strand"]])

  # Map of alignment ranges
  mar <- GRanges(m[["chr"]], IRanges(m[["aln_start"]],
                                     m[["aln_end"]]))
  
  aln_cuts <- vector("list", length(ranges))
  names(aln_cuts) <- names(ranges)
  for (i in seq_along(ranges)) {
    hit <- subjectHits(findOverlaps(ranges[i], ranges(mgr), type="any"))
    if (length(hit) == 0) {
      aln_cuts[i] <- BStringSet()
    } else {
      hit_range <- mgr[hit,]
      hit_order <- order(ranges(hit_range))
      hit <- hit[hit_order]
      hit_range <- hit_range[hit_order,]
      
      ### generate new_x out of hit_range and x[i] range
      new_x <- hit_range
      if (start(ranges[i])[1] > start(new_x)[1]) {
        start(new_x)[1] <- start(ranges[i])[1]
      }
      
      if (end(ranges[i])[1] < end(new_x)[length(new_x)]) {
        end(new_x)[length(new_x)] <- end(ranges[i])[1]
      }
      
      new_x <- update_reverse_pos(hit_range, cut_range=new_x)        
      updated_g_range <- update_genome_pos(new_x, hit_range, mar[hit,])
      
      # cut range with gaps
      crl <- Map(make_gapped_range, start(updated_g_range),
                 end(updated_g_range), list(r))
      cr <- IRanges(vapply(crl, `[`, 1, FUN.VALUE=integer(1)),
                    vapply(crl, `[`, 2, FUN.VALUE=integer(1)))
      
      aln_cuts[i] <- cut_alignment(aln, cr, updated_g_range)
    }
  }
  aln_cuts
}


update_genome_pos <- function(cut_range, hit_range, pwog_hit_range) {
  
  start <- start(cut_range) - start(hit_range) + start(pwog_hit_range)
  end <- end(cut_range) - start(hit_range) + start(pwog_hit_range)
  
  GRanges(seqlevels(cut_range), IRanges(start = start, end = end),
          strand = GenomicRanges::strand(cut_range))
}


update_reverse_pos <- function(hit_range, cut_range) {
  idx <- which(GenomicRanges::strand(cut_range) == "-")
  if (length(idx) == 0) {
    return(cut_range)
  }
  
  new_start <- start(cut_range)
  new_end <- end(cut_range)
  new_start[idx] <- end(hit_range)[idx] - (end(cut_range)[idx] -
                                             start(hit_range)[idx])
  new_end[idx] <- end(hit_range)[idx] - (start(cut_range)[idx] -
                                           start(hit_range)[idx])
  GRanges(seqnames(cut_range), IRanges(new_start, new_end))
}


cut_alignment <- function (aln, cr, updated_g_range) {
  
  strand <- as.character(GenomicRanges::strand(updated_g_range))
  seqlevel <- seqlevels(updated_g_range)
  cuts <- vector("list", length(strand))
  
  for (i in seq_along(cuts)) {
    if (strand[i] == "-") {
      cuts[[i]] <- reverse(subseq(alignment(aln), start(cr)[i], end(cr)[i]))
      metadata(cuts[[i]]) <- list(alignment_position = GRanges(seqlevel, cr[i], strand = strand[i]))
    } else {
      cuts[[i]] <- subseq(alignment(aln), start(cr)[i], end(cr)[i])
      metadata(cuts[[i]]) <- list(alignment_position = GRanges(seqlevel, cr[i], strand = strand[i]))
    }
  }

  if (length(cuts) == 1) {
    metadata(cuts[[1]]) <- list(alignment_position = metadata(cuts[[1]])[["alignment_position"]],
                                genomic_position = make_ungapped_genomic_pos(cr, map=getMap(aln), gaps=getGaps(aln)))
    return(unlist(cuts))
  } else {
    md <- lapply(cuts, function (x) metadata(x)[["alignment_position"]])
    md <- do.call(c, md)
    cuts <- do.call(xscat, cuts)
    metadata(cuts) <- list(alignment_position = md,
                           genomic_position = make_ungapped_genomic_pos(cr, map=getMap(aln), gaps=getGaps(aln)))
    return(cuts)
  }
}

