
#chps_gff <- import_annotation_from_gff("~/daten/alignment/Chlamydophila_psittaci_02DC15_uid159521/NC_017292.gff")

#' Generate intergenic regions
#' 
#' @param gff A gff annotation file to calculate the intergenic regions.
#'
#' @details Zero length ranges and unknown locus_tags are excluded.
#' 
#' @return An IRanges object including the start and end plus the intergenic 
#' name of each intergenic region.
#' @export

get_intergenic_ranges <- function(gff=chps_gff){
  
  ### extract information of gff file
  locus_tag <- gff$synonym
  genome_size <- gff@seqinfo@seqlengths
  gene_ranges <- ranges(gff)
  
  ### generate intergenic ranges
  gene_start <- c(start(gene_ranges)-1, genome_size)
  gene_end <- c(1, end(gene_ranges))
  
  t <- Map(function(s, e) {
    if(e < s){
      IRanges(1, 0)   
    }else{
      IRanges(s, e)
    }
  }, s = gene_end, e = gene_start)
  inter <- do.call(c, t)
  
  ### generate intergenic names
  internames <- do.call(c, strsplit(paste(locus_tag, locus_tag, sep=":"), ":"))  
  if(start(gene_ranges[1]) != 1){
    internames <- c("begin", internames)
  }
  if(end(gene_ranges[length(gene_ranges)]) != genome_size){
    internames <- c(internames, "end")
  }
  internames <- paste(internames[seq(1,length(internames),2)], internames[seq(2,length(internames),2)], sep=":")

  ### delete intergenic regions of zero length
  zero_width <- which(width(inter) == 0)
  if(!all_empty(zero_width)){
    inter <- inter[-zero_width]
    internames <- internames[-zero_width]
  }
  
  ### delete intergenic regions of unknown locus_tag
  which_na <- which(grepl("NA",internames, ignore.case=T, perl=T))
  if(!all_empty(which_na)){
    inter <- inter[-which_na]
    internames <- internames[-which_na]
  }
 
  
  names(inter) <- internames
  inter
}

### strand option not implemented ###
cutAlignment <- function (aln = chlam, start, end, names = NULL, 
                          strand = NULL, genome = 2) {
  
  if (length(genome) > 1) {
    warning("More than one genome specified. Only the first will be used")
    genome <- genome[1]
  }
  
  if (is.numeric(genome)) {
    genome <- names(aln)$alignment[genome]
    if (is.na(genome)) {
      stop("No annotation for genome ", which(is.na(genome))) 
    }
  } else {
    if (!genome %in% seqlevels(aln)) {
      stop("Genome ", sQuote(genome), " not present in annotation.")
    }
  }
  
  #   if (!is.null(names)) {
  #     names <- ifelse(is.na(names), "", names)
  #   }
  cores <- detectCores()
  
  res <- mcmapply(function (s, e, n) {
    cut_range <- GRanges(genome, IRanges(start=s, end=e, names=n))
    get_aln_range(aln, cut_range)
  }, s=start, e=end, n=names, mc.cores = cores - 1)
  
  
  #   cut_range <- GRanges(genome, IRanges(start=start, end=end, names=names))
  #   res <- get_aln_range(aln, cut_range)
  
  if (!is.null(strand)) {
    idx <- which(strand == "-")
    res[idx] <- lapply(res[idx], function (s) {
      reverseComplement(as(s, "DNAStringSet"))
    })
  }
  
  res
}
