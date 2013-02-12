
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
  if(not_empty(zero_width)){
    inter <- inter[-zero_width]
    internames <- internames[-zero_width]
  }
  
  ### delete intergenic regions of unknown locus_tag
  which_na <- which(grepl("NA",internames, ignore.case=T, perl=T))
  if(not_empty(which_na)){
    inter <- inter[-which_na]
    internames <- internames[-which_na]
  }
 
  
  names(inter) <- internames
  inter
}
