#' @importFrom genomeIntervals readGff3 getGffAttribute
NULL

gff2anchors <- function (gff, anchor, id_as_anchor_num=FALSE) {
  gff <- readGff3(file=gff)
  
  getPosition <- function (gff) {
    start <- vector(mode="integer", length=nrow(gff))
    end <- vector(mode="integer", length=nrow(gff))
    pos <- which(strand(gff) == "+")
    neg <- which(strand(gff) == "-")
    start[pos] <- gff[pos,1] - 1 # + gff[pos]$phase
    end[pos] <- start[pos] + ((gff[pos,2] - start[pos]) / 3) * 3
    end[neg] <- gff[neg,2] # - gff[neg]$phase
    start[neg] <- end[neg] - ((end[neg] - (gff[neg,1] - 1)) / 3) * 3
    return(list(start=start, end=end))
  }
  
  cds <- gff[attr(gff, "annotation")$type == "CDS" | attr(gff, "annotation")$type == "pseudogenic_exon",]
  pos <- getPosition(cds)
  anchor.df <- data.frame(anchorNum=if (id_as_anchor_num) {
                            escape(getGffAttribute(cds, attribute="ID"))
                          } else {
                            seq_len(nrow(cds))
                          },
                          seqname=seq_name(cds),
                          strand=strand(cds),
                          start=pos$start,
                          end=pos$end,
                          n=rep(1, nrow(cds)))
  
  write.table(x=anchor.df, file=anchor, sep="\t", quote=FALSE,
              col.names=FALSE, row.names=FALSE)
}

pairwiseCombine <- function (x) {
  first_combin <- utils::combn(x, 2)
  recip_combin <- rbind(first_combin[2,], first_combin[1,])
  cbind(first_combin, recip_combin)
}


escape <- function (s) {
  s <- gsub(" ", "%20", s)
  s
}


unescape <- function (s) {
  s <- gsub("%20", " ", s)
  s
}


reciprocal_blat <- function (genomes, sdb, gff, out, sep = "-",
                             mercator = FALSE) {
  if (length(genomes) < 2)
    stop("At least two genomes must be specified")
  
  ## Dependencies
  has_dependencies(c("sdbList", "gffRemoveOverlaps", "gff2anchors",
                     "anchors2fa", "blat", "blat2hits"))
  
  # Make sure that a gff and sdb file are present for all genomes
  for (genome in genomes) {
    gff_file <- file.path(gff, paste0(genome, ".gff"))
    sdb_file <- file.path(sdb, paste0(genome, ".sdb"))
    
    if (!file.exists(gff_file))
      stop(sprintf("gff file for %s does not exist", sQuote(genome)))
    
    if (!file.exists(sdb_file))
      stop(sprintf("SDB file for %s does not exist", sQuote(genome)))
  }
  
  # Make anchor, chrom, and protein sequence file
  # genome <- genomes[1]
  for (genome in genomes) {
    gff_file <- file.path(gff, paste0(genome, ".gff"))
    sdb_file <- file.path(sdb, paste0(genome, ".sdb"))
    
    # make chrom file
    cat(sprintf("Making chromosome file for %s...\n", sQuote(genome)))
    chrom_file <- file.path(out, paste0(genome, ".chroms"))
    system(sprintf("sdbList -l %s > %s", sdb_file, chrom_file))
    cat("done\n")
    
    # make anchor file
    cat(sprintf("Making anchors for %s...\n", sQuote(genome)))
    anchor_file <- file.path(out, paste0(genome, ".anchors"))
    if (mercator) {
      system(sprintf("cat %s | gffRemoveOverlaps -fCDS -l | gff2anchors > %s",
                     gff_file, anchor_file))
    } else {
      gff2anchors(gff_file, anchor_file)
    }
    cat("done\n")
    
    # Make protein file
    cat(sprintf("Extracting protein sequences from %s...\n", sQuote(genome)))
    prot_file <- file.path(out, paste0(genome, ".proteins.fa"))
    system(sprintf("anchors2fa %s < %s > %s", sdb_file, anchor_file, prot_file))
    cat("done\n")
  }
  
  # BLAT protein anchors pairwise
  cat("BLATing anchors ...\n")
  genome_pairs <- t(pairwiseCombine(genomes))
  for (i in seq_len(nrow(genome_pairs))) {
    genome1 <- genome_pairs[i,1]
    genome2 <- genome_pairs[i,2]
    
    cat(sprintf("%s vs %s\n", genome1, genome2))
    protein1 <- file.path(out, paste0(genome1, ".proteins.fa"))
    protein2 <- file.path(out, paste0(genome2, ".proteins.fa"))
    blat_output <- file.path(out, paste0(genome1, sep, genome2, ".blat"))
    # protein2 is the database, protein1 is the query
    system(sprintf("blat -t=prot -q=prot -out=blast8 %s %s %s",
                   protein2, protein1, blat_output))
    if (mercator) {
      hit_output <- replace_ext(blat_output, "hits", level=1)
      system(sprintf("blat2hits < %s > %s", blat_output, hit_output))
    } else {
      blat2hits(blat_output)
    }
  }
  
  return(invisible(TRUE)) 
}

