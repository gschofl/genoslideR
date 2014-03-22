#' @importFrom genomeIntervals readGff3 getGffAttribute
NULL


get_position_from_gff <- function(gff) {
  start <- vector(mode="integer", length=nrow(gff))
  end <- vector(mode="integer", length=nrow(gff))
  pos <- which(gff@annotation$strand == "+")
  neg <- which(gff@annotation$strand == "-")
  start[pos] <- gff[pos, 1] - 1
  end[pos]   <- start[pos] + ((gff[pos, 2] - start[pos])/3)*3
  end[neg]   <- gff[neg, 2]
  start[neg] <- end[neg] - ((end[neg] - (gff[neg, 1] - 1))/3)*3
  list(start=start, end=end)
}


escape <- function (s) {
  s <- gsub(" ", "%20", s)
  s
}


unescape <- function (s) {
  s <- gsub("%20", " ", s)
  s
}


gff2anchors <- function(gff_file, anchor_file, id_as_anchor_num = FALSE) {
  gff <- readGff3(gff_file)
  cds <- gff[gff@annotation$type == "CDS", ]
  pos <- get_position_from_gff(cds)
  anchor.df <- data.frame(anchorNum = if (id_as_anchor_num) {
                            escape(getGffAttribute(cds, attribute="ID"))
                          } else {
                            seq_len(nrow(cds))
                          },
                          seqname=cds@annotation$seq_name,
                          strand=cds@annotation$strand,
                          start=pos$start,
                          end=pos$end,
                          n=rep(1, nrow(cds)))
  write.table(x=anchor.df, file=anchor_file, sep="\t", quote=FALSE,
              col.names=FALSE, row.names=FALSE)
}


pairwiseCombine <- function (x) {
  first_combin <- utils::combn(x, 2)
  recip_combin <- rbind(first_combin[2,], first_combin[1,])
  cbind(first_combin, recip_combin)
}


reciprocal_blat <- function(genomes, merc_sdb, merc_gff, merc_out,
                            sep = "-", removeOverlappingCDS = FALSE,
                            opts = list()) {
  if (length(genomes) < 2) {
    stop("At least two genomes must be specified")
  }
  # Make sure that a gff and sdb file are present for all genomes
  for (genome in genomes) {
    gff_file <- file.path(merc_gff, paste0(genome, ".gff"))
    sdb_file <- file.path(merc_sdb, paste0(genome, ".sdb"))
    if (!file.exists(gff_file)) {
      stop(sprintf("gff file for %s does not exist", sQuote(genome)))
    }
    if (!file.exists(sdb_file)) {
      stop(sprintf("SDB file for %s does not exist", sQuote(genome)))
    }
  }
  # Make anchor, chrom, and protein sequence file
  # genome <- genomes[1]
  for (genome in genomes) {
    gff_file <- file.path(merc_gff, paste0(genome, ".gff"))
    sdb_file <- file.path(merc_sdb, paste0(genome, ".sdb"))
    
    # make chrom file
    cat(sprintf("Making chromosome file for %s...\n", sQuote(genome)))
    chrom_file <- file.path(merc_out, paste0(genome, ".chroms"))
    system(sprintf("sdbList -l %s > %s", sdb_file, chrom_file))
    cat("done\n")
    
    # make anchor file
    cat(sprintf("Making anchors for %s...\n", sQuote(genome)))
    anchor_file <- file.path(merc_out, paste0(genome, ".anchors"))
    if (removeOverlappingCDS) {
      system(sprintf("cat %s | gffRemoveOverlaps -fCDS -l | gff2anchors > %s",
                     gff_file, anchor_file))
    } else {
      gff2anchors(gff_file, anchor_file, id_as_anchor_num = FALSE)
    }
    cat("done\n")
    
    # Make protein file
    cat(sprintf("Extracting protein sequences from %s...\n", sQuote(genome)))
    prot_file <- file.path(merc_out, paste0(genome, ".proteins.fa"))
    system(sprintf("anchors2fa %s < %s > %s", sdb_file, anchor_file, prot_file))
    cat("done\n")
  }
  
  # BLAT protein anchors pairwise
  cat("BLATing anchors ...\n")
  genome_pairs <- t(pairwiseCombine(genomes))
  for (i in seq_len(nrow(genome_pairs))) {
    genome1 <- genome_pairs[i, 1]
    genome2 <- genome_pairs[i, 2]
    
    cat(sprintf("%s vs %s\n", genome1, genome2))
    protein1 <- file.path(merc_out, paste0(genome1, ".proteins.fa"))
    protein2 <- file.path(merc_out, paste0(genome2, ".proteins.fa"))
    blat_output <- file.path(merc_out, paste0(genome1, sep, genome2, ".blat"))
    # protein2 is the database, protein1 is the query
    opts <- paste0(sprintf("-%s=%s", names(opts), unlist(opts)), collapse = " ")
    blat_exec <- sprintf("blat -t=prot -q=prot -out=blast8 %s %s %s %s",
                         opts, protein2, protein1, blat_output)
    system(blat_exec)
    if (removeOverlappingCDS) {
      hit_output <- replace_ext(blat_output, "hits", level = 1)
      system(sprintf("blat2hits < %s > %s", blat_output, hit_output))
    } else {
      blat2hits(blat_output)
    }
  }
  
  return(invisible(TRUE)) 
}

