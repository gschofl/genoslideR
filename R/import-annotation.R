# We read the following features for each gbRange
#
# (Required) start, width, strand
# name, gene, synonym(=locus_tag), geneID, proteinID/GI, type, product
#
# Metadata neccessary for annotationList
# accession/identifier, description, length

#' @importFrom rentrez efetch
#' @importFrom rentrez esearch
#' @importFrom rentrez esummary
#' @importFrom rentrez docsum
#' @importFrom rentrez content
#' @importFrom rmisc strip_ext
#' @importFrom rmisc xvalue
#' @importFrom rmisc strsplitN
#' @importFrom rmisc %||%
#' @importFrom biofiles as.gbLocation
NULL

import_annotation_from_ptt <- function(file = anno[1], seqid = NULL) {
  
  skip <- sum(count.fields(file, sep="\t") < 9)
  ptt <- scan(file, skip = skip + 1, quote="",
              quiet=TRUE, sep="\t",
              what = list(Location = character(),
                          Strand = character(),
                          Length = character(),
                          PID = character(),
                          Gene = character(),
                          Synonym = character(),
                          Code = character(),
                          COG = character(),
                          Product = character()))
  
  location <- strsplit(ptt[["Location"]], "..", fixed=TRUE)
  start <- as.integer(vapply(location, "[", 1L, FUN.VALUE=character(1),
                             USE.NAMES=FALSE))
  width <- as.integer(vapply(location, "[", 2L, FUN.VALUE=character(1),
                             USE.NAMES=FALSE)) - start + 1L
  strand <- vapply(ptt[["Strand"]], switch, '+'=1L, '-'=-1L, NA_integer_,
                   FUN.VALUE=integer(1), USE.NAMES=FALSE)
  
  gene <- sub('^-$', NA, ptt[["Gene"]])
  synonym <- sub('^-$', NA, ptt[["Synonym"]])
  names <- ifelse(is.na(gene), synonym, gene)
  proteinID <- ptt[["PID"]]
  product <- ptt[["Product"]]

  ## missing from ptt files
  l <- length(start)
  geneID <- character(l)
  type <- rep("CDS", l)

  if (is.null(seqid)) {
    # attempt to fetch accession number
    x <- content(efetch(proteinID[1], "protein", "gp", "xml")) 
    x <- xvalue(x, '//GBQualifier[contains(GBQualifier_name, "coded_by")]/GBQualifier_value')
    if (all_empty(x)) {
      stop("No accession number could be retrieved as identifier for these annotations. Provide an identifier.")
    } else {
      seqid <- strip_ext(accession(as.gbLocation(x)))
    }
  }
  
  create_GRanges(seqid, start, width, strand, names,
                 type = type, gene = gene, synonym = synonym,
                 geneID = geneID, proteinID = proteinID, product = product)
}


import_annotation_from_tbl <- function(file = files[2], seqid = NULL, sep="\t") {
  
  header <- unlist(strsplit(readLines(file, n =  1), sep))
  if (all(c("start", "end") %ni% header)) {
    stop("The annotation file must contain 'start' and 'end' values")
  }
  
  table <- read.table(file, sep=sep, header=TRUE, quote="", as.is=TRUE)
  start <- as.integer(table[["start"]])
  width <- as.integer(table[["end"]]) - start + 1L
  table[["start"]] <- table[["end"]] <- NULL
  
  if (is.null(table[["strand"]])) {
    strand <- rep(1L, length(start)) 
  } else {
    strand <- as.integer(table[["strand"]])
    table[["strand"]] <- NULL
  }
  
  if (is.null(table[["names"]])) {
    names <- paste0("feature", seq_along(start))
  } else {
    names <- as.character(table[["names"]])
    table[["names"]] <- NULL
  }
  
  create_GRanges(seqid, start, width, strand, names, table)
}


import_annotation_from_ftb <- function(file) {
  
  l <- readLines(file, n=1)
  seqId <- strsplitN(l, split="\\s+", 2)
  m <- regexpr("[A-Za-z]{2}([A-Za-z_])?\\d+(\\.\\d)?", seqId)
  seqid <- strip_ext(regmatches(seqId, m))
  
  ft <- scan(file, sep="\t", comment.char=">", quiet=TRUE,
             quote="", fill=TRUE,
             what=list(start = character(),
                       end = character(),
                       key = character(),
                       qualifier = character(),
                       value = character()))
  
  pos_idx <- which(nzchar(ft[["start"]]))
  keys <- ft[["key"]][pos_idx]
  gene_idx <- pos_idx[keys == "gene"]
  start <- as.integer(gsub('<|>', '', ft[["start"]][gene_idx]))
  end <- as.integer(gsub('<|>', '', ft[["end"]][gene_idx]))
  strand <- ifelse(start > end, '-', '+')
  tmp_start <- ifelse(strand == '+', start, end)
  end <- ifelse(strand == '+', end, start)
  start <- tmp_start
  width <- end - start + 1L
  
  feature_idx <- Map(seq.int, gene_idx + 1, c(gene_idx[-1] - 1, length(ft[["start"]])))
  type <- product <- synonym <- gene <- character(length(feature_idx))
  proteinID <- geneID <- names <- character(length(feature_idx))
  for (i in seq_along(feature_idx)) {
    j <- feature_idx[[i]]   
    quals <- ft[["qualifier"]][j]
    vals <- ft[["value"]][j]
    product[i] <- paste(vals[quals == "product"] %||% '', collapse=", ")
    synonym[i] <- vals[quals == "locus_tag"] %||% ''
    gene[i] <- vals[quals == "gene"] %||% ''
    proteinID[i] <- strsplitN(vals[quals == "protein_id"] %||% '', "\\|", 2)
    db_xref <- vals[quals == "db_xref"] %||% ''
    geneID[i] <- strsplitN(db_xref[grepl("^GeneID:", db_xref)] %||% '', ':', 2) 
    names[i] <- gene[i] %||% synonym[i]
    type[i] <- ft[["key"]][j][which(nzchar(ft[["key"]][j]))[1]]
  }
  
  create_GRanges(seqid, start, width, strand, names,
                 type = type, gene = gene, synonym = synonym,
                 geneID = geneID, proteinID = proteinID, product = product)
}


import_annotation_from_gff <- function(file = anno[1], features = c("CDS", "RNA")) {
  
  l <- readLines(file)
  nlines <- -1L
  fline <- 0L
  if (any(l == "##FASTA")) {
    fline <- which(l == "##FASTA")[1L]
    nlines <- fline - 1L
  }
  
  con <- textConnection(l[if (nlines > 0) seq_len(nlines) else seq_along(l)])
  on.exit(close(con))
  nfields <- count.fields(con, sep = "\t", comment.char="#")
  
  if (!all(nfields == 9))
    stop("Malformed line(s) ", which(fields != 9))
  
  rm(l, nfields) # free memory
  
  # load data frame
  gff <- scan(file, nlines = nlines, sep = "\t",
              comment.char = "#", na.strings = ".",
              quote = "", quiet = TRUE,
              what = list(seqid = character(),
                          source = character(),
                          type = character(),
                          start = integer(),
                          end = integer(),
                          score = numeric(),
                          strand = character(),
                          phase = integer(),
                          attributes = character()))
  
  f_idx <- grep(paste0(features, collapse="|"), gff[["type"]], ignore.case=TRUE)
  ranges <- IRanges(gff[["start"]], gff[["end"]])
  ovl <- as(findOverlaps(ranges[f_idx,], IntervalTree(ranges), type="equal"),
            "data.frame")
  p_idx <- tapply(ovl[[2]], as.factor(ovl[[1]]), "[", 1)
  
  if (length(f_idx) != length(p_idx))
    stop("Probably malformed gff")

  start <- gff[["start"]][f_idx]
  width <- gff[["end"]][f_idx] - start + 1L
  type <- gff[["type"]][f_idx]
  strand <- vapply(gff[["strand"]][f_idx], switch, '+'=1L, '-'=-1L, NA_integer_,
                   FUN.VALUE=integer(1), USE.NAMES=FALSE)
  attr <- parse_attr(gff$attributes)
  
  ## misc_features have type region in gff; for the Seqinfo we want only the
  ## gff region with the GenBank key Source
  src_idx <- which(vapply(attr, "[", "gbkey", FUN.VALUE=character(1)) == "Src")
  if (length(src_idx) > 1) {
    warning("More then 1 'Source' field in ", basename(file))
    src_idx <- src_idx[1]
  } else if (length(src_idx) < 1) {
    warning("No 'Source' field in ", basename(file))
  }
  seqid <- strip_ext(gff[['seqid']][src_idx])
  
  f_attr <- attr[f_idx]
  p_attr <- attr[p_idx]
  product <- vapply(f_attr, "[", "product", FUN.VALUE=character(1))
  proteinID <- vapply(f_attr, "[", "protein_id", FUN.VALUE=character(1))
  geneID <- vapply(f_attr, "[", "Dbxref", FUN.VALUE=character(1))
  geneID <- vapply(lapply(strsplit(geneID, ","), function (x) {
    setNames(strsplitN(x, ":", 2), strsplitN(x, ":", 1))
  }), "[", "GeneID", FUN.VALUE=character(1))

  names <- vapply(p_attr, "[", "Name", FUN.VALUE=character(1))
  synonym <- vapply(p_attr, "[", "locus_tag", FUN.VALUE=character(1))
  gene <- vapply(p_attr, "[", "gene", FUN.VALUE=character(1))
 
  create_GRanges(seqid, start, width, strand, names,
                 type = type, gene = gene, synonym = synonym,
                 geneID = geneID, proteinID = proteinID, product = product)
}


get_attr <- function (gff, attribute) {
  anno <- annotation(gff)
  split <- strsplit(anno[["attributes"]], ";")
  sapply(attribute, extract_attr, string=split)
}


extract_attr <- function (id, string) {
  split <- paste0(id, "=")
  idx <- vapply(string, pmatch, x=split, FUN.VALUE=integer(1))
  good <- !is.na(idx)
  id <- mapply(function (s, i) strsplit(s, "=")[[c(i, 2L)]],
               s=string[good], i=idx[good], USE.NAMES=FALSE)
  idx[good] <- id
  idx
}


parse_attr <- function (anno) {
  key_val_pairs <- strsplit(anno, ";")
  
  kvp <- key_val_pairs[[1]]
  lapply(key_val_pairs, function (kvp) {
    kv <- strsplit(kvp, "=")
    setNames(object=vapply(kv, "[", 2L, FUN.VALUE=character(1)),
             nm=vapply(kv, "[", 1L, FUN.VALUE=character(1)))
  })
}


create_GRanges <- function (seqid, start, width, strand, names, ...) {
  
  seqinfo <- tryCatch({
    x <- docsum(esummary(esearch(seqid, "nuccore")))
    Seqinfo(seqnames = unlist(x["Caption"], use.names=FALSE),
            seqlengths = as.numeric(unlist(x["Length"], use.names=FALSE)),
            genome = unlist(x["Title"], use.names=FALSE))
  }, error = function (e) {
    Seqinfo(seqnames=seqid)
  })
  
  GRanges(seqnames=Rle(seqid),
          ranges=IRanges(start = start, width = width, names = names),
          strand=Rle(strand), ..., seqinfo = seqinfo) 
}


