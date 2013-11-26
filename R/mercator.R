#' @importFrom ape read.dna cbind.DNAbin dist.dna clustal muscle tcoffee nj bionj fastme.bal write.tree
#' @importFrom parallel detectCores mcmapply
NULL

#' Run mercator to build a homology map and create orthologous segments
#' that can be aligned using FSA
#' 
#' @details
#' \code{mercator} will create a hidden directory \sQuote{.mercator}
#' in the parent directory of the submitted sequence files and place
#' all results of the intermediate computational steps in this directory.
#' 
#' The results of a \code{mercator} run are stored in a directory
#' \sQuote{segments} in the parent directory of the submitted sequence files.
#' The \sQuote{segments} directory should contain numbered subdirectories
#' containing orthologous genomic segments in multi fasta format, a file
#' named \sQuote{map}, and a file named \sQuote{genomes}.
#' 
#' \code{mercator} segements can be aligned using the
#' \code{\link{alignSegments}} function, merged into a single multi-fasta
#' alignment file and imported into R as a \code{\linkS4class{DNAStringSet}}
#' using \code{\link{importAlignment}} or imported to R as an
#' \code{\linkS4class{annotatedAlignment}} using \code{\link{annotatedAlignment}}.
#' 
#' @param seq_files Path to sequence files ('fasta' or 'genbank').
#' @param anno_files (Optional) Path to annotation files. If no
#' annotation files are provided an \emph{ab initio} annotation
#' using \code{\link{glimmer3}} is performed.
#' @param anno_type Type of annotation ('glimmer3', 'genbank', 'gff', 
#' 'ptt', or 'ftable').
#' @param glimmeropts Options passed to \code{\link{glimmer3}} if an
#' \emph{ab initio} annotation is performed.
#' @param removeOverlappingCDS Exclude overlapping CDS for orthology mapping.
#' @param mask Softmask sequence before aligning using
#' \code{\link{maskSequence}}.
#' @param wd Working directory. Defaults to the parent directory of
#' the provided sequence files.
#' @seealso \code{\link{alignSegments}}, \code{\link{importAlignment}},
#' \code{\link{annotatedAlignment}}.
#' @export
mercator <- function (seq_files, anno_files = NULL, anno_type = "glimmer3",
                      glimmeropts = list(o=50, g=110, t=30), removeOverlappingCDS = TRUE,
                      mask = TRUE, wd = NULL) {
  annotation <- match.arg(anno_type, c("glimmer3", "genbank", "gbk", "gff", "ptt",
                                       "ftable"))
  ## Check dependencies for BLAT
  has_dependencies(c("sdbList", "gffRemoveOverlaps", "gff2anchors", "anchors2fa",
                     "blat", "blat2hits"))
  ## Check dependencies for mercator
  has_dependencies(c("fa2sdb", "mercator", "sdbAssemble", "phits2constraints",
                     "makeAlignmentInput", "sdbExport", "muscle", "omap2hmap",
                     "makeBreakpointGraph", "makeBreakpointAlignmentInput",
                     "findBreakpoints", "breakMap", "hmap2omap", "omap2coordinates"))
  # if not specified set the working directory to the parent directory of all
  # sequence files
  if (is.null(wd)) {
    wd <- Reduce(function(lhs, rhs) compactNA(rhs[match(lhs, rhs)]),
                 strsplit(dirname(seq_files), .Platform$file.sep))
    wd <- normalizePath(paste(wd, collapse=.Platform$file.sep))
  }
  if (!all(is_fasta(seq_files) | is_gbk(seq_files))) {
    stop("Sequences must be provided in FASTA or GenBank format")
  }
  ## Extract FASTA from GenBank files and reassign the gbk files as anno_files
  ## if annotation = 'genbank' and anno_files = NULL
  gbk_files <- which(is_gbk(seq_files))
  if (length(gbk_files) > 0) {
    if ((annotation == 'genbank' || annotation == 'gbk') && is.null(anno_files)) {
      anno_files <- seq_files[gbk_files]
    }
    outdirs <- dirname(seq_files)
    seq_files[gbk_files] <- gbk2fna(seq_files[gbk_files], outdirs[gbk_files])
  }
  # mask sequence files
  if (mask) {
    seq_files <- maskSequence(normalizePath(seq_files))
  }
  seq_files <- normalizePath(seq_files)
  
  # generate annotation if necessary
  if (is.null(anno_files)) {
    if (annotation != "glimmer3")
      stop("No annotation files provided")
    else
      anno_files <- glimmer3(normalizePath(seq_files))
  }
  anno_files <- normalizePath(anno_files)

  if (length(anno_files) != length(seq_files))
    stop("Unequal number of sequence and annotation files")
  
  if (any(strip_ext(basename(anno_files)) %ni% strip_ext(basename(seq_files))))
    stop("Names of annotation and sequence files must match.")
  
  merc <- file.path(wd, ".mercator")
  if (file.exists(merc)) {
    warning("A '.mercator' directory exists in ", wd, immediate.=TRUE)
    ans <- readline("Overwrite [y/n]? ")
    if (ans == "y") {
      unlink(merc, recursive=TRUE)
    } else {
      return(NULL)
    }
  }
  
  merc_gff <- file.path(merc, "gff")
  merc_fas <- file.path(merc, "fasta")
  merc_sdb <- file.path(merc, "sdb")
  merc_hit <- file.path(merc, "hits")
  for (dir in c(merc_gff, merc_fas, merc_sdb, merc_hit)) {
    create_if_not_exists(dir, type="dir", recursive=TRUE)
  }
  # generate mercator-readable gff files
  gff_files <- gff_for_mercator(anno_files, annotation, merc_gff)
  
  # generate mercator-readable fasta files
  fna_files <- fna_for_mercator(seq_files, merc_fas)
  
  sdb_files <- file.path(merc_sdb, replace_ext(basename(fna_files), "sdb"))
  invisible(mapply(function(x, y) {
    system(paste0("fa2sdb ", x, " < ", y))
  }, x = sdb_files, y = fna_files))
  
  # compare all of the exon sequences pairwise and generate the proper
  # input files for Mercator.
  reciprocal_blat(genomes = strip_ext(basename(fna_files)),
                  merc_sdb = merc_sdb, merc_gff = merc_gff, merc_out = merc_hit,
                  sep = "-", removeOverlappingCDS = removeOverlappingCDS)
  
  # run Mercator on the input files.
  # This generates a directory 'segments' with subdirectories
  # containing the sequences to be aligned for each orthologous segment set     
  # identified by the orthology map.
  segment_dir <- run_mercator(wd)
  
  message("Next use 'alignSegments()' to generate alignments for each of the orthologous segments produced by mercator")
  return(segment_dir)
}


gff_for_mercator <- function (f, type, wd) {
  if (missing(wd)) {
    stop("No working directory provided")
  }
  type <- match.arg(type, c("gff", "ptt", "genbank", "gbk", "ftable","glimmer3"))
  out <- switch(type,
                gff=gff2gff(f, wd),
                ptt=ptt2gff(f, wd),
                genbank=gbk2gff(f, wd),
                gbk=gbk2gff(f, wd),
                ftable=ftb2gff(f, wd),
                glimmer3=glimmer2gff(f, wd))
  return(invisible(out))
}


gff2gff <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    l <- readLines(f[i])
    nlines <- -1L
    if (any(l == "##FASTA")) {
      nlines <- which(l == "##FASTA")[1L] - 1L
    }
    rm(l) # free memory
    
    # load data frame
    gff <- scan(f[i], nlines = nlines, sep = "\t", comment.char = "#", 
                quote = "", quiet = TRUE, na.strings = '.',
                what = list(seqid = character(),
                            source = character(),
                            type = character(),
                            start = integer(),
                            end = integer(),
                            score = numeric(),
                            strand = character(),
                            phase = integer(),
                            attributes = character()))
    
    cds_idx <- which(gff$type == "CDS")
    len <- length(cds_idx)
    seqid <- strip_ext(basename(f[i]))
    gff <- data.frame(stringsAsFactors=FALSE,
                      seqid = rep(seqid, len), 
                      source = gff[["source"]][cds_idx],
                      type = gff[["type"]][cds_idx],
                      start = gff[["start"]][cds_idx],
                      end = gff[["end"]][cds_idx],
                      score = rep(".", len),
                      strand = gff[["strand"]][cds_idx],
                      phase = rep("0", len),
                      attributes = gff[["attributes"]][cds_idx])
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    write.table(gff, outfile, quote=FALSE, sep="\t", row.names=FALSE, 
                col.names=FALSE)
    outfiles <- c(outfiles, outfile)
    rm(gff) ## free memory
  }
  
  return(invisible(outfiles))
}


ptt2gff <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    skip <- sum(count.fields(f[i], sep="\t") < 9)
    ptt <- scan(f[i], skip = skip + 1, quote="",
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
    
    seqid <- strip_ext(basename(f[i]))
    loc <- strsplit(ptt[["Location"]], "..", fixed=TRUE)
    len <- length(loc)
    gff <- data.frame(stringsAsFactors=FALSE,
                      seqid = rep(seqid, len), 
                      source = rep("ptt", len),
                      type = rep("CDS", len),
                      start = vapply(loc, "[", 1L, FUN.VALUE=character(1)),
                      end = vapply(loc, "[", 2L, FUN.VALUE=character(1)),
                      score = rep(".", len),
                      strand = ptt[["Strand"]],
                      phase = rep("0", len),
                      attributes = paste0("ID=cds", seq.int(0, len - 1),
                                          ";product=", ptt[["Product"]]))
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    write.table(gff, outfile, quote=FALSE, sep="\t", row.names=FALSE, 
                col.names=FALSE)
    outfiles <- c(outfiles, outfile)
    rm(ptt,gff) ## free memory
  }
  
  return(invisible(outfiles))
}


ftb2gff <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    l <- readLines(f[i], n=1)
    seqid <- strsplitN(l, split="\\s+", 2)
    m <- regexpr("[A-Za-z]{2}([A-Za-z_])?\\d+(\\.\\d)?", seqid)
    seqid <- strip_ext(regmatches(seqid, m))
    
    ft <- scan(f[i], sep="\t", comment.char=">", quiet=TRUE,
               quote="", fill=TRUE,
               what=list(start = character(),
                         end = character(),
                         key = character(),
                         qualifier = character(),
                         value = character()))
    
    pos_idx <- which(nzchar(ft[["start"]]))
    gene_idx <- pos_idx[ft[["key"]][pos_idx] == "CDS"]
    start <- as.integer(gsub('<|>', '', ft[["start"]][gene_idx]))
    end <- as.integer(gsub('<|>', '', ft[["end"]][gene_idx]))
    strand <- ifelse(start > end, '-', '+')
    tmp_start <- ifelse(strand == '+', start, end)
    end <- ifelse(strand == '+', end, start)
    start <- tmp_start
    
    len <- length(start)
    gff <- data.frame(stringsAsFactors=FALSE,
                      seqid = rep(seqid, len), 
                      source = rep("ftb", len),
                      type = rep("CDS", len),
                      start = start,
                      end = end,
                      score = rep(".", len),
                      strand = strand,
                      phase = rep("0", len),
                      attributes = paste0("ID=cds", seq.int(0, len - 1)))
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    write.table(gff, outfile, quote=FALSE, sep="\t", row.names=FALSE, 
                col.names=FALSE)
    outfiles <- c(outfiles, outfile)
    rm(ftb,gff) ## free memory
  }
   
  return(invisible(outfiles))
}


gbk2gff <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    gbk <- readLines(f[i])
    seqid <- strsplit(gbk[1], split = "\\s+")[[1]][2]
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    features_start <- which(substr(gbk, 1, 8) == "FEATURES") + 1  
    features_end <- which(grepl("ORIGIN|CONTIG", substr(gbk, 1, 6))) - 1 %||% length(gbk) - 1
    features <- gbk[features_start:features_end]
    cds_idx <- which(substr(features, 6, 8) == "CDS")
    loc <- trim(strsplitN(features[cds_idx], "CDS", 2))
    ## discard joins
    loc <- grep("join", loc, value=TRUE, invert=TRUE)
    strand <- ifelse(grepl("complement", loc), '-', '+')
    loc <- strsplit(regmatches(loc, regexpr('<?\\d+\\.\\.>?\\d+', loc)),
                    split="..", fixed=TRUE)
    start <- gsub('<|>', '', vapply(loc, `[`, 1, FUN.VALUE=character(1)))
    end <- gsub('<|>', '', vapply(loc, `[`, 2, FUN.VALUE=character(1)))
    len <- length(loc)
    gff <- data.frame(stringsAsFactors=FALSE,
                      seqid = rep(seqid, len), 
                      source = rep("genbank", len),
                      type = rep("CDS", len),
                      start = start,
                      end = end,
                      score = rep(".", len),
                      strand = strand,
                      phase = rep("0", len),
                      attributes = paste0("ID=cds", seq.int(0, len - 1)))
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    write.table(gff, outfile, quote=FALSE, sep="\t", row.names=FALSE, 
                col.names=FALSE)
    outfiles <- c(outfiles, outfile)
    rm(gbk,gff) ## free memory
  }
  return(invisible(outfiles))
}


glimmer2gff <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    glim <- scan(f[i], comment.char = ">", quote="",
                quiet=TRUE, sep="",
                what = list(id = character(),
                            start = integer(),
                            end = integer(),
                            frame = character(),
                            score = numeric()))
    seqid <- strip_ext(basename(f[i]))
    len <- length(glim[["id"]])
    dir <- glim[["end"]] > glim[["start"]]
    start <- ifelse(dir, glim[["start"]], glim[["end"]])
    end <- ifelse(dir, glim[["end"]], glim[["start"]])
    strand <- ifelse(dir, "+", "-")
    
    max_gene_len <- ceiling((max(end) - min(start))*0.9)
    idx <- end - start > max_gene_len
    gff <- data.frame(stringsAsFactors=FALSE,
                      seqid = rep(seqid, len)[!idx], 
                      source = rep("glimmer", len)[!idx],
                      type = rep("CDS", len)[!idx],
                      start = start[!idx],
                      end = end[!idx],
                      score = rep(".", len)[!idx],
                      strand = strand[!idx],
                      phase = rep("0", len)[!idx],
                      attributes = paste0("ID=", glim[["id"]],
                                          ";frame=", glim[["frame"]],
                                          ";score=", glim[["score"]])[!idx])
    
    outfile <- file.path(wd[i], paste0(seqid, ".gff"))
    write.table(gff, outfile, quote=FALSE, sep="\t", row.names=FALSE, 
                col.names=FALSE)
    outfiles <- c(outfiles, outfile)
  }
  
  return(invisible(outfiles))
}


fna_for_mercator <- function (f, wd) {
  if (missing(wd)) {
    stop("No working directory provided")
  }
  type <- if (all(is_fasta(f))) {
    "fna"
  } else if (all(is_gbk(f))) {
    "gbk"
  } else {
    "other"
  }
  out <- switch(type,
                fna=fna2fna(f, wd),
                gbk=gbk2fna(f, wd),
                other=stop("Sequence files must be either in FASTA or GBK format"))
  
  return(invisible(out))
}


fna2fna <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    fna <- readLines(f[i])
    seqid <- strip_ext(basename(f[i]))
    outfile <- file.path(wd[i], paste0(seqid, ".fna"))
    outcon <- file(outfile, open="w")
    on.exit(close(outcon))
    fna[1] <- sprintf(">%s", seqid)
    writeLines(fna, outcon)
    outfiles <- c(outfiles, outfile)
  }
  return(invisible(outfiles))
}


gbk2fna <- function (f, wd) {
  outfiles <- character()
  wd <- if (length(wd) == 1) rep(wd, length(f)) else wd
  for (i in seq_along(f)) {
    gbk <- readLines(f[i])
    seqid <- strsplit(gbk[1], split = "\\s+")[[1]][2]
    header <- sprintf(">%s", seqid)
    outfile <- file.path(wd[i], paste0(seqid, ".fna"))
    outcon <- file(outfile, open="w")
    writeLines(header, outcon)
    ori_start <- which(substr(gbk, 1, 6) == "ORIGIN") + 1
    ori_end <- which(substr(gbk, 1, 2) == "//") - 1
    gbk <- gbk[seq.int(ori_start, ori_end)]
    fna <- vapply(gbk, function(x) {
      toupper(paste0(substr(x, 11, 20), substr(x, 22, 31),
                     substr(x, 33, 42), substr(x, 44, 53),
                     substr(x, 55, 64), substr(x, 66, 75),
                     collapse = ""))
    }, character(1), USE.NAMES = FALSE)
    writeLines(fna, outcon)
    close(outcon)
    outfiles <- c(outfiles, outfile)
  }
  return(invisible(outfiles))
}


run_mercator <- function (wd) {
  merc     <- file.path(wd, ".mercator")
  merc_fna <- file.path(merc, "fasta")
  merc_sdb <- file.path(merc, "sdb")
  merc_hit <- file.path(merc, "hits")
  merc_out <- file.path(merc, "out")
  dir.create(merc_out)
  genomes <- strip_ext(dir(merc_fna))
  system(sprintf("mercator -i %s -o %s %s", 
                 merc_hit, merc_out, paste(genomes, collapse=' ')))
  
  # comparative scaffolding of genome sequences
  sdb_files <- dir(merc_sdb, full.names=TRUE)
  for (sdb in sdb_files) {
    system(sprintf("sdbAssemble %s %s < %s", sdb,
                   file.path(merc_out, basename(sdb)),
                   file.path(merc_out, paste0(strip_ext(basename(sdb)), ".agp"))))
  }
  
  # Run phits2constraints to generate the "constraints" file
  system(sprintf("phits2constraints -i %s -m %s < %s > %s",
                 merc_hit, merc_out,
                 file.path(merc_out, "pairwisehits"),
                 file.path(merc_out, "constraints")))
    
  # Generate a guide tree in Newick format
  guide_tree(wd)
  
  # refine breakpoints
  find_breakpoints(wd)
  
  # generate alignment input for FSA
  map_path <- dir(merc_out, "\\<better\\.map\\>")
  segment_dir <- file.path(wd, "segments")
  if (file.exists(segment_dir)) {
    unlink(segment_dir, recursive=TRUE)
  }
  dir.create(segment_dir)
  system(sprintf("makeAlignmentInput --map=%s %s %s",
                 map_path, merc_out, segment_dir))
  
  return(invisible(segment_dir))
}


guide_tree <- function(wd) {
  merc <- file.path(wd, ".mercator")
  merc_hits <- file.path(merc, "hits")
  merc_out <- file.path(merc, "out")
  # check if all necessary files are present  
  f <- dir(merc_out, full.names=TRUE)
  if (sum(runs_pos <- grepl(pattern="\\<runs\\>", basename(f))) != 1L) {
    stop("'runs' file missing")
  }
  if (!any(grepl(pattern="\\.sdb\\>", basename(f)))) {
    stop("sdb-files missing")
  }
  merc_tree <- file.path(merc_out, "tree")
  create_if_not_exists(merc_tree, type="dir")
  
  runs <- read.table(f[runs_pos], header=FALSE, colClasses="integer")
  runs <- runs[complete.cases(runs),]
  names(runs) <- scan(f[grepl("genomes$", f)], what="character", quiet=TRUE)
  
  if (nrow(runs) == 0) {
    stop("There are no orthologous genes to generate a guide tree")
  }
  if (nrow(runs) > 20) {
    idx <- sample(seq_len(nrow(runs)), 20, replace = FALSE)
    runs <- runs[idx,]
  } 

  match_genes(runs, merc)
  
  tree <- make_tree(merc_tree, "muscle", "K80", "nj")
  ape::write.tree(tree, file=file.path(merc_out, "treefile"))
  
  return(invisible(tree))
}


match_genes <- function(runs, merc) {
  
  merc_out <- file.path(merc, "out")
  merc_hits <- file.path(merc, "hits")
  tree_dir <- file.path(merc_out, "tree")
  genomes <- names(runs)
  
  anchor <- list()
  for (i in seq_along(genomes))
    anchor[[i]] <- read.table(file.path(merc_hits, paste0(genomes[i], ".anchors")),
                              sep="\t")
  
  sdbs <- dir(merc_out, "\\.sdb$", full.names=TRUE)
  for (i in seq_along(genomes)) {
    gene_idx <- match(runs[,i], anchor[[i]][,1])
    seqname <- as.character(anchor[[i]][gene_idx, 2])
    strand <- as.character(anchor[[i]][gene_idx, 3])
    start <- anchor[[i]][gene_idx, 4]
    end <- anchor[[i]][gene_idx, 5]
    fasta <- list()
    for (j in seq_along(gene_idx)) {
      cat(system(sprintf("sdbExport -r --name=%s %s %s %s %s %s",
                         paste0(genomes[i], ":", runs[j,i]),
                         sdbs[i], seqname[j], start[j], end[j],
                         strand[j]),
                 intern=TRUE),
          file=file.path(tree_dir, paste0("orf", j, ".fa")),
          sep="\n", append=TRUE)
    }
  }
  
  return(invisible(TRUE))
}


#' [INTERNAL] conctruct NJ or BIONJ tree from one or more multi-fasta files
#' 
#' @param fasta_dir Directory containing (a) multi-FASTA file(s).
#' @param align One of 'clustal', 'muscle', or 'tcoffee'.
#' @param dist.model See \code{\link[ape]{dist.dna}}
#' @param tree One of 'nj', 'bionj', or 'fastme'.
#' 
#' @keywords internal
make_tree <- function (fna_dir, align="muscle", dist.model="K80", tree="bionj") {
  
  fna <- dir(fna_dir, "\\.fa$", full.names=TRUE)
  align.fun <- match.fun(match.arg(align, c('muscle', 'clustal', 'tcoffee')))
  tree.fun <- match.fun(match.arg(tree, c('bionj', 'nj', 'fastme.bal')))
  
  alignment <- list()
  for (i in seq_along(fna)) {
    dna <- read.dna(fna[i], format="fasta")
    alignment[[i]] <- align.fun(dna)
    cat(sprintf("Aligning %s ...\n", basename(fna[i])))
  }
  
  o <- NULL
  for (aln in alignment) {
    o <- c(o, list(order(labels(aln))))
  }
  for (i in seq_along(alignment)) {
    alignment[[i]] <- alignment[[i]][o[[i]],]
  }
  
  concat_aln <- do.call("cbind.DNAbin", c(alignment, check.names=FALSE))
  rownames(concat_aln) <- 
    sapply(strmatch(pattern="[^>].[^:]+", rownames(concat_aln), capture=FALSE),
           "[", 1)

  dist <- dist.dna(concat_aln, model=dist.model)
  tree <- tree.fun(dist)
  return(tree)
}

find_breakpoints <- function (wd) {
  
  merc <- file.path(wd, ".mercator")
  merc_out <- file.path(merc, "out")
  merc_break <- file.path(merc_out, "breakpoints")
  dir.create(merc_break)
  
  f <- dir(merc_out, full.names=TRUE)
  if (sum(map_pos <- grepl(pattern="\\<pre\\.map\\>", basename(f))) != 1L)
    stop("'pre.map' is missing")
  
  if (sum(tree_pos <- grepl(pattern="\\<treefile\\>", basename(f))) != 1L)
    stop("'treefile' is missing")
  
  pwd <- getwd()
  setwd(merc_out)
  on.exit(setwd(pwd))
  system("omap2hmap genomes < pre.map > pre.h.map")
  system("makeBreakpointGraph --remove-colinear pre.h.map treefile")
  system("makeBreakpointAlignmentInput --out-dir=breakpoints")
  fsaAlignSegmentDirs(initdir="breakpoints", seqfile="seqs.fasta",
                      outfile="fsa.mfa", fsa.opts=list(logfile="fsa.log"),
                      ncores=detectCores() - 1)
  system("findBreakpoints -r 200 --alignmentFile=fsa.mfa pre.h.map treefile edges breakpoints > myBreakpoints")
  system("breakMap myBreakpoints < pre.h.map > better.h.map")
  system("hmap2omap genomes < better.h.map > better.map")
  system("omap2coordinates < better.map > coordinates")
  
  return(invisible(TRUE))
}


fsaAlignSegmentDirs <- function (initdir = ".", seqfile = "seqs.fasta",
                                 outfile = "fsa.mfa", constraints = "cons",
                                 skip.completed = TRUE, fsa.opts = list(),
                                 ncores = detectCores()) {
  assert_that(has_command('fsa')) 
  segments <- normalizePath(dir(initdir, "^\\d+$", full.names=TRUE))
  segments <- segments[order(as.numeric(split_path(segments)))]
  segdirs <- mcmapply(alignSegmentDir, segment = segments,
                      MoreArgs = list(seqfile = seqfile, outfile = outfile,
                                      constraints = constraints, skip.completed = skip.completed,
                                      fsa.opts = fsa.opts),
                      mc.cores=ncores, USE.NAMES=FALSE)
  
  return(invisible(segdirs))
}


alignSegmentDir <- function(segment, seqfile, constraints, outfile = "fsa.mfa",
                            skip.completed = TRUE, fsa.opts = list()) {
  
  if (skip.completed && file.exists(file.path(segment, outfile)) && 
      file.info(file.path(segment, outfile))$size > 0) {
    return(NULL)
  }
  
  # check for seqfile
  if (!file.exists(file.path(segment, seqfile))) {
    stop("Directory ", sQuote(basename(segment)), " does not have the required seqfile ",
         sQuote(seqFile), .call=FALSE)
  }
  
  # check for optional constraints file
  if (!file.exists(file.path(segment, constraints))) {
    warning("Directory ", sQuote(basename(segment)), " does not have constraints",
            call.=FALSE, immediate.=TRUE)
  } else {
    fsa.opts <- merge_list(fsa.opts, list(mercator=constraints))
  }
  
  cwd <- getwd()
  setwd(segment)
  fsa(seqfile, outfile, fsa.opts)
  setwd(cwd)
}


#' Sequence alignment with fsa
#' 
#' @param seqfile Sequence files in fasta format.
#' @param outfile Multifasta alignment file.
#' @param opts a named list of options for fsa.
#' @param ... Named values interpreted as options for fsa.
#' 
#' @export
fsa <- function (seqfile, outfile = "fsa.mfa", opts = list(), ...) {
  
  if (missing(seqfile)) {
    system("fsa --help")
    return(invisible(NULL))
  }
  
  if (!all(file.exists(seqfile)))
    stop("Can not open input files")
  
  if (length(seqfile) > 1)
    infiles <- paste(infiles, collapse=" ")
  
  args <- merge_list(opts, list(...))
  SysCall("fsa", args = args, stdin = seqfile, stdout = outfile,
          redirection = FALSE, style = "gnu")
}

