#' @importFrom stringr str_trim str_extract str_count
#' @importFrom Matrix Matrix summary
#' @importFrom plyr dlply .
NULL


#' Find Reciprocal Best Hits using Blat 
#' 
#' @details
#' \code{rbh} will create a hidden directory \sQuote{.rbh}
#' in the parent directory of the submitted sequence files and place
#' all results of the intermediate computational steps in this directory.
#' 
#' The results of a \code{rbh} run are stored in a directory
#' \sQuote{rbh} in the parent directory of the submitted sequence files.
#' 
#' @param seq_files Path to sequence files ('fasta' or 'genbank').
#' @param anno_files (Optional) Path to annotation files. If no
#' annotation files are provided an \emph{ab initio} annotation
#' using \code{\link{glimmer3}} is performed.
#' @param anno_type Type of annotation ('glimmer3', 'genbank', 'gff', 
#' 'ptt', or 'ftable').
#' @param blatopts Options passed to \code{blat}.
#' @param glimmeropts Options passed to \code{\link{glimmer3}} if an
#' \emph{ab initio} annotation is performed.
#' @param wd Working directory. Defaults to the parent directory of
#' the provided sequence files.
#' @export
rbh <- function(seq_files,
                anno_files = NULL,
                anno_type = "glimmer3",
                blatopts = list(),
                glimmeropts = list(o=50, g=110, t=30),
                wd = NULL) {
  annotation <- match.arg(anno_type, c("glimmer3", "genbank", "gbk", "gff", "ptt", "ftable"))
  ## Check dependencies for BLAT
  has_dependencies(c("blat", "fa2sdb", "sdbList", "anchors2fa"))
  # if not specified set the working directory to the parent directory of all
  # sequence files
  if (is.null(wd)) {
    wd <- Reduce(function(lhs, rhs) compactNA(rhs[match(lhs, rhs)]),
                 strsplit(dirname(seq_files), .Platform$file.sep))
    wd <- normalizePath(paste(wd, collapse = .Platform$file.sep))
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
  seq_files <- normalizePath(seq_files)
  
  # generate annotation if necessary
  if (is.null(anno_files)) {
    if (annotation != "glimmer3")
      stop("No annotation files provided")
    else
      anno_files <- glimmer3(seq_files, opts = glimmeropts)
  }
  anno_files <- normalizePath(anno_files)
  assert_that(length(anno_files) == length(seq_files))
  if (any(strip_ext(basename(anno_files)) %ni% strip_ext(basename(seq_files)))) {
    stop("Names of annotation and sequence files don't match.")
  }
  rbh <- file.path(wd, ".rbh")
  if (file.exists(rbh)) {
    warning("A '.rbh' directory exists in ", wd, immediate. = TRUE, call. = FALSE)
    ans <- readline("Overwrite [y/n]? ")
    if (ans == "y") {
      unlink(rbh, recursive = TRUE)
    } else {
      return(NULL)
    }
  }
  
  rbh_gff <- file.path(rbh, "gff")
  rbh_fas <- file.path(rbh, "fasta")
  rbh_sdb <- file.path(rbh, "sdb")
  rbh_hit <- file.path(rbh, "hits")
  for (dir in c(rbh_gff, rbh_fas, rbh_sdb, rbh_hit)) {
    create_if_not_exists(dir, type = "dir", recursive = TRUE)
  }
  # generate mercator-readable gff files
  gff_files <- gff_for_mercator(f = anno_files, type = annotation, wd = rbh_gff)
  
  # generate mercator-readable fasta files
  fna_files <- fna_for_mercator(f = seq_files, wd = rbh_fas)
  
  sdb_files <- file.path(rbh_sdb, replace_ext(basename(fna_files), "sdb"))
  invisible(mapply(function(x, y) {
    system(paste0("fa2sdb ", x, " < ", y))
  }, x = sdb_files, y = fna_files))
  
  # compare all of the CDS sequences
  reciprocal_blat(genomes = strip_ext(basename(fna_files)),
                  merc_sdb = rbh_sdb, merc_gff = rbh_gff, merc_out = rbh_hit,
                  sep = "-", removeOverlappingCDS = FALSE, opts = blatopts)
  
  rbh_dir <- file.path(wd, "rbh")
  if (file.exists(rbh_dir)) {
    unlink(rbh_dir, recursive=TRUE)
  }
  dir.create(rbh_dir)
  file.copy(from = dir(rbh_hit, full.names = TRUE), to = rbh_dir)
  message("Next use 'orthologyMatrix()' to contruct an orthology matrix")
  return(rbh_dir)
}


#' Construct orthology matrices
#' 
#' Constructs pairwise orthology matrices for all combinations of species
#' as provided by the \code{base_names} argument based on reciprocal best
#' BLAT hits. Use \code{\link{mergeOrthologyMatrix}} to merge these
#' pairwise comparison into a single data frame.
#' 
#' @param blat_dir Directory where the anchor and hit files produced by
#' \code{link{mercator}} live.
#' @param base_names Species for which an orthology matrix should be
#' constructed. Must correspond to the relevant filenames minus extension.
#' If \code{NULL} all are selected
#' @param pident_threshold Threshold for (cumulative) percentage of sequence
#' identity
#' @param paln_threshold Threshhold for (cumulative) percentage of alignment
#' length (query coverage).
#' @param name_sep String that separates two compared species in filenames.
#' @param log Write logfiles
#' 
#' @return A list of class \sQuote{orthoMatrix}, which is essentially a 
#' list of \code{\link[Matrix]{dgCMatrix}}es with attributes
#' \sQuote{labels}, \sQuote{nGenes}, and \sQuote{combinations} attached.
#' 
#' @export
orthologyMatrix <- function(blat_dir,
                            base_names = NULL,
                            pident_threshold = 60,
                            paln_threshold = 40,
                            name_sep = "-",
                            log = TRUE) {
  if (is.null(base_names)) {
    base_names <- strsplitN(basename(dir(blat_dir, "*anchors$", full.names = TRUE)),
                            split = "\\.", n = 1)
  }
  hitfiles <- dir(blat_dir, "*hits$", full.names = TRUE)
  anchorfiles <- normalizePath(paste0(file.path(blat_dir, base_names), ".anchors"))
  anchors <- count_lines(anchorfiles)
  gene_num <- as.numeric(str_trim(str_extract(string = anchors, " *\\d+")))
  species <- data.frame(stringsAsFactors = FALSE, base_names = base_names, gene_num = gene_num)
  
  Omat <- list()
  combinations <- utils::combn(seq_len(nrow(species)), 2)
  for (combi in seq_len(ncol(combinations))) {
    #combi <- 1
    if (log) {
      log_file <- file.path(blat_dir, "logs", sprintf("orthology_matrix_%s.log", combi))
      create_if_not_exists(dirname(log_file), type = "dir")
      create_if_not_exists(log_file, type = "file")
    }
    i <- combinations[1, combi] 
    j <- combinations[2, combi]
    file_ij   <- file.path(blat_dir, paste0(species[i, 1], name_sep, species[j, 1], ".hits"))
    file_ji   <- file.path(blat_dir, paste0(species[j, 1], name_sep, species[i, 1], ".hits"))
    n_genes_i <- species[i, 2]
    n_genes_j <- species[j, 2]
    label_i   <- species[i, 1]
    label_j   <- species[j, 1]
    hits_ij   <- read.table(file_ij, as.is=TRUE) 
    hits_ji   <- read.table(file_ji, as.is=TRUE)
    names(hits_ij) <- names(hits_ji) <- c("query", "subject", "pident", "paln", "score", "eval")
    ## remove all hits with an e-value below the cutoff
    if (sum(bad_hits <- hits_ij$pident < pident_threshold |
                        hits_ij$paln < paln_threshold) > 0L) {
      hits_ij <- hits_ij[!bad_hits, ]  
    }  
    if (sum(bad_hits <- hits_ji$pident < pident_threshold |
                        hits_ji$paln < paln_threshold) > 0L) {
      hits_ji <- hits_ji[!bad_hits, ]
    }
    # write to logfile
    valid_queries_i   <- sort(unique(hits_ij$query))
    valid_queries_j   <- sort(unique(hits_ji$query))
    match_found_for_i <- numeric(0)
    match_found_for_j <- numeric(0)
    
    # split by query i and subject j
    ij_q_split <- dlply(hits_ij, .(query))
    ij_s_split <- dlply(hits_ij, .(subject))
    # split by subject i and query j
    ji_s_split <- dlply(hits_ji, .(subject))
    ji_q_split <- dlply(hits_ji, .(query))
    omat <- Matrix(data = FALSE, nrow = n_genes_i, ncol = n_genes_j,
                   dimnames = list(seq_len(n_genes_i), seq_len(n_genes_j)))
    cat(sprintf("Extracting orthologs from:\n\t%s (%s valid queries) vs. %s (%s valid queries) ...\n\n",
                label_i, length(valid_queries_i), label_j, length(valid_queries_j)))
    
    if (log) {
      cat(sprintf("Logging %s vs %s\n\n", label_i, label_j), file = log_file)
    }
    for (i in seq_along(ij_q_split)) {
      if (log) {
        cat(sprintf("Split #%s:\n", i), file = log_file, append = TRUE)
      }
      ## catch an error thrown if a reciprocal match is missing because
      ## it has been unilaterally filtered by the pident and paln cutoffs
      tryCatch({
        qi <- ij_q_split[[i]]
        si_pat <- paste0(paste0("\\<", unique(qi$subject), "\\>"), collapse="|")
        
        sj <- ij_s_split[grep(si_pat, names(ij_s_split))]
        sj <- do.call("rbind", sj)
        rownames(sj) <- NULL
        qi_pat <- paste0(paste0("\\<", unique(sj$query), "\\>"), collapse="|")
        
        qj <- ji_q_split[grep(si_pat, names(ji_q_split))]
        qj <- do.call("rbind", qj)
        rownames(qj) <- NULL
        sj_pat <- paste0(paste0("\\<", unique(c(sj$query, qj$subject)), "\\>"), collapse="|")
        
        si <- ji_s_split[grep(sj_pat, names(ji_s_split))]
        si <- do.call("rbind", si)
        rownames(si) <- NULL
        
        si <- si[, c("subject", "query", "pident", "paln", "score", "eval")]
        qj <- qj[, c("subject", "query", "pident", "paln", "score", "eval")]
        names(qi) <- names(sj) <- names(si) <- names(qj) <- c("i", "j", "pident", "paln", "score", "eval")
        matches <- do.call("rbind", list(qi, sj, si, qj))  
        if (length(unique(matches$i)) > 1L || length(unique(matches$j)) > 1L) {
          x <- matches[order(matches[, "eval"], -matches[, "score"], -matches[, "pident"]),]
          while (nrow(x) > 0L) {         
            if (log) {
              write.table(x, file = log_file, row.names = FALSE, col.names = TRUE,
                          append = TRUE, quote = FALSE, sep = "\t")
            }
            i <- x[1, "i"]
            j <- x[1, "j"] 
            x <- x[!x[, "i"] == i, ]
            x <- x[!x[, "j"] == j, ]
            if (log) {
              cat(sprintf("i = %s, j = %s\n\n", i, j), file=log_file, append=TRUE)
            }
            omat[i, j] <- TRUE
            match_found_for_i <- c(match_found_for_i, i)
            match_found_for_j <- c(match_found_for_j, j)
          }
        } else {
          i <- unique(matches$i)
          j <- unique(matches$j)
          if (log) {
            write.table(matches, file=log_file, row.names=FALSE, col.names=TRUE,
                        append=TRUE, quote=FALSE, sep="\t")
            cat(sprintf("i = %s, j = %s\n\n", i, j), file=log_file, append=TRUE)
          }
          omat[i,j] <- TRUE
          match_found_for_i <- c(match_found_for_i, i)
          match_found_for_j <- c(match_found_for_j, j)
        }
      }, 
               error = function(e) {
                 s <- paste(sprintf("%s: %s", names(ij_q_split[[i]]), ij_q_split[[i]]), collapse = ", ")
                 message(sprintf("No reciprocal match for %s\n", s))
                 if (log) {
                   cat(sprintf("No reciprocal match for \n%s\n\n", s), file = log_file,
                       append = TRUE)
                 }
               })
    }
    
    if (log) {
#       cat(sprintf("Valid queries from species i:\n%s\n\n",
#                   linebreak(paste(valid_queries_i, collapse=" "), width=19)),
#           file=log_file, append=TRUE)
#       cat(sprintf("Matches found for species i:\n%s\n\n",
#                   linebreak(paste(sort(match_found_for_i), collapse=" "),
#                             width=19)), file=log_file, append=TRUE)
#       cat(sprintf("Valid queries from species j:\n%s\n\n",
#                   linebreak(paste(valid_queries_j, collapse=" "), width=19)),
#           file=log_file, append=TRUE)
#       cat(sprintf("Matches found for species f:\n%s\n\n",
#                   linebreak(paste(sort(match_found_for_j), collapse=" "),
#                             width=19)), file=log_file, append=TRUE)
    }
    
    Omat[[combi]] <- omat
    rm(omat)
  }
  structure(
    Omat,
    labels = species$base_names,
    ngenes = species$gene_num,
    combinations = combinations,
    anchors = anchorfiles,
    class = "OrthoMatrix"
  )
}

#' Plot an orthology matrix
#' 
#' @param omat A list of orthology matrices.
#' @param i Index.
#' @method plot OrthoMatrix
#' @export
plot.OrthoMatrix <- function(x, i=NULL, ...) {
  plot.omat <- function(x, i) {
    xlab <- attr(x, "labels")[attr(x, "combinations")[2,i]]
    ylab <- attr(x, "labels")[attr(x, "combinations")[1,i]]
    im <- image(as(x[[i]], "dgCMatrix"), xlab=xlab, ylab=ylab, lwd=.1,
                sub=sprintf("%d x %d genes", dim(x[[i]])[1], dim(x[[i]])[2]),
                main=sprintf("Combination %s vs. %s",
                             attr(x, "combinations")[1,i],
                             attr(x, "combinations")[2,i]))
    im
  }
  
  if (is.null(i)) {
    stop("No i provided")
    #     l <- length(omat)
    #     l <- 4
    #     n_col <- ceiling(sqrt(l))
    #     n_row <- ceiling(sqrt(l)) - if (l %% n_col >= n_col) 1 else 0
    #     left <- seq(0, n_col - 1)/n_col
    #     right <- seq(1, n_col)/n_col
    #     bottom <- seq(0, n_row - 1)/n_row
    #     top <- seq(1, n_row)/n_row
    # 
    #     for (i in seq_along(omat))
    #       im <- plot.omat(omat, i)
    #       print(im, c(0,.5,.5,1), more=TRUE) }
  } else if (i <= length(omat)) {
    im <- plot.omat(x, i)
    print(im)
  } else stop("i is out of range")
}

##' Merge a list of pairwise orthology matrices into an OrthoFrame
##' 
##' Merges all pairwise comparisons of genomes as provided by
##' \code{\link{orthologyMatrix}} into a single overall orthology Matrix.
##' (an object of class \sQuote{\bold{OrthoFrame}}). It checks for
##' conflicts and provides those as an attribute \sQuote{\emph{conflicts}}.
##' The attribute \sQuote{\emph{labels}} holds the species names.
##' 
##' @param omat An \code{orthoMatrix} object
##' @return An \code{OrthoFrame} object.
##' @export
mergeOrthologyMatrix <- function(omat) {
  merge_list <- list()
  for (i in seq_along(omat)) {
    combi <- attr(omat, "combination")[,i]
    species1 <- Matrix::summary(omat[[i]])$i
    species2 <- Matrix::summary(omat[[i]])$j
    ortho.df <-  structure(data.frame(species1, species2), names=combi)
    ortho.df <-ortho.df[order(ortho.df[,1]),]
    merge_list[[i]] <- ortho.df
  }
  
  mergedMatrix <- merge_all(x=merge_list)
  
  ## construct OrthoFrame
  oFrame <- mergedMatrix[,order(as.numeric(colnames(mergedMatrix)))]
  ngenes <- attr(omat, "ngenes")
  for (i in seq(ncol(oFrame))) {
    orphans <- setdiff(seq_len(ngenes[i]), unique(oFrame[,i]))
    orph_mat <- matrix(rep(NA, length(orphans)*ncol(oFrame)), ncol=ncol(oFrame))
    orph_mat[,i] <- orphans
    colnames(orph_mat) <- colnames(oFrame)
    oFrame <- rbind(oFrame, orph_mat, deparse.level=0)
  }
  oFrame <- oFrame[do.call(order, oFrame),]
  rownames(oFrame) <- NULL
  structure(
    oFrame,
    labels    = attr(omat, "labels"),
    conflicts = check_conflicts(oFrame),
    ngenes    = attr(omat, "ngenes"),
    anchors   = attr(omat, "anchors"),
    class     = c("OrthoFrame", "data.frame")
  )
}

#### Helper functions for mergeOrthologyMatrix ####

merge_all <- function(x) {
  merged.df <- Reduce(function(x, y) {
    collapse_duplicates(merge(x, y, all=TRUE))
  }, x, accumulate=FALSE)
  merged.df <- merged.df[, order(names(merged.df))]
  merged.df <- merged.df[do.call(order, merged.df), ]
  rownames(merged.df) <- NULL
  merged.df
}

collapse_duplicates <- function(m) {
  SKIP <- 0
  #print(sprintf("SKIP = %s", SKIP))
  if (is.na(dup <- get_duplicate(m, SKIP))) {
    #print(sprintf("dup = %s", dup))
    return(m)
  } else {
    #print(sprintf("dup = %s", dup))
    #print(sprintf("SKIP = %s", SKIP))
    while (!is.na(dup)) {
      m <- merge_duplicate(m, dup, SKIP)
      dup <- get_duplicate(m, SKIP)
      #print(sprintf("dup = %s", dup))
      #print(sprintf("SKIP = %s", SKIP))
    }
    return(m)
  }
}

get_duplicate <- function(m, SKIP) {
  ## reset rownames
  rownames(m) <- NULL
  dup <- as.numeric(rownames(m[duplicated(m[,1], incomparables=NA),]))
  #print(sprintf("full dup = %s", paste(dup, collapse=" ")))
  if (length(dup) < 1L + SKIP)
    return(NA_integer_)
  else
    return(dup[1L + SKIP])
}

merge_duplicate <- function(m, dup, SKIP) {
  dup.row <- m[c(dup - 1, dup),]
  if (dup.row[1,1] != dup.row[2,1]) {
    stop("Didn't get duplicate")
  }
  ## check if in all columns one row has an NA entry.
  ## if that is the case collapse
  if (all(apply(is.na(dup.row)[,-1L], 2, any))) {
    dup.row <- c(dup.row[1,1], colSums(dup.row[,-1L], na.rm=TRUE))
    dup.row[dup.row == 0L] <- NA
    m[dup[1] - 1,] <- t(as.matrix(dup.row))
    m <- m[-dup[1],]
  } else {
    # SKIP <<- SKIP + 1
    assign("SKIP", SKIP + 1, pos=sys.frame(-1))
  }
  return(m)
}

check_conflicts <- function(merged_df) {
  conflicts <- apply(merged_df, 2, function(x) {
    runs <- rle(sort(x, na.last=TRUE))
    conflicts <- runs$values[runs$lengths > 1]
    if (length(conflicts) == 0) return(NULL)
    attr(conflicts, "times") <- runs$lengths[runs$lengths > 1] 
    conflicts
  })
  
  cc <- list() ## collection of conflicts
  for (i in seq_along(conflicts)) {
    if (is.null(conflicts[[i]])) {
      next
    } else {
      for (j in seq_along(conflicts[[i]])) {
        x <- merged_df[which(merged_df[,i] == conflicts[[i]][j]),]
        cc <- c(cc, list(x))
      }
    }
  }
  dups <- which(duplicated(unlist(lapply(cc, rownames))))
  cc <- do.call(rbind, cc)[-dups,]
  return(cc)
}

#### Methods for OrthoFrame objects ####

#' @method print OrthoFrame
#' @export
print.OrthoFrame <- function(x, ...) {
  NextMethod()
  cat("\nConflicts that may need to be solved manually:\n")
  print(attr(x, "conflicts"))
}

#' Access labels
#' @rdname labels
#' @export
labels <- function(object, ...) {
  UseMethod("labels")
}

#' @method labels OrthoFrame
#' @export
labels.OrthoFrame <- function(object, ...) {
  attr(object, "labels")
}

#' @rdname labels
#' @export
`labels<-` <- function(x, value) {
  UseMethod("labels<-")
}

#' @method "labels<-" OrthoFrame
#' @export
`labels<-.OrthoFrame` <- function(x, value) {
  if (ncol(x) != length(value)) {
    stop("Number of labels must corresspond to the number of columns in the OrthoFrame")
  }
  attr(x, "labels") <- as.character(value)
  return(x)
}

#' Access anchor files 
#' @export
anchors <- function(object, ...) {
  UseMethod("anchors")
}

#' @method anchors OrthoFrame
#' @export
anchors.OrthoFrame <- function(object) {
  res <- list()
  anchors <- attr(object, "anchors")
  for (a in anchors) {
    res <- c(res, list(read.table(a, header=FALSE, stringsAsFactors=FALSE)))
  }
  names(res) <- labels(object)
  res
}

#' @method $ OrthoFrame
#' @export
`$.OrthoFrame` <- function(x, name) {
  if (nchar(name) != ncol(x) || !all(strsplit(name, "")[[1]] %in% c("0","1","."))) {
    NextMethod()
  } else {
    pat <- str_trim(name)
    idx <- venn_groups(x)
    n.dots <- str_count(pat, "\\.")
    if (n.dots > 0L) {
      l <- list()
      for (i in seq_len(n.dots)) l[[i]] <- c(0,1)
      l <- expand.grid(l)
      spat <- strsplit(pat, "")[[1]]
      smat <- matrix(rep(spat, nrow(l)), nrow=nrow(l), byrow=TRUE)
      smat[, which(spat == ".")] <- as.matrix(l)
      patterns <- apply(smat, 1, paste0, collapse="")
    } else {
      patterns <- pat
    }
    do.call(rbind, lapply(patterns, function(pat) x[idx[[pat]], ]))
  }
}

#' @method [[ OrthoFrame
#' @export
`[[.OrthoFrame` <- function(x, ..., exact = TRUE) {
  arg <- list(...)[[1L]]
  if (is.character(arg) && nchar(arg) == ncol(x) && all(strsplit(arg, "")[[1]] %in% c("0","1","."))) {
    `$.OrthoFrame`(x, arg)
  } else {
    NextMethod()
  }
}

## construct a list of overlapping indices that
## can be used to find orthology groupings
orthoList <- function(x) {
  if (!is(x, "OrthoFrame")) {
    stop("x must be of class OrthoFrame")
  }
  rownames(x) <- NULL
  labels <- labels(x)
  cols <- ncol(x)
  idx <- as.integer(rep(rownames(x), cols))
  x <- unlist(x)
  names(x) <- NULL
  x[!is.na(x)] <- idx[!is.na(x)]
  x <- as.list(data.frame(stringsAsFactors=FALSE, matrix(x, ncol=cols)))
  x <- lapply(x, function(elem) {
    if (any(is.na(elem))) elem[!is.na(elem)] else elem 
  })
  structure(x, names=seq_along(x), labels=labels, class=c("OrthoList", "list"))
}


venn_groups <- function(x) {
  if (is(x, "OrthoFrame")) {
    x <- orthoList(x)
  }
  if (!is(x, "OrthoList")) {
    stop("Neither OrthoFrame nor OrthoList provided")
  }
  all_values <- unique(unlist(x))
  all_values <- all_values[!is.na(all_values)]
  grp <- character(0)
  for (i in all_values) {
    k <- vapply(x, function(j) i %in% j, logical(1))
    grp <- c(grp, paste(as.numeric(k), collapse=""))
  }
  df <- data.frame(group=grp, value=as.numeric(all_values))
  group_list <- NULL
  for (level in levels(df$group)) {
    group_list <- c(group_list, list(df[df$group == level, "value"]))    
  }
  structure(
    group_list,
    names = levels(df$group),
    labels = attr(x, "labels"),
    class = "VennList"
  )
}

