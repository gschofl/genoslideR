#' @importFrom rmisc replace_ext
#' @importFrom rmisc not.null
#' @importFrom rmisc all_empty
NULL

#' Ab initio genome annotation using glimmer3
#'
#' @param fasta Path to a genome file(s) in fasta format
#' @param opts a named list of options for glimmer3
#' @param ... named values interpreted as options for glimmer3
#' @param cleanup Clean up intermediate files.
#' 
#' @return Character vector. Path to glimmer3 file(s).
#' @export
glimmer3 <- function(fasta, opts = list(o=50, g=110, t=30),
                     ..., cleanup = TRUE) {
  
  ## check dependencies
  hasDependencies(c("glimmer3", "elph", "long-orfs", "extract",
                    "build-icm", "start-codon-distrib"))
  opts <- merge_list(opts, list(...))
  ncores <- detectCores() - 1 
  glimmer_files <- mclapply(fasta, glimmer, opts=opts, cleanup=cleanup,
                            mc.cores=ncores)
  return(unlist(glimmer_files))
}

glimmer <- function (f, opts, cleanup) {
  
  glimmer_dir <- file.path(dirname(f), paste0(".glimmer_", strip_ext(basename(f))))
  if (file.exists(glimmer_dir))
    unlink(glimmer_dir, recursive=TRUE)
  dir.create(glimmer_dir)
  
  orfs <- longorfs(f, glimmer_dir)
  train <- extract(f, orfs)
  icm <- build_icm(train)
  run1 <- run_glimmer(f, icm, tag="run1", opts)
  
  # Get the training coordinates from the first predictions
  coords <- replace_ext(run1[["glimmer3"]], "coords", level=2)
  writeLines(readLines(run1[["glimmer3"]])[-1], coords)
  
  # Create a position weight matrix (PWM) from the regions
  # upstream of the start locations in coords
  upstream <- extract_upstream(f, coords, len=25, sep=0)
  pwm <- pwm(upstream, 6)
  motif <- replace_ext(upstream, "motif", level=1)
  
  write.table(ncol(pwm), file=motif, row.names=FALSE, col.names=FALSE)
  write.table(matrix(pad(pwm, 7), ncol=ncol(pwm)), 
              file=motif, row.names=rownames(pwm), col.names=FALSE,
              quote=FALSE, append=TRUE)

  # Determine the distribution of start-codon usage in coords
  startuse <- system(paste("start-codon-distrib -3", f, coords), intern = TRUE)
  
  # Run second glimmer
  run <- run_glimmer(f, icm, "", opts=opts, b=motif, P=startuse)
  
  if (cleanup)
    unlink(c(orfs, train, icm, run1, coords, motif, upstream))
  
  return(invisible(run[["glimmer3"]]))
}


# Find long, non-overlapping orfs to use as a training set
longorfs <- function (f, outdir = glimmer_dir) {
  message("Finding long orfs for training")
  longorfs <- file.path(outdir, replace_ext(basename(f), "longorfs", level=1))
  st <- system(paste("long-orfs -n -t 1.15", f, longorfs),
               intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (not.null(attr(st, "status")))
    stop("Failed to find long-orf training set")
  
  return(invisible(longorfs))
}


# Extract the training sequences from the genome file
extract <- function(f, orfs) {
  message("Extracting training sequences")
  train <- replace_ext(orfs, "train", level=1)
  st <- system(paste("extract -t",f, orfs, ">", train), intern = TRUE)
  if (not.null(attr(st, "status")))
    stop("Failed to extract training sequences")
  
  return(invisible(train))
}


# Build the icm from the training sequences
build_icm <- function(train) {
  message("Building ICM")
  icm <- replace_ext(train, "icm", level=1)
  st <- system(paste("build-icm -r", icm, "<", train),
               intern = TRUE)
  if (not.null(attr(st, "status")))
    stop("Failed to build ICM")
  
  return(invisible(icm))
}


# Run Glimmer
run_glimmer <- function(f, icm, tag = "",
                        opts = list(o=50, g=110, t=30), ...) {
  
  message("Running Glimmer3")
  args <- merge_list(opts, list(...))
  run <- if (all_empty(tag)) {
    strip_ext(icm, level=1)
  } else {
    replace_ext(icm, tag, level=1)
  }
  
  if (any(nchar(names(args)) != 1)) {
    stop("Use short options with glimmer3")
  }
  
  st <- SysCall("glimmer3", args=args, style="unix", redirection=FALSE,
                stdin=paste(f, icm, run), intern=TRUE)
  
  if (not.null(attr(st, "status")))
    stop("Failed to run Glimmer3")
  
  file.rename(from=paste0(run, ".detail"), to=paste0(run, ".glimmer3.detail"))
  file.rename(from=paste0(run, ".predict"), to=paste0(run, ".glimmer3"))
  return(invisible(c(detail=paste0(run, ".glimmer3.detail"),
                     glimmer3=paste0(run, ".glimmer3"))))
}
  
  
extract_upstream <- function (f, coords, len = 25, sep = 0) {
  
  max_gene_len <- 100000
  out <- replace_ext(coords, "upstream", level=1)
  
  coords <- scan(file=coords, flush = TRUE, quiet = TRUE,
                 what=list(tag=character(),
                           start=numeric(),
                           end=numeric()))
  
  gene_len <- ifelse(coords[["start"]] < coords[["end"]],
                     1 + coords[["end"]] - coords[["start"]],
                     1 + coords[["start"]] - coords[["end"]])
  
  strand <- ifelse(coords[["start"]] < coords[["end"]], 1, -1)
  strand[gene_len > max_gene_len] <- strand[gene_len > max_gene_len]*-1
  
  upstream <- paste0(pad(coords[["tag"]], n=2, where="right"),
                     pad(coords[["start"]] - strand * (sep + len), n=8),
                     pad(coords[["start"]] - strand * (sep + 1), n=8))
  
  st <- system(paste("extract", f, "- >", out), input=upstream,
               intern = TRUE)
  
  if (not.null(attr(st, "status")))
    stop("Failed to extract upstream regions")
  
  return(invisible(out))
}


pwm <- function (mfa = upstream, len = 6, ...) {
  
  args <- list(...)
  elph <- SysCall(paste0("elph ", mfa, " LEN=", len), b=TRUE,
                  args=args, intern=TRUE)
  
  idx <- grep(pattern="^Motif counts:", elph) + 1:4
  pwm <- as.matrix(read.table(textConnection(elph[idx]),
                              header=FALSE, row.names=1))
  rownames(pwm) <- strsplitN(rownames(pwm), ":", 1)
  pwm
}


