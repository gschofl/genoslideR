#' @importFrom RCurl getURL
#' @importFrom RCurl curlOptions
#' @importFrom RCurl getCurlHandle
#' @importFrom rmisc not_empty
NULL

#' Retrieve bacterial genomes from NCBI FTP site
#' 
#' Get genomes from \url{ftp://ftp.ncbi.nih.gov/genomes/Bacteria/}
#'
#' @param which A regular expression to match the directories
#' of completed bacterial genomes.
#' @param what A regular expression to match file extensions. Possible
#' values are \emph{gbk} (Genbank), \emph{gff} (GFF3 file), \emph{fna}
#' (FASTA), \emph{ptt} (protein table).
#' @param where Local target directory for downloading
#' @param ignore.case Case sensitive pattern matching if \code{FALSE}
#'
#' @export
ncbi_bacteria <- function(which, what="gbk|gff|fna", where="~/Bacteria", 
                          ignore.case=TRUE) {
  
  if (missing(which)) {
    stop("Provide a regexp or an index to delimit genomes")
  }
  
  ncbi_url <- "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/"
  bact_dirs <- strsplit(getURL(ncbi_url, ftplistonly=TRUE), "\\n")[[1]]
  
  if (is.numeric(which)) {
    target <- bact_dirs[which]
  } else {
    target <- bact_dirs[grep(which, bact_dirs, ignore.case=ignore.case)] %||% NA
  }
  
  print(target)
  idx <- readline("Download from these directories (y/n/index) [y]: ")
  if (not_empty(idx) && idx != "y") {
    if (idx == "n") {
      target <- NULL
    } else {
      idx <- tryCatch(eval(parse(text=idx)), error = function (e) {
        print(sprintf("%s is not a valid R expression", sQuote(idx)))
      })
      target <- target[idx]
    }
  }
  
  ## generate target dirs
  target_dirs <- path.expand(file.path(where, target))
  target_dirs <- target_dirs[!file.exists(target_dirs)]
  x <- sapply(target_dirs, dir.create, recursive=TRUE)
  
  ## fetch data
  opts <- curlOptions(timecondition=TRUE, ftp.use.epsv=FALSE,
                      forbid.reuse=TRUE, filetime=TRUE)
  curl <- getCurlHandle(.opts=opts)
  urls <- sprintf("%s%s/", ncbi_url, target)
  for (url in urls) {
    files <- strsplit(getURL(url, ftplistonly=TRUE), "\\n")[[1]]
    files <- files[grep(what, files)]
    f_urls <- sprintf("%s%s", url, files) 
    
    for(f in f_urls) {
      to <- strsplitN(f, .Platform$file.sep, c(1,2), "end")
      out <- path.expand(file.path(where, to))
      timevalue <- unclass(file.info(out)$ctime)
      time_val <- curlOptions(timevalue=timevalue)
      contents <- getURL(f, curl=curl, .opts=time_val, verbose=TRUE)
      if (nchar(contents) > 0) {
        cat(contents, file=out)
      }
    }
  }
}

#' Install external dependencies for genoslideR
#' 
#' Works only on Ubuntu (and possibly Debian)
#' 
#' @param sudo if \code{TRUE} install with sudo
#'
#' @export
install_genoslider_dependencies <- function (sudo = TRUE) {
  
  if (!hasCommand("zsh")) {
    stop("Install 'zsh' and then try again ...")
  }
  
  if (is_empty(Sys.getenv("TERM"))) {
    term_emul <- c("gnome-terminal", "konsole", "xterm")
    term_emul <- term_emul[hasCommand(term_emul)][1]
    if (is.na(term_emul)) {
      stop("No terminal emulator found")
    }
    term <- paste(term_emul, "-e")
  } else {
    term <- ""
  }

  exec <- system.file("src", "install_dependencies_ubuntu.sh", package="genoslideR")
  if (sudo)
    system(paste(term, "'sudo", exec, "'"))
  else
    system(paste(term, "'", exec, "'"))

}

#' Quickly look up the size of a file
#' 
#' @param files Path(s) to files.
#' @param unit Return file size in bytes (B, default), kilobytes (kB),
#' megabytes (MB), gigabytes (GB), or terabytes (TB). 
#' 
#' @export
file_size <- function(files, unit = NULL) {
  unit <- match.arg(unit, c("B", "kB", "MB", "GB", "TB"))
  file.info(files)$size / switch(unit,
                                 B=1,
                                 kB=1024,
                                 MB=1024^2,
                                 GB=1024^3,
                                 TB=1024^4)
}
