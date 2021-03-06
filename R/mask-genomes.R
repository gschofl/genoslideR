#' Softmask repeats and low-complexity regions using
#' RepeatScout and RepeatMasker
#' 
#' \code{maskSequence} will create a hidden directory \sQuote{.masked}
#' in the parent directory of the submitted fasta files and place all
#' all intermediate files as well as the final softmasked sequence file
#' into this directory.
#' 
#' If multiple fasta files share the same parent directory, the softmasked
#' files are placed in individual subdirectories within \sQuote{.masked}.
#' 
#' @param fasta Path to fasta file(s).
#' 
#' @return Character vector. Path to softmasked file(s).
#' @export
maskSequence <- function(fasta) {
  ## check dependencies
  has_dependencies(c("awk", "build_lmer_table", "RepeatScout", "filter-stage-1.prl",
                     "RepeatMasker", "nmerge", "faSoftMask"))
  ncores <- detectCores() - 1 
  masked <- mclapply(fasta, repeatmasker, mc.cores=ncores)
  return(unlist(masked))
}


repeatmasker <- function (f) {
  pwd <- dirname(f)
  fn <- basename(f)
  ## check if a softmasked file is already exists
  mask_dir <- file.path(pwd, paste0(".masked_", strip_ext(fn)))
  softmasked <- file.path(mask_dir, paste0(strip_ext(fn), '.softmasked.fa'))
  if (file.exists(softmasked)) {
    message(sQuote(basename(softmasked)), " already exists.\nRemove the ",
            sQuote(basename(mask_dir)), " directory to start over.")
    return(softmasked)
  }
  dir.create(mask_dir)
  # generate a temporary infile with capital letters
  tmp_f <- paste0(f, "~")
  system(paste("awk 'FNR == 1'", f, ">", tmp_f))
  system(paste("awk 'FNR > 1'", f, "| tr '[a,c,g,t,X]' '[A,C,G,T,x]' >>", tmp_f))
  
  # move into work directory and softlink to infile
  setwd(mask_dir)
  on.exit(setwd(pwd))
  unmasked <- replace_ext(fn, "unmasked.fa", level=1)
  file.link(tmp_f, unmasked)
  
  # RepeatScout
  #
  # tabulate the frequency of all l-mers in the sequence
  # to be analysed. use: build_lmer_table
  freq <- replace_ext(fn, "freq", level=1)
  cmd <- paste("build_lmer_table -sequence", unmasked, "-freq", freq, "-v")
  system(cmd)
  
  rep <- replace_ext(fn, "rep", level=1)
  log <- replace_ext(fn, "log", level=1)
  cmd <- paste("RepeatScout -sequence", unmasked, "-output", rep, "-freq",
               freq, "-stopafter 500 -vv >", log)
  system(cmd)
  
  # filter out low-complexity and tandem elements
  filtered <- replace_ext(fn, "rep.filtered", level=1)
  cmd <- paste("cat", rep, "| filter-stage-1.prl >", filtered)
  system(cmd)
  
  # Repeatmasker
  #
  # scan for interspersed repeats using file.rep.filtered (if nonempty)
  # as library of sequences to be masked
  interspersed <- replace_ext(fn, "interspersed.fa", level=1)
  if (file.info(filtered)$size > 0) {
    file.link(unmasked, interspersed)
    cmd <- paste("RepeatMasker -no_is -nolow -lib", filtered, interspersed)
    system(cmd)
  }
  
  # scan for low complexity repeats
  lowcomp <- replace_ext(fn, "lowcomp.fa", level=1)
  file.link(unmasked, lowcomp)
  cmd <- paste("RepeatMasker -no_is -noint", lowcomp)
  system(cmd)
  
  # merge mask into one hardmasked file
  # nmerge is part of the ABBlast distribution /path/to/abblast/filter
  isp_masked <- paste0(interspersed, ".masked")
  lwc_masked <- paste0(lowcomp, ".masked")
  hardmasked <- replace_ext(fn, "hardmasked.fa", level=1)
  
  if (file.exists(interspersed)) {
    cmd <- paste("nmerge", isp_masked, lwc_masked, ">", hardmasked)
    system(cmd)
  } else if (file.info(lowcomp)$size == 0) {
    file.copy(lwc_masked, hardmasked)
  } else {
    file.link(unmasked, hardmasked)
  }
  
  # create softmasked file
  # faSoftMask is part of Colin Dewey's source code package
  # /path/to/cndsrc-###/utils/
  softmasked <- replace_ext(fn, "softmasked.fa", level=1)
  cmd <- paste("faSoftMask", unmasked, hardmasked, ">", softmasked)
  system(cmd)
  
  unlink(c(tmp_f, unmasked, interspersed, lowcomp))
  return(normalizePath(softmasked))
}




