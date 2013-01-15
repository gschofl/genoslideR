require(devtools)
require(digest)
genoslider <- as.package("~/R/Projects/Devel/genoslideR/")
document(genoslider, clean=TRUE)
load_all(genoslider)






ncbi_bacteria("chlamyd.*", ".*", "~/local/workspace/Chlamydia")











# require("BiocGenerics")
# require("IRanges")
# require("GenomicRanges")
# require("Biostrings")
# require("stringr")
# require("rentrez")
# require("biofiles")
# require("parallel")
# 
# source("R/align-genomes.r")
# source("R/aligner.r")
# source("R/gapped_ranges.r")
# source("R/glimmer.r")
# source("R/import-alignment.r")
# source("R/import-annotation.r")
# source("R/mask-genomes.r")
# source("R/mercator.r")
# source("R/reciprocal-blat.r")
# source("R/slice-alignment.r")
# source("R/annotation.R")
# source("R/get_aln_range.R")
# source("R/utils.r")

# alignment für Jochen
# aln_dir <- normalizePath("~/bioinf/var/aln_jochen/")
# accn <- c("CP002806", "NC_000117", "NC_004552", "NC_003361",
#           "NC_007899", "NC_002620", "NC_015408", "CP001713") 
# 
# for (a in accn) {
#   out <- file.path(aln_dir, paste0(a, ".fna"))
#   write(efetch(a, "nuccore", rettype="fasta"), file=out)
# }
# 
# for (a in accn) {
#   out <- file.path(aln_dir, paste0(a, ".ftp"))
#   write(efetch(a, "nuccore", rettype="ft"), file=out)
# }
# 
# seq_files <- file.path(aln_dir, paste0(accn, ".fna"))
# anno_files <- file.path(aln_dir, paste0(accn, ".ftb"))
# 
# seg_dir <- mercator(seq_files, seq_type="fasta", anno_files,
#                     anno_type="ftable", mask=TRUE)
# 
# aln <- align_genomes(seg_dir=aln_dir)
# 
# chlam <- annotatedAlignment(aln, anno=anno_files, type="ftable")
# 
# save(chlam, file=file.path(aln_dir, "chlamAln.Rda"))

load(file.path(aln_dir, "chlamAln.Rda"))


cutAlignment <- function (aln = chlam, start, end, names = NULL, strand = NULL,
                          genome = 2) {
  
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
  
  if (!is.null(names)) {
    names <- ifelse(is.na(names), "", names)
  }
  
  cut_range <- GRanges(genome, IRanges(start=start, end=end, names=names))
  res <- get_aln_range(aln, cut_range)
  
  if (!is.null(strand)) {
    idx <- which(strand == "-")
    res[idx] <- lapply(res[idx], function (s) {
      reverseComplement(as(s, "DNAStringSet"))
    })
  }
  
  res
}

findAnnotation <- function (query, subject) {
  subj <- annotation(subject)
  query <- metadata(query)[["genomic_position"]]
  genomes <- seqlevels(query)
  anno <- GRangesList()
  for (g in genomes) {
    hits <- subjectHits(findOverlaps(query[[g]], subj[[g]], type="any"))
    anno[[g]] <- subj[[g]][hits,]
  }
  
  anno
}

s = start(chlam)[["NC_002620"]][340:350]
e = end(chlam)[["NC_002620"]][340:350]
n = locusTag(chlam)[["NC_002620"]][340:350]
x <- cutAlignment(aln=chlam, start=s, end=e, names=n, genome = "NC_002620")
findAnnotation(query=x[[1]], subject=chlam)



