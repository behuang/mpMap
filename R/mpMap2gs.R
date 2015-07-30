#' Export objects to Richard Mott's genome_scan
#'
#' Given an mpcross object, export it in files suitable for input
#' to Richard Mott's reconstruction and genome_scan (http://mus.well.ox.ac.uk/19genomes/magic.html) algorithms
#'
#' @export 
#' @importFrom utils write.table
#' @param mpcross Object to export 
#' @param cm2bp Factor for converting a genetic map to physical map (how many bp per centiMorgan)
#' @param removeMono Remove monomorphic markers - not required
#' @return Two directories of files ready for input to reconstruction; VARIANT.TABLES (containing founder info) and LINES (containing progeny info). Note that the file reconstruction_lib.c may need modifying if the number of chromosomes simulated is nonstandard.

mpMap2gs <- function(mpcross, cm2bp=250000, removeMono=FALSE) {
  ## create directories to hold parent and final files
  system(paste("mkdir ", dir, "VARIANT.TABLES", sep=""))
  system(paste("mkdir ", dir, "LINES", sep=""))

  if (is.null(mpcross$map)) stop("Need a map to output files\n")
  nchr <- length(mpcross$map)

  nall <- apply(mpcross$founders, 2, function(x) length(table(x)))

  if (removeMono)
	mpcross <- subset(mpcross, markers=which(nall>1))

  ## need to use some consistent coding for genotypes - otherwise how to
  ## know how to recode? 
  ## assume no heterozygotes in founders
  alleles <- c("A", "C", "G", "T")
  founders <- mpcross$founders
  if (length(table(as.vector(founders)))>4) stop("Can only deal with four different alleles in founders currently\n")
  nam <- names(table(as.vector(founders)))
  for (i in 1:length(nam)) founders[founders==nam[i]] <- alleles[i]
  founders <- as.matrix(founders)
           
  nfounders <- nrow(founders)
  nmrk <- vector(length=nchr)
  for (i in 1:nchr) {
      nmrk[i] <- length(mpcross$map[[i]])
      mat <- data.frame(chr=rep(i,nmrk[i]), pse.bp=round(mpcross$map[[i]]*cm2bp), bp=round(mpcross$map[[i]]*cm2bp), nalleles=apply(founders, 2, function(x) length(table(x))), maf=apply(founders, 2, function(x) max(table(x))))
      fou <- founders[, match(names(mpcross$map[[i]]), colnames(mpcross$founders))]
      mat <- cbind(mat, t(fou))
      write.table(mat, file=file(paste(dir, "VARIANT.TABLES/chr", i, ".alleles.txt", sep=""), "wb"), row.names=F, quote=F, sep="\t")
  }
 
  ## recode to match founders - SNPs where they don't match properly (i think)
  ## just get discarded.
  finals <- mpcross$finals 
  for (i in 1:length(nam)) finals[finals==nam[i]] <- alleles[i]
  finals <- as.matrix(finals)
           
  nlines <- nrow(finals)

  for (i in 1:nlines) {
     mat <- data.frame(chr=paste("Chr", rep(1:nchr, nmrk), sep=""), bp=round(mpcross$map[[i]]*cm2bp), allele1=as.vector(founders[1,]), allele2=as.vector(finals[i,]))
     write.table(mat, file=file(paste(dir, "LINES/L", i, sep=""), "wb"), row.names=F, quote=F, col.names=F)
  }

  ## if phenotypes exist, also output them in correct format for genome_scan
  if (!is.null(mpcross$pheno)) {
     pheno <- mpcross$pheno[match(rownames(mpcross$finals), rownames(mpcross$pheno)),]
     mat <- data.frame(SUBJECT.NAME=paste("L", 1:nrow(pheno), sep=""), trait=pheno)
     write.table(mat, file=file(paste(dir, "phenoGS.txt", sep=""), "wb"), row.names=F, quote=F)
  }
}


