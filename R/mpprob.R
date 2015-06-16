#' Compute founder probabilities for multi-parent crosses
#'
#' Using haplotype probabilities, computes the probability that each location on a genome was inherited from each founder. Locations are run either at markers only, at the midpoints of all intervals or at step sizes of x cM. Probabilities can be computed using internally, or with R/happy.hbrem or R/qtl.
#' @export mpprob
#' @aliases mpprob calcmpprob print.mpprob
#' @param object Object of class \code{mpcross}
#' @param chr Subset of chromosomes
#' @param step Step size (in cM) to create grid of positions at which to compute probabilities. At default value of 0, probabilities are calculated at marker positions only
#' @param impmap Full set of positions at which to compute probabilities. Overrides all other parameters for determining positions
#' @param threshold Threshold for calling founder probabilities
#' @param mapfx Map function used to convert map to recombination fractions
#' @param ibd Flag to indicate whether to compute probabilities using IBD genotypes
#' @param mrkpos Flag to indicate whether to compute probabilities at both marker positions and step size or just step size. Is overridden for step size of 0. Not used for progam = "mpMIM"
#' @param program R package to use to compute probabilities
#' @param tempfiledirectory Directory in which to output temporary files. Default is current working directory
#' @param generations Number of generations to assume in HAPPY. see \code{happy}
#' @param est Flag for whether to impute founder alleles
#' @param geprob Probability of genotyping error
#' @details If \code{program=="mpMap"} then probabilities are computed using flanking markers at positions across the genome and represent 3-point haplotype probabilities. If \code{program=="happy"} then probabilities are computed using default values in R/happy.hbrem, which calculates ancestral haplotypes without using pedigree information. This only allows for probabilities to be computed at midpoints of intervals. If \code{program=="qtl"} then probabilities are computed from multipoint founder probabilities in R/qtl. 
#' 
#' If \code{step<0} for R/mpMap or R/qtl or \code{step==0} for R/happy.hbrem, then probabilities are computed at the midpoints of marker intervals. However, if \code{step==0} for R/qtl or R/mpMap, probabilities are computed only at marker locations.
#' @return The input mpcross object is returned with two additional components:
#' \item{prob}{A list with founder probabilities for each chromosome. Format is a matrix with n.founders * n.markers columns and n.lines rows. Each group of n.founders columns will add up to 1. Founder probabilities are in the order of founders in the input founder matrix. }
#' \item{estfnd}{A list with estimated founders for each chromosome. Format is a matrix with n.markers columns and n.lines rows. Missing values indicate where no founder probability exceeded the input threshold. Numeric values for founders indicate the row in the input founder matrix corresponding to the estimated founder.}
#' @seealso \code{\link[mpMap]{plot.mpprob}}, \code{\link[mpMap]{summary.mpprob}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' plot(mpp.dat)
#' summary(mpp.dat)

mpprob <- function(object, chr, step=0, impmap, mrkpos=TRUE, mapfx=c("haldane", "kosambi"), ibd=FALSE, threshold=0.7, program=c("mpMap", "qtl", "happy", "mpMIM"), tempfiledirectory="", generations=5, est=TRUE, geprob=0.0001)
{
	program <- match.arg(program)
	if (missing(mapfx)) mapfx="haldane"

  if (missing(chr))
  {
	chr <- names(object$map)
	print("No chromosomes specified, will default to all")
  }

  if (is.numeric(chr)) chr <- names(object$map)[chr]
  if (step<0) mrkpos <- FALSE ## override for midpoints with mpMap
  if (step==0) mrkpos <- TRUE ## override mrkpos flag in this case

  n.founders <- nrow(object$founders)
  n.finals <- nrow(object$finals)
  out <- object
 
  out$prob <- list()
  out$estfnd <- list()
  attr(out$estfnd, "threshold") <- threshold
  ## as a default, convert the probabilities to estimated founders as well.
  ffx <- function(x) 
	if (is.na(max(x)) | max(x) <= threshold) return(NA) else return(which.max(x))
 
  chr1 <- vector()
  for (i in chr) 
	if (length(object$map[[i]])==1) chr1 <- c(chr1, i)
  chr <- setdiff(chr, chr1)

  if (length(chr1)>0) 
	mrkpos <- TRUE

  if (program=="mpMap")
    out$prob <- calcmpprob(object, chr, step, mapfx, ibd, mrkpos)
  if (program=="qtl")
  {
    if(!("Design" %in% names(object$pedigree)))
    {
    if (!inherits(object$pedigree, "data.frame")) object$pedigree <- as.data.frame(object$pedigree)
		object$pedigree <- cbind(object$pedigree, identifyDesign(object$pedigree), stringsAsFactors = FALSE)
		colnames(object$pedigree)[ncol(object$pedigree)] <- "Design"
    }
	removeSelfing <- function(designName)
	{
		searchResult <- regexpr("self", designName, fixed=T)
		if(searchResult != -1)
		{
			return(substr(designName, 1, searchResult-1))
		}
		return(designName)
	}
	designColumn <- object$pedigree[, "Design"]
	observedDesigns <- unique(sapply(designColumn[object$pedigree[, "Observed"] == 1], removeSelfing))
	if(length(observedDesigns) != 1)	
	{
		stop("Only one type of design can currently be used with program qtl")
	}
	if(observedDesigns == "8wayG3aic1")
	{
		qtlType <- "ri8selfIRIP1"
	}
	else if(observedDesigns == "8wayG3aic2")
	{
		qtlType <- "ri8selfIRIP2"
	}
	else if(observedDesigns == "8wayG3aic3") 
	{
		qtlType <- "ri8selfIRIP2"
	}
	else if(observedDesigns == "8wayG3aic10")
	{
		qtlType <- "ri8selfIRIP10"
	}
	else if(observedDesigns == "8wayG3")
	{
		qtlType <- "ri8self"
	}
	else if(observedDesigns == "4wayG2")
	{
		qtlType <- "ri4self"
	}
	else
	{
		stop("Unable to convert design to one known by R/qtl")
	}

	if (tempfiledirectory != "") write2cross(object, filestem=paste(tempfiledirectory, "/tmp", sep=""))
	else write2cross(object, filestem="tmp")
   	cr <- qtl:::readMWril(tempfiledirectory, "tmp.ril.csv", "tmp.founder.csv", type=qtlType)
 	cr <- subset(cr, chr=chr)
	if (!missing(impmap)) 
	gp <- calc.genoprob2(cr, pos=impmap, error.prob=geprob)
	else if (step >= 0 & mrkpos==F)
	{
	  impmap <- list()
	  for (i in chr)  {
		impmap[[i]] <- c(seq(0, max(cr$map[[i]]), step), max(cr$map[[i]]))
	 	names(impmap[[i]]) <- c(paste("loc", impmap[[i]][1:(length(impmap[[i]])-1)],sep=""), names(cr$map[[i]])[length(cr$map[[i]])])
	  } 
	  gp <- calc.genoprob2(cr, pos=impmap, error.prob=geprob)
   	}
	else if (step >= 0)
	gp <- calc.genoprob(cr, step=step, error.prob=geprob)
	else if (step < 0) {
		pos <- list()
		for (i in chr)
		{ 
		  pos[[i]] <- vector(length=length(object$map[[i]])*2-1)
		  pos[[i]][seq(1, length(pos[[i]]), 2)] <- object$map[[i]]
		  pos[[i]][seq(2, length(pos[[i]]), 2)] <- object$map[[i]][1:(length(object$map[[i]])-1)]+diff(object$map[[i]])/2
		  names(pos[[i]])[seq(1, length(pos[[i]]), 2)] <- names(object$map[[i]])
		  names(pos[[i]])[seq(2, length(pos[[i]]), 2)] <- paste("loc", round(pos[[i]][seq(2, length(pos[[i]]), 2)], 1), sep="")
		} 
		class(pos) <- "map"
	gp <- calc.genoprob2(cr, pos=pos, error.prob=geprob)
	}
	if (length(chr1)>0)
	for (i in chr1) {
	m <- match(names(object$map[[i]]), colnames(object$finals))
	gp$geno[[i]]$prob <- array(dim=c(nrow(object$finals), 1, nrow(object$founders)))
	gp$geno[[i]]$prob[,1,] <- t(sapply(object$finals[,m], function(x) return(1*(x==object$founders[,m])/(sum(object$founders[,m]==x)))))
	attr(gp$geno[[i]]$prob, "map") <- object$map[[i]]
	}

	prob1 <- lapply(gp$geno, function(x) return(x$prob))
	prob <- lapply(prob1, function(x) {
		mat <- matrix(nrow=dim(x)[1], ncol=dim(x)[2]*dim(x)[3])
		for (i in 1:dim(x)[3])
			mat[, seq(i, ncol(mat), dim(x)[3])] <- as.matrix(x[,,i])
		return(mat)})
	crmap <- lapply(gp$geno, function(x) return(attr(x$prob, "map")))

	for (i in chr)
	{
		#Break into markers and step locations
		isMarker <- which(names(crmap[[i]]) %in% names(object$map[[i]]))
		markerNames <- names(object$map[[i]])[names(object$map[[i]]) %in% names(crmap[[i]])]
		
		if(step > 0) 
		{
			stepSequence <- seq(min(object$map[[i]]), max(object$map[[i]]), step)
			#Numerical inaccuracy means this search needs some numeric tolerance
			isStep <- which(sapply(crmap[[i]], function(x) min(abs(x - stepSequence)) < 1e-6))
		}
		else isStep <- c(1:length(crmap[[i]]))
		#for simulated data we could have a marker at the same location as a step - position. In this case count it as a marker
		isStep <- setdiff(isStep, isMarker)
		stepNames <- names(crmap[[i]])[isStep]
		#Everything has to be either a step location or a marker location
		if(length(c(isStep, isMarker)) != dim(prob[[i]])[2]/n.founders) stop("Internal error")
		
		colnames(prob[[i]]) <- 1:ncol(prob[[i]])
		colnames(prob[[i]])[(rep(isStep, each=n.founders) - 1) * n.founders + rep(1:n.founders, length(isStep))] <- paste(rep(stepNames, each=n.founders), ", Founder ", rep(1:n.founders, times=length(isStep)), sep="")
		colnames(prob[[i]])[(rep(isMarker, each=n.founders) - 1) * n.founders + rep(1:n.founders, length(isMarker))] <- paste(rep(markerNames, each=n.founders), ", Founder ", rep(1:n.founders, times=length(isMarker)), sep="")
		rownames(prob[[i]]) <- rownames(object$finals)
	}
	#If there aren't meant to be any markers, remove them? Not clear why this is here, unless the R/qtl code always generates data at marker locations and it has to be removed later (as in here). 
	if (length(isMarker)>0 && mrkpos == FALSE)
	{
	  m <- match(names(object$map[[i]]), colnames(object$finals))
    m2 <- (rep(isMarker, each=n.founders)-1)*n.founders+rep(1:n.founders, length(m))
		crmap[[i]] <- crmap[[i]][-isMarker]
		prob[[i]] <- prob[[i]][,-m2]
	}
	attr(prob, "map") <- crmap

	out$prob <- prob
	names(out$prob) <- names(attr(out$prob,"map")) <- c(chr, chr1)
  }

  if (program=="happy"){
	require(happy.hbrem)
	prob <- list()
	map <- list()
	write2happy(object, filestem="tmp")
 	hin <- happy.hbrem::happy("tmp.data", "tmp.alleles", generations=generations, haploid=TRUE)
	# number of marker intervals
	nint <- unlist(lapply(object$map, length))-1
	cnint <- c(0, cumsum(nint))
	mrkint <- 1
	for (i in 1:length(chr)) {
	  prob[[i]] <- matrix(nrow=nrow(object$finals), ncol=nint[i]*n.founders)
	  while(mrkint <= cnint[i+1]) {
	    prob[[i]][, (mrkint-cnint[i]-1)*n.founders+1:n.founders] <- happy.hbrem::hdesign(hin, mrkint)
	    mrkint <- mrkint+1
	  }
  	  map[[i]] <- object$map[[i]][1:(length(object$map[[i]])-1)] + diff(object$map[[i]])/2
	  colnames(prob[[i]]) <- paste("C",i,"P ", rep(1:nint[i], each=n.founders), ", Founder ", 1:n.founders, sep="")
	  rownames(prob[[i]]) <- rownames(object$finals)
	} # end of chr loop
    class(map) <- "map"
    attr(prob, "map") <- map
    out$prob <- prob
	names(out$prob) <- names(attr(out$prob,"map")) <- c(chr, chr1)
  } # end of happy loop
	if(program == "mpMIM")
	{
		object <- subset(object, chr=chr)
		require(mpMIM)
		data.file <- tempfile()
		probabilities.file <- tempfile()
		locations.file <- tempfile()
		#We don't care about the phenotypic data, this is only for probabilities anyway
		if(!("pheno" %in% names(object)))
		{
			warning("Creating dummy phenotype");
			object$pheno <- data.frame(rep(0, nrow(object$finals)))
		}
		mim.export(object, outputfile = data.file, phenocolumn = 1)
		if(step == 0)
		{
			mmGenerateLocations(dataFile = data.file, locationsFile = locations.file, perinterval = 0)
		}
		else
		{
			mmGenerateLocations(dataFile = data.file, locationsFile = locations.file, equalspacing = step)
		}
		generatedMap <- mim.import(locations.file)
		mmCalculateProbabilities(dataFile = data.file, locationsFile = locations.file, probabilitiesFile = probabilities.file)
		probabilities <- mim.import(probabilities.file)
		restructuredProbabilities <- list()
		for(currentChromosome in names(generatedMap)) 
		{
			currentChromosomeProbabilities <- probabilities[names(generatedMap[[currentChromosome]])]
			restructuredProbabilities[[currentChromosome]] <- do.call(cbind, as.list(currentChromosomeProbabilities))
			colnames(restructuredProbabilities[[currentChromosome]]) <- paste(rep(names(currentChromosomeProbabilities), each = n.founders), ", Founder ", 1:n.founders, sep="")
		}
		unlink(data.file)
		unlink(probabilities.file)
		unlink(locations.file)
		out$prob <- restructuredProbabilities
		attr(out$prob, "map") <- generatedMap
	}
	
	if (est) 
	{
		for (ii in names(attr(out$prob, "map")))
		{
			haps <- out$prob[[ii]]    
			haps[is.nan(haps)] <- NA

			fmat <- matrix(nrow=nrow(haps), ncol=ncol(haps)/n.founders)
			for (kk in seq(1, ncol(haps), n.founders))
			{
				fmat[,(kk-1)/n.founders+1] <- apply(haps[,kk:(kk+n.founders-1)], 1, ffx)
				fmat[,(kk-1)/n.founders+1] <- factor(fmat[,(kk-1)/n.founders+1], levels=1:n.founders)
			}
			out$estfnd[[ii]] <- fmat
			colnames(out$estfnd[[ii]]) <- names(attr(out$prob, "map")[[ii]])
			rownames(out$estfnd[[ii]]) <- rownames(out$prob[[ii]])
		}
	}
  
  attr(out$prob, "step") <- step
  attr(out$prob, "program") <- program
  attr(out$prob, "mapfx") <- mapfx
  attr(out$prob, "mrkpos") <- mrkpos

	#name all the rows of probabilities  
	for(i in names(out$prob))
	 	rownames(out$prob[[i]]) = paste("L", 1:n.finals, sep="")
	
  class(out) <- unique(c("mpprob", class(object)))
  return(out)
}
