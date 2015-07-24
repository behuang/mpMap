#' (Composite) Interval Mapping with linear models for QTL detection in multi-parent crosses
#'
#' Interval mapping in multi-parent crosses using multi-stage linear model 
#' approach including ability to use cofactors (CIM)
#'
#' @export
#' @param object Object of class \code{mpcross}
#' @param threshold Significance threshold for QTL p-values
#' @param chr Subset of chromosomes for which to compute QTL profile
#' @param step Step size at which to compute the QTL profile. See \code{\link[mpMap]{mpprob}} for further description of default values
#' @param responsename Response name for testing 
#' @param ncov Number of marker covariates to search for - default is to search for as many as possible using stepAIC (forward/backward selection)
#' @param window Window of cM on each side of markers where we exclude covariates in CIM
#' @param dwindow Window of markers to use for smoothing in QTL detection 
#' @param mrkpos Flag for whether to consider both marker positions and step positions or just steps. Is overridden if step=0
#' @param fixed If input, vector of fixed effects for each individual to be included in model with main effect and interaction with founder probability 
#' @param foundergroups If input, groups of founders for which to cluster parental alleles together at every marker. Currently overrides mrkpos and step arguments. Note that this is not currently working with ncov>0
#' @param ... Additional arguments
#' @return The original input object with additional component QTLresults containing the following elements:
#' \item{pheno}{Input phenotype data}
#' \item{pvalue}{Each component contains estimated p-values at each position on a given chromosome}
#' \item{wald}{Each component contains Wald statistics at each position on a given chromosome}
#' \item{fndrfx}{Each component contains founder effects estimated at each position on a given chromosome}
#' \item{qtl}{Each component contains the position and effects of a detected QTL}
#' \item{fixedmain}{Each component contains wald statistics for main effect of fixed variable (if input) at each position on a given chromosome}
#' \item{fixedintx}{Each component contains wald statistics at each position on a given chromosome for gene x fixed interaction (if input)}
#' \item{fixedintdf}{Each component contains the df for the gene x fixed interaction at each position on a given chromosome}
#' \item{call}{Input arguments to function} 
#' and with attributes describing the number of QTL detected, and the threshold used for detection. Note: Now uses the function findqtl to find all QTL peaks, see \code{\link[mpMap]{findqtl}} for more information. 
#' @details
#' Modified version of mpIM to just fit linear models for QTL mapping. 
#' This has been separated off from the main function as some functionality
#' has been recently implemented just for the linear model fitting - this 
#' allows for rapid development of the linear model functionality while
#' maintaining the mixed model functionality in a stable state more easily. 
#' No argument names have been changed, hence mpIMlm will be called in the
#' same way as the previous version of mpIM. The only requirement is that
#' phenotypes are included in the mpcross object rather than allowing
#' for separate input of a phenotype matrix. 
 
#' Note that no weights are used in this analysis, which may result in a loss of efficiency compared to a single-stage approach. 
#'
#' If fixed is input will add terms to the model to test for a fixed effect of 
#' the input vector (so make sure the class is correct) and for an interaction
#' between the input vector and the founder haplotypes. Note that only a single
#' fixed covariate can currently be included to avoid overparametrization. 
#'
#' If foundergroups is input, then probabilities at each location will be collapsed within the groups of founders
#' in fitting the model.  
#'
#' @seealso \code{\link[mpMap]{plot.mpqtl}}, \code{\link[mpMap]{summary.mpqtl}}, \code{link[mpMap]{fit.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' ## Two-stage simple interval mapping 
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")


mpIMlm <- function(object, threshold=1e-3, chr, step=0, responsename="predmn", ncov=1000, window=10, dwindow=5, mrkpos=TRUE, fixed, foundergroups, ...)
{
  ### Initial setup for all approaches
  lines <- rownames(object$finals)
  n.founders <- nrow(object$founders) 

  if (missing(chr)) chr <- 1:length(object$map) 
  else if (is.character(chr)) chr <- match(chr, names(object$map))

  ## Set this to the default of computation at every marker. 
  if (!missing(foundergroups)) { 
	mrkpos <- TRUE
	step <- 0
  }
 
  if (!(inherits(object, "mpprob") && attr(object$prob, "step")==step 
  && attr(object$prob, "mrkpos")==mrkpos))  
  object <- mpprob(object, program="qtl", step=step, chr=chr, mrkpos=mrkpos)


  if (!missing(fixed)) { 
    ## check whether fixed values can be matched up to the genotyped lines
    if (length(setdiff(names(fixed), rownames(object$finals)))>0) 
	stop("Observations have fixed effects recorded which have not been genotyped. Please check names and remove lines if necessary\n")
    if (length(setdiff(rownames(object$finals), names(fixed)))>0)
	cat("You do not have fixed effects for all individuals, may want to check this \n")
	vec <- vector(length=nrow(object$pheno))
	names(vec) <- rownames(object$finals)
	vec[match(names(fixed), names(vec))] <- fixed
	if (class(fixed)=="factor") vec <- as.factor(vec)
	fixed <- vec
  }
  fmap <- attr(object$prob, "map")

  ## If this argument is not input - all founders form their own groups; i.e. nothing changes. 
  if (missing(foundergroups)) {
	foundergroups <- matrix(rep(1:n.founders, length(unlist(fmap))), nrow=n.founders)
	colnames(foundergroups) <- unlist(lapply(fmap, names)) 
  }
  output <- object

  fixedmain <- list()
  fixedintx <- list()
  fixedintdf <- list()
  wald <- list()
  pval <- list()
  fndrfx <- list()
  se <- list()
  degf <- list()
  cofactors <- NULL

  require(aods3)
  if (!missing(fixed)) 
    output$pheno <- cbind(output$pheno, fixed)

  idname <- "id"
  pheno <- as.data.frame(output$pheno)
  pheno[[idname]] <- rownames(output$pheno)
  if (is.null(rownames(output$pheno))) pheno[[idname]] <- paste("L", 1:nrow(pheno), sep="")

  ## create a vector of responses just for use in finding marker cofactors
  predmn <- pheno[[responsename]][na.exclude(match(lines, pheno[[idname]]))]

  ### Setting up formulas in the case of covariates
  ### In this case need to just use the markers for selection
  ### set up dataframe with just those values
  if (ncov>0) 
  {
      	require(MASS)
      	if (inherits(object, "mpprob") && attr(object$prob, "step")==0)
          mrkobj <- object else	mrkobj <- mpprob(object, step=0, program="qtl", chr=chr)
	mrkgen <- do.call("cbind", mrkobj$prob)

	## recenter before selecting marker covariates
      	mrkgen <- scale(mrkgen, scale=F)

 	#largest possible formula
    	formula <- as.formula(paste("predmn~", paste(paste("mrkgen[,", seq(1, ncol(mrkgen), n.founders), ":", seq(n.founders, ncol(mrkgen), n.founders), "]", sep=""), collapse="+")))
    	mod <- lm(as.formula(paste("predmn~1", sep="")))
    	modc <- stepAIC(mod, scope=formula, steps=ncov)  

    	covar <- attr(modc$terms, "variables")
      
    	ncov <- length(covar)-2
    	rmv <- list(length=ncov)
    	# the marker index at which each chromosome begins
	gmap <- attr(mrkobj$prob, "map")
	ngmap <- unlist(lapply(gmap, names))
    	nchrmrk <- unlist(lapply(gmap, length))
    	csnchrmrk <- c(0, cumsum(nchrmrk))

	## format the covariates which will be included in the model
      	if (length(covar)>2){
          cofactors <- matrix(nrow=length(covar)-2, ncol=3)
  	  for (k in 3:length(covar))
	  {
	      #as.character(covar[[k]][n.founders]) = index into gen3 for this covariate, e.g. 13:16, 21:24
	      #strsplit(as.character(covar[[k]][n.founders]), ":") = split this into two numbers
	      #start = number of the marker that these four covariates represent
	      ###modification, if start is intedend to be the marker index represented by this covariate then it's off by 4, was 7 before. Overflowed
	        start <- (as.numeric(strsplit(as.character(covar[[k]])[4], ":")[[1]][2]))/n.founders
          #which chromosome is this marker on
	        chrn <- which.max(csnchrmrk[csnchrmrk<start])
	        cpos <- match(ngmap[start], names(fmap[[chrn]]))
  	      # rmvpos is the list of positions for which covariate k should be excluded (with the chromosome stored at the front)
	        rmvpos <- c(chrn, which(fmap[[chrn]] <= fmap[[chrn]][cpos]+window & fmap[[chrn]] >= fmap[[chrn]][cpos] - window))
	        rmv[[k-2]] <- rmvpos
      
	      ## process cofactor names for saving later
              cofactors[k-2,] <- c(ngmap[start], chrn, round(fmap[[chrn]][cpos], 3))
	  }
        cofactors <- as.data.frame(cofactors)
        names(cofactors) <- c("MarkerName", "Chr", "Pos")
  	}
   	    
    	#chromosome numbers of the genetic covariates
    	rmvchr <- unlist(lapply(rmv, function(x) return(x[1])))

    	### set up a list with a matrix for each component representing chr 
    	### one row per genetic covariate, one column per location 
    	### FALSE indicates to EXCLUDE the covariate for that location
    	terms <- fmap
    	for (i in 1:length(terms)) 
    	{
	  terms[[i]] <- matrix(data=TRUE, nrow=ncov, ncol=length(terms[[i]]))
	  ## if cofactors are on this chromosome, set array accordingly
	  if (i %in% rmvchr)
	  for(k in 1:length(rmv))
	  if(rmv[[k]][1] == i) 	terms[[i]][k, rmv[[k]][-1]] <- FALSE
    	}

    	## alter fixed terms in model depending on position
    	formall <- modc$call[[2]]
    	termlab <- attr(modc$terms, "term.labels")
    }
    # matrix has all TRUE values with no covariates
    else 
    {
    	terms <- fmap
    	for (i in 1:length(terms)) 
 	  terms[[i]] <- matrix(data=TRUE, nrow=ncov, ncol=length(terms[[i]]))
    	termlab <- ""
    }

  for (j in chr)
  {
    nam <- names(object$map)[j]
    cat("------Analyzing Chr ",nam,"----------\n")
    gen <- object$prob[[nam]]
    fgc <- foundergroups[, match(names(fmap[[j]]), colnames(foundergroups)), drop=F] ## need to make sure of naming scheme

    ## All set to NA by default
    wald[[nam]] <- rep(NA, ncol(gen)/n.founders)
    pval[[nam]] <- rep(NA, ncol(gen)/n.founders)
    degf[[nam]] <- rep(NA, length=ncol(gen)/n.founders)
    fndrfx[[nam]] <- matrix(nrow=n.founders, ncol=ncol(gen)/n.founders)
    se[[nam]] <- matrix(nrow=n.founders, ncol=ncol(gen)/n.founders)
    if (!missing(fixed)) {
	fixedmain[[nam]] <- rep(NA, ncol(gen)/n.founders)
	fixedintx[[nam]] <- rep(NA, ncol(gen)/n.founders)
	fixedintdf[[nam]] <- rep(NA, ncol(gen)/n.founders)
    }
    df <- matrix(nrow=nrow(pheno), ncol=ncol(gen))
  
    colnames(gen) <- paste("P", rep(1:(ncol(gen)/n.founders), each=n.founders), "G", LETTERS[1:n.founders], sep="")

    genid <- vector()
    for (k in 1:nrow(pheno)) 
      genid[k] <- match(as.character(pheno[[idname]][k]), as.character(lines))
    
#    genid <- vector()
#    pheid <- vector()
#    for (k in 1:length(lines)) {
#     matchid <- which(as.character(pheno[[idname]])==as.character(lines[k]))
#     genid <- c(genid, rep(k, length(matchid)))
#     pheid <- c(pheid, matchid) 
#    }
    df <- as.matrix(gen[genid,])
    df <- cbind(pheno, df)
    df <- as.data.frame(df)
##    df <- df[match(lines, df$id),] ## this is causing problems for multivar analysis
    ## need to check whether it's required for some other analysis. 
    
    for (k in 1:ncol(pheno))
    {
     fx <- paste("as.", class(pheno[,k]), sep="")
     df[,k] <- do.call(fx, list(df[,k]))
    }
    names(df) <- c(names(pheno), colnames(gen))

    # recenter probabilities to 0
    df[, (ncol(pheno)+1):ncol(df)] <- scale(df[, (ncol(pheno)+1):ncol(df)], scale=FALSE)

    for (k in seq(1, ncol(gen), n.founders))
    {
	index <- (k-1)/n.founders+1
	df2 <- df[, 1:ncol(pheno)]

	### Need to update the probabilities so that founders within groups are summed
	### Should be able to do this in some way with matrix multiplication

	mat <- as.matrix(df[, ncol(pheno)+k:(k+n.founders-1)])
	mm <- model.matrix(~factor(fgc[, index])-1)
 	mat <- mat %*% mm
	if (ncol(mat)<19) {
	varc <- apply(mat, 2, function(x) var(x, na.rm=T))
	if (any(varc<1e-6)) {
	  mat <- mat[, -(which(varc<1e-6))]
	  mat <- t(apply(mat, 1, function(x) x/sum(x, na.rm=T)))
	}
 	}
	colnames(mat) <- paste("P", index, "G", LETTERS[1:ncol(mat)], sep="")
	df2 <- cbind(df2, mat)
	ngrps <- ncol(mat)

	if (ncov>0) 
	# include all necessary covariates
	form <- as.formula(paste("predmn~", paste(termlab[which(terms[[j]][,index]==1)], collapse="+"), "+", paste(names(df)[grep(paste("P", index, "G", sep=""), names(df))], collapse="+"), sep="")) 
	else form <- as.formula(paste("predmn~", paste(names(df2)[(ncol(pheno)+1):ncol(df2)], collapse="+"), sep=""))
	  # fit the model
	
	if (!missing(fixed)) 
	  form <- as.formula(paste("predmn~", paste(paste("fixed*", names(df2)[(ncol(pheno)+1):ncol(df2)], sep=""),collapse="+"), sep="")) 
	    
    	### Note: unlikely to need this when founders are compressed
    	qq <- which(cor(df2[, (ncol(pheno)+1):ncol(df2)])>.95, arr.ind=T)
    	qq <- qq[qq[,1]<qq[,2],, drop=F]
    	if (nrow(qq)>0) {
      	  for (pp in 1:nrow(qq))
          df2[, ncol(pheno)+qq[pp,1]] <- df2[,ncol(pheno)+qq[pp,2]] <- (df2[,ncol(pheno)+qq[pp,1]]+df2[,ncol(pheno)+qq[pp,2]])/2
    	}
    
	mod <- lm(form, data=df2)
    
	# Collinearity means some coefficients will be NA
	coe <- coef(mod)[which(!is.na(coef(mod)))]	

	# if we did not get any non NA values, window was not large enough
	if(length(grep(paste("P", index, "G", sep=""), names(coe))) != 0)
	{
	    if (!missing(fixed)) {
	    	# Need to do terms separately
	    	index1 <- grep(":P", names(coe))
	    	wt <- wald.test(varb=vcov(mod), b=coe, Terms=index1)
	    	fixedintx[[nam]][index] <- wt$result$chi2[1]
		fixedintdf[[nam]][index] <- wt$result$chi2[2]
	   	
		## now do main effect
		index2 <- grep("fixed", names(coe))
		index2 <- setdiff(index2, index1)
		wt <- wald.test(varb=vcov(mod), b=coe, Terms=index2)
		fixedmain[[nam]][index] <- wt$result$chi2[1]
	    } else index1 <- NULL
		
	    index3 <- grep(paste("P", index, "G", sep=""), names(coe))
	    index3 <- setdiff(index3, index1) 

	    #test the current location for significance (actually a joint test)
	    wt <- wald.test(varb=vcov(mod), b=coe, Terms=index3)

	    pval[[nam]][index] <- wt$result$chi2[3] 
	    wald[[nam]][index] <- wt$result$chi2[1]
	    degf[[nam]][index] <- wt$result$chi2[2]

	    a <- summary(mod)$coefficients
	    index4 <- grep(paste("P", index, "G", sep=""), rownames(a))
	    index4 <- setdiff(index4, grep(":P", rownames(a)))
	    fndrfx[[nam]][,index] <- c(a[index4,1], rep(NA, n.founders-length(index4))) 
	    se[[nam]][,index] <- c(a[index4,2], rep(NA, n.founders-length(index4))) 
	  }
	  else
	    cat("Error when testing location ", index, " along chromosome ", j, ". Window may be too small.\n ")
     }
  ind <- which(pval[[nam]]==0)
  if (length(ind)>0)
  for (m in ind) pval[[nam]][m] <- -pchisq(wald[[nam]][m], degf[[nam]][m], log.p=TRUE)
  }
 
  minp <- unlist(lapply(pval, function(x) min(x, na.rm=TRUE)))
  maxw <- unlist(lapply(wald, which.max))
  sigchr <- which(minp < threshold)
  map <- attr(object$prob, "map")

  ## select out most significant QTL
  pos <- vector(length=length(fndrfx))
  pos[1:length(pos)] <- NA
  pos[sigchr] <- maxw[sigchr]

  ## positions in cM
  posqtl <- vector(length=length(pos))
  for (i in 1:length(pos))
    posqtl[i] <- (!is.na(pos[i]))*map[[i]][maxw[i]]

  qtl <- list()
  for (i in which(!is.na(pos)))
  {
    qtl[[names(map)[i]]] <- t(as.matrix(c(posqtl[i], fndrfx[[names(map)[i]]][, pos[i]], se[[names(map)[i]]][, pos[i]])))
    attr(qtl[[names(map)[i]]], "index") <- pos[i]
  }
  
  results <- list()
  if (!missing(baseModel))   results$baseModel <- baseModel
  else results$baseModel <- lm(predmn~1, data=pheno)
  results$pheno <- pheno
  results$pvalue <- pval
  results$wald <- wald
  results$se <- se
  results$degf <- degf
  results$fndrfx <- fndrfx
  results$qtl <- qtl
  results$call <- match.call()
  attr(results$qtl, "threshold") <- threshold
#  attr(results$qtl, "nqtl") <- sum(!is.na(pos)) 
  attr(results$qtl, "ncov") <- ncov
  attr(results$qtl, "window") <- window
  attr(results, "method") <- method
  if (!missing(fixed)) {
	results$fixedmain <- fixedmain
	results$fixedintx <- fixedintx
	results$fixedintdf <- fixedintdf
  }

  output$QTLresults <- results
  output$QTLresults$cofactors <- cofactors
  class(output) <- c("mpqtl", class(object))

  output <- findqtl(output, dwindow, threshold=-log10(threshold))
  output
}
