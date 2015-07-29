#' Estimate pairwise recombination fractions between markers
#'
#' Estimates pairwise recombination fractions by maximizing the likelihood for a multi-parent cross over a grid of possible values. Theta values and corresponding LOD scores are returned for each pair of markers in the object.
#' @export
#' @import Rcpp
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @param r Grid of potential recombination values. If missing the function will maximize over (0, .005, .01, .015, ... , .095, .1, .11, .12, ... .49, .5). 
#' @param gpu Boolean value, true indicates that a GPU should be used if available
#' @param lineWeights In some cases of segregation distortion it can be useful to weight the contribution of each line to the likelihood
#' @param mpi Flag for whether to parallelize the computation
#' @param \dots Additional arguments to be passed on to mpestrfMpi
#' @return Returned object is of the class 'mpcross' with the additional component \code{rf}. If n.mrk is the number of markers genotypes, this is a list with components:
#' \item{rf$theta}{ n.mrk x n.mrk matrix of estimated recombination fractions between each pair of loci}
#' \item{rf$lod}{ n.mrk x n.mrk matrix of LOD scores at the estimated recombination values}
#' \item{rf$lkhd}{ n.mrk x n.mrk matrix of likelihood values at the estimated recombination values}
#' @seealso \code{\link[mpMap]{mpcross}}, \code{\link[mpMap]{plot.mpcross}}
#' @examples
#' map <- qtl::sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, 
#'	  qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), 
#'	  seed=1)
#' dat.rf <- mpestrf(sim.dat)

mpestrf <- function(object, r, gpu, lineWeights, mpi=FALSE, ...)
{
	if(mpi) 
	{
		if (!requireNamespace("Rmpi", quietly = TRUE)) 
    		stop("Rmpi needed for MPI mpestrf to work. Please install it.\n",
      		call. = FALSE)
		if (Rmpi::mpi.comm.size()>0) 
		{
			if(inherits(object, "mpcross")) 
			{
				tryCatch({
				return(mpestrfMpi(list(object), r, gpu, lineWeights, ...))
				}, error = function(err) {
				stop(paste("mpestrfMpi failed: ", err))
				})
			}
			else
			{
				return(mpestrfMpi(object, r, gpu, lineWeights, ...))
			}
		}
		else stop("Attempted to use MPI outside mpirun")
	} 
	else 
	{
		if(inherits(object, "mpcross"))
		{
			if(missing(lineWeights))
			{
				lineWeights <- rep(1, nrow(object$finals))
			}
			if(!is.numeric(lineWeights) || length(lineWeights) != nrow(object$finals)) stop("Invalid input for argument lineWeights")
			return(mpestrfSubset(objects = list(object), r=r, gpu=gpu, lineWeights = list(lineWeights)))
		}
		else
		{
			if(missing(lineWeights))
			{
				lineWeights <- lapply(object, function(x) rep(1, nrow(x$finals)))
			}
			if(length(lineWeights) != length(object)) stop("Invalid input for argument lineWeights")
			for(i in 1:length(object))
			{
				if(!is.numeric(lineWeights[[i]]) || length(lineWeights[[i]]) != nrow(object[[i]]$finals)) stop("Invalid input for argument lineWeights")
			}
			return(mpestrfSubset(objects = object, r=r, gpu=gpu, lineWeights = lineWeights))
		}
	}
}

stopifany <- function(...) { stopifnot(!any(...)) }

# hack which I'm yet to understand
custom.bcast.Robj2slave <- function(object) {
	tryCatch({
	customrecv <- function() {
		  customobjects <- NULL
	  	  Rmpi::mpi.send.Robj(0,0,1)
		  customobjects <<- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag())
	}
	Rmpi::mpi.bcast.Robj2slave(customrecv)
	Rmpi::mpi.bcast.cmd(customrecv())
	closed_slaves <- 0 
	n_slaves <- Rmpi::mpi.comm.size()-1 
	while (closed_slaves < n_slaves) {		
			message <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag()) 
			message_info <- Rmpi::mpi.get.sourcetag() 
			slave_id <- message_info[1] 
			tag <- message_info[2] 
			Rmpi::mpi.send.Robj(object, slave_id, 1); 
			closed_slaves <- closed_slaves + 1 
	}
	}, error = function(err) {
	   stop(paste("custom.bcast.Robj2slave failed: ", err))
	})
}

masterStoreTile <- function(rect,rf, theta, lod, lkhd)
{
	theta[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$theta[1:rect$sizex,1:rect$sizey]
	lod[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$lod[1:rect$sizex,1:rect$sizey]
	lkhd[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$lkhd[1:rect$sizex,1:rect$sizey]
}

mpi.run.slavempestrf <- function(gridDimX, gridDimY, theta, lod, lkhd)
{
	print("In mpi.run.slavempestrft")
	tryCatch({
	Rmpi::mpi.bcast.cmd(slavempestrf())
	print("broadcast done")
	closed_tiles <- 0
        n_tiles <- gridDimX * gridDimY
        while (closed_tiles < n_tiles) {
	      		print("mpi.run.slavempestrft waiting on a nibble..")
	                tileDim <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag())
                        message_info <- Rmpi::mpi.get.sourcetag()
                        slave_id <- message_info[1]
                        tag <- message_info[2]
	      		print(paste("mpi.run.slavempestrft waiting on a tile from ",slave_id))
			res <- Rmpi::mpi.recv.Robj(slave_id, tag)
			masterStoreTile(tileDim, res$rf, theta, lod, lkhd)
                        closed_tiles <- closed_tiles + 1
			print(paste(closed_tiles," of ", n_tiles, " complete"))
	}
	}, error = function(err) {
	   stop(paste("mpi.run.slavempestrf failed: ", err))
	})
	return(res)
}

mpestrfMpi <- function(objects, r, gpu, lineWeights,leaveAsFileBacked=FALSE, onlyMasterWrites=TRUE, dir_base, passObjectsAsFile = FALSE)
{
	nmrks <- ncol(objects[[1]]$founders)
	if(missing(dir_base)) dir_base <- "bigdata/"
	# create the empty output matrices
	if (!requireNamespace("bigmemory", quietly = TRUE)) 
    	  stop("bigmemory needed for mpestrfMPI to work. Please install it.\n",
      	  call. = FALSE)

	# decide how we want to decompose the data
        gridDimX <- 1
      	gridDimY <- Rmpi::mpi.comm.size() - 1

      tryCatch({
	dir.create(dir_base, showWarnings=FALSE)
	file_base <- Sys.getenv("PBS_JOBID")
	if (file_base == "") 
	{
		if (!requireNamespace("R.utils", quietly = TRUE)) 
    		stop("R.utils needed for mpestrfMpi to work. Please install it.\n",
      		call. = FALSE)
		file_base <- paste(R.utils::System$getHostname(),Sys.getpid(),sep="-")
	}
	thetabf=paste(file_base,"theta.bin",sep=".")
	thetadf=paste(file_base,"theta.desc",sep=".")
	lodbf=paste(file_base,"lod.bin",sep=".")
	loddf=paste(file_base,"lod.desc",sep=".")
	lkhdbf=paste(file_base,"lkhd.bin",sep=".")
	lkhddf=paste(file_base,"lkhd.desc",sep=".")
	theta <- bigmemory::filebacked.big.matrix(nmrks, nmrks, backingpath=dir_base, backingfile=thetabf, type="double", descriptorfile=thetadf)
	lod <- bigmemory::filebacked.big.matrix(nmrks, nmrks,  backingpath=dir_base, backingfile=lodbf, type="double", descriptorfile=loddf)
	lkhd <- bigmemory::filebacked.big.matrix(nmrks, nmrks,  backingpath=dir_base, backingfile=lkhdbf, type="double", descriptorfile=lkhddf)
	thetadesc <- bigmemory::describe(theta)
	loddesc <- bigmemory::describe(lod)
	lkhddesc <- bigmemory::describe(lkhd)
	}, error = function(err) {
	   stop(paste("mpestrfMpi failed to setup big matrices: ", err))
	})

	# upload required data to the slaves
	tryCatch({
		if (passObjectsAsFile) 
		{
			filename <- paste(file_base,'mpcrossobjects',sep=".")
			save(objects,file=paste(dir_base,filename,sep="/"))
		} 
		else 
		{
			custom.bcast.Robj2slave(objects)
		}

	if (!missing(r)) Rmpi::mpi.bcast.Robj2slave(r)
	Rmpi::mpi.bcast.Robj2slave(gpu)
	if (!missing(lineWeights)) Rmpi::mpi.bcast.Robj2slave(lineWeights)
	Rmpi::mpi.bcast.Robj2slave(dir_base)
	Rmpi::mpi.bcast.Robj2slave(file_base)
	Rmpi::mpi.bcast.Robj2slave(thetadesc)
	Rmpi::mpi.bcast.Robj2slave(loddesc)
	Rmpi::mpi.bcast.Robj2slave(lkhddesc)
	Rmpi::mpi.bcast.Robj2slave(slavempestrf)
	Rmpi::mpi.bcast.Robj2slave(nmrks)
	Rmpi::mpi.bcast.Robj2slave(gridDimX)
	Rmpi::mpi.bcast.Robj2slave(gridDimY)
	Rmpi::mpi.bcast.Robj2slave(onlyMasterWrites)
	Rmpi::mpi.bcast.Robj2slave(passObjectsAsFile)

	}, error = function(err) {
	   stop(paste("mpestrfMPI failed to bcast data to slaves: ", err))
	})
	
	# process the data
	tryCatch({

	if (onlyMasterWrites) {
	# use the method where slaves send the tile back to the master to write out
        res <- mpi.run.slavempestrf(gridDimX, gridDimY, theta, lod, lkhd)

	} else {
	# use the method where slaves write their own tiles to disk
  	reslist <- Rmpi::mpi.remote.exec(slavempestrf())
	if (all(unlist(lapply(reslist,is.atomic)))) {
	   # one or more slaves returned an error message
	   stop(paste("Slave failed:",reslist))
	}
	res <- reslist$slave1$result

	}
	}, error = function(err) {
	   stop(paste("A slavempestrft call failed: ", err))
	})

	if (is.list(res) == FALSE) {
	  stop(paste("slavempestrf failed",reslist))
        }

	# extract the results
	# it seems we need to reattach if the slaves have made changes..
	theta <- bigmemory::attach.big.matrix(thetadesc,backingpath=dir_base)
	lkhd <- bigmemory::attach.big.matrix(lkhddesc,backingpath=dir_base)
	lod <- bigmemory::attach.big.matrix(loddesc,backingpath=dir_base)
	if (leaveAsFileBacked) {
	       res$rf$thetadesc <- thetadesc
	       res$rf$loddesc <- loddesc
	       res$rf$lkhddesc <- lkhddesc
	       res$rf$theta <- theta
	       res$rf$lkhd <- lkhd
	       res$rf$lod <- lod
	} else {
	      tryCatch({
	       # may not fit into memory...
	       res$rf$theta <- theta[,]
	       res$rf$lkhd <- lkhd[,]
	       res$rf$lod <- lod[,]
	      }, error = function(err) {
		 stop(paste("mpestrfMPI failed, consider leaveAsFileBacked. Error: ", err))
	      })
	}
	stopifany(is.null(res$rf$theta), is.null(res$rf$lkhd), is.null(res$rf$lod))

	# regenerate the dimnames
	markerNames <- colnames(res$finals)
	dimnames(res$rf$theta) <- list(markerNames,markerNames)
	dimnames(res$rf$lod) <- list(markerNames,markerNames)
	dimnames(res$rf$lkhd) <- list(markerNames,markerNames)

	return(res)
}

slavempestrf <- function() {

	if (missing(passObjectsAsFile)) passObjectsAsFile <- NULL
	if (missing(gpu)) gpu <- NULL
	if (missing(customobjects)) customobjects <- NULL
	if (missing(lineWeights)) lineWeights <- NULL
	if (missing(dir_base)) dir_base <- NULL
	if (missing(file_base)) file_base <- NULL
	if (missing(thetadesc)) thetadesc <- NULL
	if (missing(loddesc)) loddesc <- NULL
	if (missing(lkhddesc)) lkhddesc <- NULL
	if (missing(nmrks)) nmrks <- NULL
	if (missing(gridDimX)) gridDimX <- NULL
	if (missing(gridDimY)) gridDimY <- NULL
	if (missing(onlyMasterWrites)) onlyMasterWrites <- NULL

	if (!requireNamespace("bigmemory", quietly = TRUE)) 
    	  stop("bigmemory needed for mpestrfMPI to work. Please install it.\n",
      	  call. = FALSE)
      slaves <- Rmpi::mpi.comm.size() - 1
      myID <- Rmpi::mpi.comm.rank()

      	tryCatch(
		{
			if (passObjectsAsFile) 
			{
				filename <- paste(file_base,'mpcrossobjects',sep=".")
				objects <- load(file=paste(dir_base,filename,sep="/"))
			}
			else 
			{
				objects <- customobjects
			}
		}, error = function(err) stop(paste("slave failed to get mpcross objects: ", err)))

      str <- paste("Slave ",myID," of ",slaves)
      str <- paste(str, "class(objects)=", class(objects))
      str <- paste(str, " length(objects)=",length(objects))
      str <- paste(str, " lapply(objects,class)=",lapply(objects,class))

      tryCatch({

      markers <- nmrks
# ncol(objects[[1]]$founders)
      tileDimX <- ceiling(markers / gridDimX)
      tileDimY <- ceiling(markers / gridDimY)

      }, error = function(err) {
        stop(paste("slave setup failed: ", err))
      })

      calcTile <- function(x,y) {
      	       stopifnot(x>0, y>0, x<=gridDimX,y<=gridDimY)

      	       rect <- new.env()
      	       rect$x1 <- max(1, (x-1) * tileDimX + 1)
 	       rect$x2 <- min(x * tileDimX, markers)

      	       rect$y1 <- max(1, (y-1) * tileDimY + 1)
      	       rect$y2 <- min(y * tileDimY, markers)
      	       rect$sizex <- rect$x2 - rect$x1 + 1
      	       rect$sizey <- rect$y2 - rect$y1 + 1
      	       rect
      }
	getLock <- function(desc) 
	{
		dir.create(".locks", showWarnings=FALSE)
		lockname <- file.path(".locks",desc@description$filename)
		success <- FALSE
		attempts <- 0
		while(success != TRUE) 
		{
			success <- dir.create(lockname, showWarnings=FALSE)
			if (success != TRUE) 
			{
				Sys.sleep(1)
				attempts <- attempts + 1
			}
		}
	}
	releaseLock <- function(desc) 
	{
		lockname <- file.path(".locks",desc@description$filename)
		unlink(lockname, recursive = TRUE)
	}
	storeTile <- function(rect, rf) 
	{
		getLock(thetadesc)
		theta <- bigmemory::attach.big.matrix(thetadesc,backingpath=dir_base)
		theta[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$theta[1:rect$sizex,1:rect$sizey]
		bigmemory::flush(theta)
		releaseLock(thetadesc)

		getLock(loddesc)
		lod <- bigmemory::attach.big.matrix(loddesc,backingpath=dir_base)
		lod[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$lod[1:rect$sizex,1:rect$sizey]
		bigmemory::flush(lod)
		releaseLock(loddesc)

		getLock(lkhddesc)
		lkhd <- bigmemory::attach.big.matrix(lkhddesc,backingpath=dir_base)
		lkhd[rect$x1:rect$x2,rect$y1:rect$y2] <- rf$lkhd[1:rect$sizex,1:rect$sizey]
		bigmemory::flush(lkhd)
		releaseLock(lkhddesc)
	}
	sendResToMaster <- function(rect, res) 
	{
		Rmpi::mpi.send.Robj(rect,0,1) 		
		Rmpi::mpi.send.Robj(res,0,1)
        }

      inc <- function(x) { eval.parent(substitute(x <- x + 1)) }

      tryresult <- tryCatch({
      x <- 1
      while(x <= gridDimX) {
      	      y <- 1
     	      while(y <= gridDimY) {
	      	      if (y %% slaves == myID-1) {
		      	 str <- paste(str, "[" , x , "," , y , "] ")

			 tryCatch({
			 tileDim <- calcTile(x,y)
			 }, error = function(err) {
			    stop(paste("calcTile failed: ", err))
			 })


			 if(!exists("lineWeights")) lineWeights <- lapply(objects, function(x) vector(mode="integer", length=0))
			 if(!exists("r")) r <- c(0:20/200, 11:50/100)	  

			 tryCatch({
			 res <- mpestrfSubset(objects=objects,
						     gpu=gpu, r=r,
						     start1=tileDim$x1, finish1=tileDim$x2+1,
						     start2=tileDim$y1, finish2=tileDim$y2+1)
			 }, error = function(err) {
			    stop(paste("mpestrfSubset failed: ", err, class(objects)))
			 })

			 if (onlyMasterWrites) {
			   tryCatch({
			   sendResToMaster(tileDim,res)
			   }, error = function(err) {
			     stop(paste("sendResToMaster failed: ", err))
			   })
			 } else {
			   tryCatch({
			   storeTile(tileDim,res$rf)
			   }, error = function(err) {
			      stop(paste("storeTile failed: ", err))
			   })
			 }
		      }
	      	      inc(y)
	      }		  
	      inc(x)
       }      
       }, error = function(err) return(list(mesg=paste(str,err),result=NULL)))
       return(list(mesg=paste(str),result=res))
}


mpestrfSubset <-
function(objects, r, gpu, lineWeights, start1, finish1, start2, finish2)
{ 
  	if(missing(lineWeights))
	{
		lineWeights <- lapply(objects, function(x) rep(1, nrow(x$finals)))
	}
	lineWeights <- lapply(lineWeights, as.numeric)
	if (missing(objects))	stop("Missing a required argument for this function")
	if (missing(r)) r <- c(0:20/200, 11:50/100)
	for(i in 1:length(objects))
	{
		if (!inherits(objects[[i]], "mpcross")) stop("Object must be of class mpcross")
		if(is.null(colnames(objects[[i]]$finals)) && !is.null(colnames(objects[[i]]$founders)))
		{
			warning("Entry finals did not have column names, replacing with column names from entry founders")
			colnames(objects[[i]]$finals) <- colnames(objects[[i]]$founders)
		}
		if(is.null(colnames(objects[[i]]$founders)) && !is.null(colnames(objects[[i]]$finals)))
		{
			warning("Entry founders did not have column names, replacing with column names from entry finals")
			colnames(objects[[i]]$founders) <- colnames(objects[[i]]$finals)
		}
		if(is.null(colnames(objects[[i]]$founders)) && is.null(colnames(objects[[i]]$finals)))
		{
			stop("One of the values founders and finals must have column names")
		}
		if(ncol(objects[[i]]$founders) != ncol(objects[[i]]$founders))
		{
			stop("Founders and finals data matrices must have the same number of columns")
		}
		if(any(colnames(objects[[i]]$founders) != colnames(objects[[i]]$finals))) 
		{
			stop("Columns names for object$founders and object$finals were inconsistent")
		}
		objects[[i]]$pedigree <- convertped(objects[[i]]$pedigree)

		n.founders <- nrow(objects[[i]]$founders)
		n.loci <- ncol(objects[[i]]$founders)
		n.finals <- nrow(objects[[i]]$finals)

		if(missing(gpu)) gpu <- FALSE
		if(class(objects[[i]]$founders) != "matrix") objects[[i]]$founders <- as.matrix(objects[[i]]$founders, rownames.force=TRUE)
		if(class(objects[[i]]$finals) != "matrix") objects[[i]]$finals <- as.matrix(objects[[i]]$finals, rownames.force=TRUE)

		if(mode(objects[[i]]$founders) != "integer") mode(objects[[i]]$founders) <- "integer"
		if(mode(objects[[i]]$finals) != "integer") mode(objects[[i]]$finals) <- "integer"

		if(class(objects[[i]]$id) != "integer") objects[[i]]$id <- as.integer(objects[[i]]$id)
		if(class(objects[[i]]$fid) != "integer") objects[[i]]$fid <- as.integer(objects[[i]]$fid)
		if(missing(start1) || missing(finish1))
		{
			marker1Range <- c(1, n.loci+1)
		}
		else marker1Range <- c(start1, finish1)
		if(missing(start2) || missing(finish2))
		{
			marker2Range <- c(1, n.loci+1)
		}
		else marker2Range <- c(start2, finish2)
		
#		if(any(apply(objects[[i]]$founders, 2, function(x) length(unique(x)) == 1))) stop("Non-segregating markers must be removed prior to calculation of recombination fractions")
	}
	rpairs <- .Call("rfhaps", objects, r, marker1Range, marker2Range, lineWeights, gpu, -2, PACKAGE="mpMap")
	output <- objects[[1]]
	if(length(objects) > 1)
	{
		markers <- colnames(objects[[1]]$founders)
		n.markers <- length(markers)
		output$founders <- matrix(0L, nrow=0, ncol=n.markers)
		colnames(output$founders) <- markers
		output$finals <- matrix(0L, nrow=0, ncol=n.markers)
		colnames(output$finals) <- markers
		output$pedigree <- NULL
		output$pheno <- NULL
	}
	output$rf <- rpairs
	return(output)
}

