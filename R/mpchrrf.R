#' Return Chromosome Recombination Fractions
#' 
#' Returns the matrix of recombination fraction estimates between markers in a specified linkage group
#' @export
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @param chr A vector of the chromosome(s) of interest
#' @return A list with a component for each element of \code{chr}. Each component is the matrix of recombination fractions for the markers in that linkage group
#' @examples
#' map <- qtl::sim.map(len=rep(10,2), n.mar=30, eq.spacing=FALSE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1)
#' dat.rf <- mpestrf(sim.dat)
#' dat.lg<-mpgroup(dat.rf, 2)
#' rf1<-mpchrrf(dat.lg, chr=1)

mpchrrf<-function(object, chr){
  if (!inherits(object, "mpcross")) 
    stop("Object must be of class mpcross")
  if (missing(object)|missing(chr))
    stop("Missing a required argument for this function")
  if (is.null(object$rf)){
    stop("Recombination fractions need to be estimated first")
  }
  if (is.null(object$lg)){
    stop("Linkage groups must be estimated first")
  }
  if (is.character(chr)){
    if (is.null(object$map)){
      stop("Object does not have a map to match chromosome names to")
    }
    tmp.chr<-match(chr, names(object$map))
    if (any(is.na(tmp.chr))){
      stop(paste("Chromosome(s)",paste(chr[which(is.na(tmp.chr))], collapse=","), "not present on object"))
    }
    chr<-tmp.chr
  }
  check.lg<-match(chr, object$lg$groups)
  if (any(is.na(check.lg))) {
    bad.chr<-chr[which(is.na(check.lg))]
    stop(paste("Chromosome(s)",paste(bad.chr, collapse=","), "not present on object"))
  }
  lg<-lapply(chr, function(x){which(object$lg$groups==x)})
  rfs<-lapply(lg, function(x){object$rf$theta[x,x]})
  if(!is.null(object$map)){
    names(rfs)<-names(object$map)[chr]
  }
  return(rfs)
}
