#' Find Largest Complete mpcross object
#' 
#' Finds the largest subset of markers for which all recombination fraction estimates are available for each linkage group. The subsetted object can be returned, or a summary of marker reduction can be given.
#' @export
#' @useDynLib mpMap
#' @param object Object of class \code{mpcross}
#' @param output.type Type of output 
#' @return If the output type is "summary",  then a 2 x n.chr matrix is returned. The first row gives the number of markers present on that linkage group. The second row gives the number of markers in the maximal subset. If the output type is "object" then an \code{mpcross} object is returned. The object is a subset of the orginial object, with the largest set of markers where all pairs of recombination fractions estimates are availbe in each linkage group. 
#' @examples
#' map <- sim.map(len=1, n.mar=30, eq.spacing=FALSE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1)
#' dat.rf <- mpestrf(sim.dat)
#' dat.lg<-mpgroup(dat.rf, groups=1)
#' clean.summary <- mpclean(dat.lg)
mpclean<-function(object, output.type=c("summary", "object")){
  if (!inherits(object, "mpcross")) 
    stop("Object must be of class mpcross")
  if (missing(object))
    stop("Missing a required argument for this function")
  if (is.null(object$rf)){
    stop("Recombination fractions must be estimated first")
  }
  if (is.null(object$lg)){
    stop("Linkage groups must be estimated first")
  }
  if (missing(output.type)){output.type<-"summary"}
  n.chr<-length(object$lg$all.groups)
  output<-matrix(0,nrow=2, ncol=n.chr)
  lg<-object$lg$groups
  keep.markers<-c()
  for (i in 1:n.chr){
    lg.group<-which(lg==i)
    rf.group<-object$rf$theta[lg.group, lg.group]
    rf.group<-remove.NA(rf.group)
    output[1,i]<-length(lg.group)
    output[2,i]<-nrow(rf.group)
    keep.markers<-c(keep.markers, colnames(rf.group))
  }
  if (output.type=="summary"){
    rownames(output)<-c("Markers","Complete Markers")
    if (!is.null(object$map)){colnames(output)<-names(object$map)}
    return(output)
  }
  if (output.type=="object"){
    kind<-which(colnames(object$finals) %in% keep.markers)
    output<-subset(object, markers=keep.markers)
    output$lg$groups<-object$lg$groups[kind]
    output$lg$all.groups<-unique(object$lg$groups[kind])
    return(output)
  }  
}

remove.NA<-function(rf)
{
  rowcount<-apply(rf,1,function(x){sum(is.na(x))})
  while(any(is.na(rf))) { 
    delrow<-which.max(rowcount)
    decrement<-which(is.na(rf[delrow,]))
    rf<-rf[-delrow, -delrow]
    rowcount[decrement]<-rowcount[decrement]-1
    rowcount<-rowcount[-delrow]
  }
  return(rf)
}
