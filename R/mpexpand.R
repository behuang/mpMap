#' Expand a binned object to one containing full marker information
#' 
#' Given a binned mpcross object, will expand the map to include all markers within bins at the same positions. If an additional mpcross object is input (containing genotype information for all markers), will also reformat that object according to the binned order. 
#' @export
#' @useDynLib mpMap
#' @param object Binned object of class \code{binmpcross} 
#' @param fullmpcross Object of class \code{mpcross} containing 
#' @return If \code{fullmpcross} is input, this function will return an object of class \code{mpcross} where bins have been expanded into the full set of markers included in each bin. Otherwise this function will return a map object where bins have been replaced by the full set of markers within a bin, all co-located. 
#' @seealso \code{\link{mpcollapse}}
mpexpand <- function(object, fullmpcross){
   if (!inherits(object, "binmpcross")) 
      stop("Must input binned object")
   if (is.null(object$map)) 
      stop("Must have created a map based on binned object, otherwise there is nothing to expand")
  if (!missing(fullmpcross))
   if (length(setdiff(colnames(fullmpcross$finals), object$bins$MarkerName))>0)  cat("Missing some markers in fullmpcross\n")

   ## make sure everything is in the correct order
   object <- subset(object, markers=1:ncol(object$finals)) 

   ## Expand out the map first; essentially putting all markers within bins at the same location as the bin marker

   oldmap <- object$map
   newmap <- list()
   binmarkers <- colnames(object$finals)
   for (i in 1:length(oldmap)) {
      binschr <- object$bins[object$bins$group==i,]
      if (!missing(fullmpcross)) 
         binschr <- binschr[binschr$MarkerName%in%colnames(fullmpcross$finals),]
      bm <- binschr$binMarkerName ## replicated binmarker names
      om <- oldmap[[i]][names(oldmap[[i]])%in%bm]
      newmap[[i]] <- rep(om, table(bm)[match(names(om), names(table(bm)))])
      names(newmap[[i]]) <- binschr$MarkerName
   }
   names(newmap) <- names(oldmap)
   class(newmap) <- "map"

   ## If no other object just return the map
   if (missing(fullmpcross)) return(newmap)
   ## Otherwise expand out the genotypes as well 
   output <- fullmpcross
   if (!is.null(object$rf)) cat("RF must be re-estimated once markers are unbinned\n")
   
   output$map <- newmap
   output <- mpgroup(output)
   class(output) <- "mpcross"
   output <- subset(output, markers=unlist(lapply(newmap, names)))
   return(output)
}
