getAllFunnels <- function(mpcross)
{
	if(class(mpcross$founders) != "matrix") mpcross$founders <- as.matrix(mpcross$founders)
	if(class(mpcross$finals) != "matrix") mpcross$finals <- as.matrix(mpcross$finals)
	
	if(mode(mpcross$fid) == "numeric") mode(mpcross$fid) <- "integer"
	if(mode(mpcross$id) == "numeric") mode(mpcross$id) <- "integer"
	if(mode(mpcross$founders) == "numeric") mode(mpcross$founders) <- "integer"
	if(mode(mpcross$finals) == "numeric") mode(mpcross$finals) <- "integer"
	if(!inherits(mpcross, "mpcross")) stop("Input must have class mpcross")
	funnels <- .Call("getAllFunnels", mpcross, PACKAGE="mpMap")
	return(funnels)
}