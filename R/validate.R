validate <- function(mpcross, checkPedigree = TRUE, checkRF = FALSE, checkLG = FALSE, checkFID = TRUE)
{
	result <- .Call("validate", mpcross, checkPedigree, checkRF, checkLG, checkFID)
	if(!is.null(result))
	{
		stop(result)
	}
}