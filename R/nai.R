#' Count how many generations of advanced intercross (AI) are in a pedigree. The number of AI generations can differ for each individual in the pedigree. 
#'
#' Counts the number of generations of breeding preceding selfing and 
#' subtracts off the number necessary to minimally mix the founders' genomes. The number of founders for the design is extracted from the pedigree, and must be a power of 2. 
#' @export
#' @param pedigree Pedigree for a multi-parent cross. Can be generated using \code{\link[mpMap]{sim.mpped}}
#' @param IDs A numeric vector containing the rows within the pedigree, for which we are interested in the number of AI generations. Defaults to all. 
#' @param nFounders Number of founders
#' @return A vector of integers, containing the number of generations of advanced intercrossing after mixing stage but before selfing, for each line in the final population
#' @seealso \code{\link[mpMap]{sim.mpped}}
#' @examples
#' sim.map <- list(Chr1=seq(0,100,10))
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' nai(sim.ped)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1, 5)
#' nai(sim.ped)

nai <- function(pedigree, IDs, nFounders)
{
	if(inherits(pedigree, "mpcross")) pedigree <- pedigree$pedigree
	if(missing(IDs)) IDs <- which(pedigree[,4] == 1)
	#It's easy to compute the total number of generations before selfing, and easy to compute the number of generations (3 or 2) at the very top of the design before the Advanced Intercrossing
	#step. So just subtract one from the other.
	if(missing(nFounders))	nFounders <- sum(pedigree[,2] == 0 & pedigree[,3] == 0)
	#If nFounders is not a power of 2, report an error
	if(as.integer(log(nFounders, base=2)) != log(nFounders, base=2))
	{
		stop("Number of founders must be a power of 2")
	}
	nTopGenerations <- log(nFounders, base=2)
	f <- function(id)
	{
		pedigreeRow <- id
		while(pedigree[pedigreeRow, 2] == pedigree[pedigreeRow, 3])
		{
			pedigreeRow <- pedigree[pedigreeRow, 2]
		}
		#detect NA values in the pedigree
		if((pedigree[pedigreeRow, 2] != pedigree[pedigreeRow, 2]) || (pedigree[pedigreeRow, 3] != pedigree[pedigreeRow, 3]))
		{
			stop("Invalid pedigree")
		}
		ngen <- 0
		while(pedigree[pedigreeRow, 2] > 0)
		{
			pedigreeRow <- pedigree[pedigreeRow, 2]
			ngen <- ngen + 1
		}
		return(ngen - nTopGenerations);
	}
	return(sapply(IDs, f))
}
