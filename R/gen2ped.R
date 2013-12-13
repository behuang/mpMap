gen2ped <-
function(nfunnels=1, nperfam=50, nssdgen=6, nseeds=1, iripgen=0)
{
  if (nfunnels != 1) 
	stop("Only able to generate a single funnel")

  obs <- vector()
  # start with founders
  ped <- rbind(c(1,0,0), c(2,0,0))

  n1 <- nrow(ped)+1  
  ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nperfam)), rep(1, nperfam), rep(2, nperfam)))
  n2 <- nrow(ped)  
  
  # at this point have done all the mixing, need to do AI and then selfing
  if (iripgen>0)
  {
	  for (i in 1:iripgen)
	  {
		for (j in n1:n2)
		{
			ped <- rbind(ped, c(nrow(ped)+1, j, sample(setdiff(n1:n2, j), 1)))
		}
		n1 <- n2+1
		n2 <- nrow(ped)
	  }
  }
  
  #now the selfing
  obs <- rep(0, nrow(ped))
  for (i in 1:nperfam)
  { 
	#pick out relevant line at end of AI step
    index <- i+n1-1

    for (j in 1:nseeds)
    {
		obs <- c(obs, rep(0, nssdgen-1), 1)
    	ped <- rbind(ped, c(nrow(ped)+1, index, index))
		if (nssdgen>1)
		{
			ped <- rbind(ped, cbind(c((nrow(ped)+1):(nrow(ped)+nssdgen-1)), c(nrow(ped):(nrow(ped)+nssdgen-2)), c(nrow(ped):(nrow(ped)+nssdgen-2))))
		}
    }
  }

  # fourth column is whether individual was genotyped
  ped <- cbind(ped, obs)
  ped <- as.data.frame(ped)
  names(ped) <- c("id", "Male", "Female", "obs")

  return(ped)
}

