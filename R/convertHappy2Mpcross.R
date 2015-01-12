## Available on http://mus.well.ox.ac.uk/magic/
#source('magic.R')
#source('happy.preCC.R')
#tmp <- load.condensed.database()

## Create an mpcross object from happy database
convertHappy2mpMap <- function(cd=tmp, phen) 
{
   nfounders = length(cd$strains)
   nfinals = length(cd$subjects)
   nmrk = length(cd$additive$matrices) 
   fid = 1:nfounders   
   nchr = length(cd$chrs)
   map <- list()
   for (i in 1:nchr) {
      map[[i]] <- cd$additive$map[which(cd$additive$chromosome==i)]
      names(map[[i]]) <- cd$additive$markers[which(cd$additive$chromosome==i)]
   }
   names(map) <- paste("Chr", 1:nchr, sep="")
   class(map) <- "map"

   ### NOTE: AT the moment assuming we just have biallelic markers - need to deal with multiple alleles differently
   genfou <- list()
   for (i in 1:nchr) {
      temp2 <- read.table(paste("chr", i, ".MAGIC.markers", sep=""), sep="\t")
      temp3 <- read.table(paste("chr", i, ".MAGIC.alleleinfo", sep=""), sep="\t")
      temp3 <- temp3[!is.na(temp3[,2]),]
      genfou[[i]] <- matrix(nrow=nfounders, ncol=nrow(temp2))
      for (j in 1:ncol(genfou[[i]])) {
         genfou[[i]][,j] <- temp3[(j*2-1:0),2][(temp3[j*2-1,3:ncol(temp3)]==0)+1]
      }
      rownames(genfou[[i]]) <- cd$strains
      colnames(genfou[[i]]) <- temp2[,2]
   }
   
   ## Note that these will be in terms of 1/2/3/4 - levels of factor ACGT

   ## 1. finals
   genfin <- list()
   for (i in 1:nchr) 
   {
      temp <- read.table(paste("chr", i, ".MAGIC.data", sep=""), row.names=1)
      temp <- temp[,6:ncol(temp)]
      temp <- temp[, seq(1, ncol(temp), 2)]
      genfin[[i]] <- matrix(nrow=nrow(temp), ncol=ncol(temp))
      rownames(genfin[[i]]) <- rownames(temp)
      colnames(genfin[[i]]) <- colnames(genfou[[i]])
      genfin[[i]][temp=="A"] <- 1
      genfin[[i]][temp=="C"] <- 2
      genfin[[i]][temp=="G"] <- 3
      genfin[[i]][temp=="T"] <- 4
      genfin[[i]] <- genfin[[i]][, match(names(map[[i]]), colnames(genfin[[i]]))]
   }
   fin <- do.call("cbind", genfin)
   fin <- fin[match(cd$subjects, rownames(fin)),]
   for (i in 1:nchr) genfou[[i]] <- genfou[[i]][, match(names(map[[i]]), colnames(genfou[[i]]))]
   fou <- do.call("cbind", genfou)

   ## Make up a pretend pedigree since it doesn't really matter
   ped <- cbind(1:19, rep(0, 19), rep(0, 19), rep(0,19))
   ped <- rbind(ped, cbind(19+1:nfinals, sample(1:19, nfinals, replace=T), sample(1:19, nfinals, replace=T), rep(1, nfinals)))
   colnames(ped) <- c("id", "Male", "Female", "Observed")
   rownames(ped) <- c(rownames(fou), rownames(fin))
   ped <- as.data.frame(ped)
   
   id <- which(ped[,4]==1)

   prob <- list()
   allprob <- do.call("cbind", lapply(cd$additive$matrices, function(x) x[[1]]))
   rownames(allprob) <- cd$subjects
   index <- 0
   for (i in 1:nchr) {
      prob[[i]] <- allprob[, index+1:(length(map[[i]])*nfounders)]
      index <- index+length(map[[i]])*nfounders
      colnames(prob[[i]]) <- paste(rep(names(map[[i]]), each=nfounders), ", Founder ", 1:nfounders, sep="")
   }
   attr(prob, "map") <- map 
   attr(prob, "step") <- 0
   attr(prob, "program") <- "happy"
   attr(prob, "mapfx") <- "haldane"
   attr(prob, "mrkpos") <- TRUE
   names(prob) <- names(map)
   mpc <- mpcross(founders=fou, finals=fin, fid=fid, id=id, pedigree=ped)
   mpc$map <- map
   mpc$prob <- prob
   class(mpc) <- c("mpprob", "mpcross")

   if(!missing(phen)) {
      mpc <- subset(mpc, lines=match(rownames(phen), rownames(fin)))
      mpc$pheno <- phen
   }
   nalleles <- apply(mpc$finals, 2, function(x) length(table(x)))
   attr(mpc$prob, "mrkpos") <- TRUE

   ## one marker is monomorphic - remove this. For some reason subset not working - do it manually
   index <- grep(colnames(mpc$finals)[which(nalleles==1)], names(mpc$map[[5]]))
   mpc$finals <- mpc$finals[,-which(nalleles==1)]
   mpc$founders <- mpc$founders[, -which(nalleles==1)]
   mpc$map[[5]] <- mpc$map[[5]][-index]
   mpc$prob[[5]] <- mpc$prob[[5]][, -(index-1)*19+1:19]
   attr(mpc$prob, "map") <- mpc$map
   mpc
}

