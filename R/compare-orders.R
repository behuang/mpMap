#' Compare potential orders for linkage groups 
#' 
#' Compares potential orderings on the number of expected crossovers or 
#' likelihood values. Rewrite of older version.  
#' @export
#' @importFrom stats window
#' @param cross Object of class \code{mpcross} or \code{cross}
#' @param chr Selected chromosomes
#' @param orders Orders to be compared
#' @param method Method for comparison. Note that "likelihood" is much more time-consuming than "countXO"
#' @details Uses functions from R/qtl in order to compare possible orderings on the basis of the number of obligate crossovers or log-likelihood values. 
#' @return The matrix of orders with crossover counts appended. 
#' @references R/qtl
#' @seealso \code{\link[mpMap]{mporder}}, \code{\link[qtl]{countXO}}
#' @examples
#' map <- qtl::sim.map(len=100, n.mar=11, eq.spacing=TRUE, include.x=FALSE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1)
#' compare_orders(sim.dat, chr=1, orders=rbind(1:11, c(1:3, 6:4, 7:11)))

compare_orders <- function(cross, chr, orders, method=c("countXO", "likelihood"))
{
    if (missing(method)) method <- "countXO"
    if (inherits(cross, "mpcross")) {
	write2cross(cross, "tmp", chr=chr)
	mpcross <- cross
	cross <- qtl::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(cross, "type"))
    }

    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    if (missing(chr)) {
        chr <- names(cross$geno)[1]
        warning("chr argument not provided; assuming you want chr ", 
            chr)
    }
    else {
        if (length(chr) > 1) 
            stop("ripple only works for one chromosome at a time.")
	if (length(intersect(chr, names(cross$geno)))==0)
            stop("Chr ", chr, "not found.")
    }
    cross <- subset(cross, chr = chr)
    chr.name <- names(cross$geno)[1]

    if (qtl::nmar(cross)[1] < 3) {
        warning("Less than three markers.")
        return(NULL)
    }
    n.mar <- totmar(cross)

 ### Essentially rewriting to subset down to the part of the cross
 ### where we're comparing orders, and then running countXO
 ### is it possible to reorder the cross in R/qtl or do we need 
 ### to do it in mpMap and retransmit each time? 

    # BEH
    if (missing(orders)) 
    	stop("Must input orders to compare\n")

    n.orders <- nrow(orders)
    if (n.orders > 49) 
        print.by <- 10
    else if (n.orders > 14) 
        print.by <- 5
    else print.by <- 2

    if (method == "likelihood") {
        loglik <- 1:n.orders
        chrlen <- 1:n.orders
        m <- seq(0, by = 5, length = n.mar)
        temcross <- cross
        if (is.matrix(cross$geno[[1]]$map)) 
            temcross$geno[[1]]$map <- rbind(m, m)
        else temcross$geno[[1]]$map <- m
        for (i in 1:n.orders) {
            if (i == 1) 
                cat("  ", n.orders, "total orders\n")
            if ((i%/%print.by) * print.by == i) 
                cat("    --Order", i, "\n")
            temcross$geno[[1]]$data <- cross$geno[[1]]$data[, 
                orders[i, ]]
            newmap <- qtl::est.map(temcross, chr, error.prob=1e-4, map.function="haldane", 
                m = 0, p = 0, maxit=4000, tol=1e-6, sex.sp=TRUE, FALSE)
            loglik[i] <- attr(newmap[[1]], "loglik")
            chrlen[i] <- diff(range(newmap[[1]]))
        }
        loglik <- (loglik - loglik[1])/log(10)
        o <- order(loglik[-1], decreasing = TRUE) + 1
        orders <- cbind(orders, LOD = loglik, chrlen)[c(1, o), 
            ]
    }
    else {
	oblxo <- vector(length=n.orders)
	for (i in 1:n.orders) {
		mpcross <- subset(mpcross, chr=chr.name)
	  if(i>1) {
		mpc <- subset(mpcross, markers=orders[i,])
		mpc$map <- list()
		mpc$lg <- NULL
		mpc$map[[chr.name]] <- rep(0, ncol(mpc$finals))
		names(mpc$map[[chr.name]]) <- colnames(mpc$finals)
		write2cross(mpc, "tmp")
        	cr <- qtl::readMWril("", "tmp.ril.csv", "tmp.founder.csv", type=attr(mpcross, "type"))	
	  	  } else cr <- cross
	  oblxo[i] <- sum(qtl::countXO(cr, chr.name, bychr=TRUE))
	}
        o <- order(oblxo[-1]) + 1

        orders <- cbind(orders, obligXO = oblxo)[c(1, o), ]
    }

    rownames(orders) <- c("Initial", paste(1:(nrow(orders) - 
        1)))
    class(orders) <- c("matrix")
    attr(orders, "chr") <- chr.name
    attr(orders, "window") <- window
    attr(orders, "method") <- method
    orders[, 1:n.mar] <- t(apply(orders[, 1:n.mar, drop = FALSE], 
        1, function(a) {
            n <- length(a)
            if ((1:n)[a == 1] > (1:n)[a == n]) 
                return(rev(a))
            else return(a)
        }))
    orders
}
