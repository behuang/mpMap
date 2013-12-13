\name{mpMap-internal}
\Rdversion{1.1}
\alias{check_qtl}
\alias{check_ped}
\alias{combine_chr}
\alias{combine_ibd}
\alias{combine_rf}
\alias{convertped}
\alias{cleanrf}
\alias{CR_gen_geno}
\alias{CR_calcLD}
\alias{generate_error}
\alias{generate_obs}
\alias{generate_pheno}
\alias{getpval}
\alias{wald.test.asreml}
\alias{cintern}
\alias{scan3pt}

\title{ Internal mpMap functions} 
\description{
Internal MPmap functions. These are generally not to be called by the user.}
\usage{
check_qtl(QTL, map, n.founders)
check_ped(pedigree, nFounders)
combine_chr(map, map.function)
combine_ibd(obsgeno)
convertped(pedigree)
CR_gen_geno(map, map.function, ped, seed, transpos, transval)
CR_calcLD(finals, founders, pedigree, pairs, rmat)
cleanrf(mpcross)
generate_error(geno, error.prob)
generate_obs(geno, map, full.prob, fg, transpos, transval, founderld)
generate_pheno(n.founders, qtlgeno, qtleffects, vare, n.ind)
getpval(model)
scan3pt <- function(mpcross, chr, step=0, mapfx=c("haldane", "kosambi"))
wald.test.asreml(object, cc)
cintern(cc, tau, vrb, sigma2)
}
\keyword{internal}
