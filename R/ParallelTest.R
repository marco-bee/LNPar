#' Profile-based testing for a Pareto tail
#'
#' This function computes the bootstrap test for the null hypothesis of
#' a pure lognormal distribution versus the alternative of a lognormal-Pareto
#' mixture, where the parameters of the latter are estimated via maximum profile
#' likelihood. Implemented via parallel computing.
#' @param nboot number of bootstrap replications.
#' @param y observed data.
#' @param obsTest value of the test statistics computed with the data under analysis.
#' @param minRank minimum possible rank of the threshold.
#' @return A list with the following elements:
#'
#' LR: nboot simulated values of the llr test under the null hypothesis.
#'
#' pval: p-value of the test.
#' @export
#' @examples
#' \donttest{minRank = 90
#' mixFit <- LPfitProf(TN2016,minRank,0)
#' ell1 <- mixFit$loglik
#' estNull <- c(mean(log(TN2016)),sd(log(TN2016)))
#' ellNull <- sum(log(dlnorm(TN2016,estNull[1],estNull[2])))
#' obsTest <- 2*(ell1-ellNull)
#' nboot = 2
#' TestRes = ParallelTest(nboot,TN2016,obsTest,minRank)}

ParallelTest = function(nboot,y,obsTest,minRank)
{
nreps.list <- sapply(1:nboot, list)
n = length(y)
chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
if (nzchar(chk) && chk == "TRUE") {
  n.cores <- 2L
} else {
  n.cores <- parallel::detectCores()
}
clust <- parallel::makeCluster(n.cores)
LRVec = rep(0,nboot)
temp <- parallel::parLapply(clust,nreps.list,LPtest,n,mean(log(y)),sd(log(y)),minRank)
parallel::stopCluster(cl=clust)
for (i in 1:nboot)
{
  LRVec[i] = as.vector(unlist(temp[[i]]))
}
obsp <- length(LRVec[LRVec>obsTest])/nboot
results <- list(LR=LRVec,pVal=obsp)
return(results)
}
