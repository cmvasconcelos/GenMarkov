
#MMC estimation with Nicolau (2014) specification through maxLik()
#' Title
#'
#' @param y matrix of categorical data sequences
#' @param initial numerical vector of initial values
#' @param nummethod Numerical maximisation method, currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR" (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. The default method is "BFGS". For more details see maxLik.
#'
#' @return The function returns the parameter estimates, standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
#' @export
#'
#' @examples
#' set.seed(1234)
#' s1 <- sample(c(1,2), 500, replace=TRUE)
#' s2 <- sample(c(1,2), 500, replace=TRUE)
#' multi.mtd_probit(y = cbind(s1,s2), initial=c(1,1,1), nummethod='bfgs')
#'
multi.mtd_probit <- function(y, initial, nummethod = "bfgs") {

  #Define number of equations to estimate
  nr.eq <- ncol(y)

  #Calculate frequency transition matrices and vectorize
  ni <- apply(ProbMatrixF(y), 3, function(x) matrixcalc::vec(x))

  #Calculate probability transition matrices and vectorize
  qi <- apply(ProbMatrixQ(y), 3, function(x) matrixcalc::vec(x))

  #Define indices to select appropriate columns of ni and qi, for each equation
  index_vector = 1:ncol(ni)

  split_result <- lapply(split(1:length(index_vector), rep(1:nr.eq, each = nr.eq)),
                         function(indices) index_vector[indices])


  #Obtain etas for each equation
  results <- lapply(split_result, FUN = function(i) {
    neq <- ni[,i]
    qeq <- qi[,i]

    #Define log-likelihood
    LogLikelihood_p <- function(eta, n = neq, qn = qeq){
      ll <- 0

      denom <- unlist(lapply(split_result, function(i) {
        stats::pnorm(sum(cbind(rep(1,nrow(qi[,i])), qi[,i]) %*% eta))
      }))

      eta_mat <- log(stats::pnorm(cbind(rep(1,nrow(qn)), qn) %*% eta)/sum(denom))

      ll <- sum(t(n)%*%eta_mat)

      return(ll)
    }

    #Maximize log-likelihood
    otim <- maxLik::maxLik(LogLikelihood_p, start = initial, method = nummethod)

    res <- summary(otim)

    return(res)
  })

  #Return results
  names(results) <- paste('Equation', 1:nr.eq)

  return(results)

}
