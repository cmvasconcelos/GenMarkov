#MMC estimation dependent on covariate
#' Title
#'
#' @param y matrix of categorical data sequences
#' @param x matrix of covariates
#' @param initial numerical vector of initial values.
#' @param ... additional arguments to be passed down to auglag()
#'
#' @return The function returns the parameter estimates, standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
#' @export
#'
#' @examples
#' set.seed(1234)
#' s1 <- sample(c(1,2), 500, replace=TRUE)
#' s2 <- sample(c(1,2), 500, replace=TRUE)
#' x <- rnorm(500)
#' mmcx(y = cbind(s1,s2), x = cbind(x), initial=c(1,1))
#'
mmcx <- function(y, x, initial, ...) {

  r.eq <- ncol(y)
  m1 <- max(y)

  q <- ProbValuesXDependent(s = y, x = x) #Calculate P(St | St-1, x) for all equations

  #Define indices to select appropriate values, for each equation
  index_vector = 1:(ncol(y) * ncol(y))

  split_result <- lapply(split(1:length(index_vector), rep(1:ncol(y), each = ncol(y))),
                         function(indices) index_vector[indices])

  j = 1
  for (i in split_result) {

    LogLikelihood_x <- function(lambda, qi = q[, i]) {

      qi_lambda <- qi %*% lambda

      qi_lambda[qi_lambda < 0] <- 1

      ll <-  sum(-log(qi_lambda))

      ll
    }

    #Impose restrictions
    he <- function(lambda) {
      h <- rep(NA, 1)
      h[1] <- sum(lambda) - 1
      h
    }

    hi <- function(lambda) {
      w <- rep(NA, 1)
      for (i in 1:length(lambda)) {
        w[i] <- lambda[i]
      }
      w
    }

    #Jacobians of restrictions to improve optimization in auglag()
    he.jac <- function(lambda) {
      j <- matrix(NA, 1, length(lambda))
      j[1,] <- rep(1, length(lambda))
      j
    }

    hi.jac <- function(lambda) {
      j <- diag(1, length(lambda), length(lambda))
      j
    }

    #Optimization through auglag() function
    opt <-
      alabama::auglag(
        par = initial,
        fn = LogLikelihood_x,
        hin = hi,
        heq = he,
        heq.jac = he.jac,
        hin.jac = hi.jac,
        control.outer = list(trace = FALSE),
        ...
      )

    #Only performs inference if the model reaches convergence
    if (any(opt$ineq < 0) | any(round(opt$par, 1) == 1.0)) {
      ll <- '.'
      hessian <- '.'
      lambdahat <- '.'
      inf <- list(warning = 2)
    } else{
      hessian <- -opt$hessian
      ll <- -opt$value
      lambdahat <- opt$par

      inf <- Inference_x(hess=hessian, lambda= lambdahat)
    }

    output.table(lambdahat, inf$se, inf$zstat, inf$pvalue, ll, j, flag=inf$warning)

    j = j + 1
  }


}
