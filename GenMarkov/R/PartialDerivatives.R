#Calculate first partial derivatives for numerical maximization - (This function is adapted from march::march.mtd.construct())
#' Title
#'
#' @param ni numerical matrix of number of transitions between states and data sequences
#' @param qi numerical matrix of transitions probabilities between states and data sequences
#' @param lambda numerical vector
#'
#' @return numerical vector with partial derivates values of each lambda
#' @export
#'
#' @examples
#' #######
PartialDerivatives <- function(ni, qi, lambda) {

  den <- apply(qi, 1, FUN = function(q){ q %*% lambda })

  pd_lambda <- colSums(ni*(qi/den))
  pd_lambda[is.infinite(pd_lambda)] <- 0

  return(pd_lambda)

}
