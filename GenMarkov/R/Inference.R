#Second partial derivatives computation to perform inference on the parameters
#' Title
#'
#' @param ni numerical matrix of transitions frequencies
#' @param qi numerical matrix of transitions probabilities
#' @param lambda numerical vector
#'
#' @return list with standard error, z-statistics and pvalues associated with each estimate
#' @export
#'
#' @examples
#' #########
Inference <- function(ni, qi, lambda) {

  den <- apply(qi, 1, FUN = function(q){ (q %*% lambda)^2 })

  num <- apply(qi, 2, FUN = function(qi_c){(qi_c*qi)}, simplify = FALSE)

  hess_list = apply(array(unlist(num), dim=c(nrow(ni),length(lambda),length(lambda))), 3,
               FUN = function(num){
                 colSums(ni*(-(num/den))) }, simplify = FALSE)

  hess = matrix(unlist(hess_list), nrow = length(lambda), ncol=length(lambda), byrow = TRUE)

  hessinv <- solve(-hess)

  var <- diag(hessinv)
  se <- sqrt(var)
  zstat <- lambda / se
  pvalue <- 2 * (1 - stats::pnorm(abs(zstat)))


  return(l = list(
    se = se,
    zstat = zstat,
    pvalue = pvalue
  ))

}
