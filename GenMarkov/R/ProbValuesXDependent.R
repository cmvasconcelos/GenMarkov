#Calculate conditional probabilities (Values of non-homogeneous Markov Chains)
#' Title
#'
#' @param s numerical matrix with categorical data sequences
#' @param x Matrix of covariates (exogeneous variables)
#'
#' @return non-homogeneous probability transition matrix
#' @export
#'
#' @examples
#' #############
ProbValuesXDependent <- function(s, x) {
  m1 <- max(s)
  n <- nrow(s)

  #Create matrix with dummies for each state
  dummies_list <-
    apply(s, 2, function(x)
      fastDummies::dummy_cols(x, remove_selected_columns = TRUE))
  dummies <- matrix(unlist(dummies_list),
                    ncol = m1 * m1,
                    nrow = nrow(s))

  # Create all possible combinations of column indices
  combinations <- expand.grid(1:ncol(s), 1:ncol(dummies))
  #Order by the first variable
  combinations <- combinations[order(combinations$Var1),]

  # Extract columns from S and S_L based on the combinations
  combined_list <- lapply(1:nrow(combinations),  function(i, x) {
    cbind(s[, combinations[i, 1]], x, dummies[, combinations[i, 2]])
  }, x = x)

  #Estimate conditional probabilities of the type: P(Sjt | Sjt-1 = k), being k = 1, ..., m1
  #If m1==2, then we only have two states, we can resort to a logistic regression and define the dependent variable
  #as 1 if S=1 and 0 if S=2.

  if (m1 == 2) {
    estimate_condprobs <- sapply(combined_list, function(data) {
      #Define dependent variable
      y = ifelse(data[, 1] == 1, 1, 0)

      #Define lagged St
      s_l = Hmisc::Lag(data[, 3])

      #Estimate logistic regression
      res = stats::glm(y[s_l == 1] ~ data[, 'x'][s_l == 1], family = stats::binomial(link = 'logit'))

      #Extract fitted values
      px1 = cbind(stats::fitted(res), 1 - stats::fitted(res))
      #Create matrix to store fitted values
      px = cbind(rep(0, length(s_l[-1])), rep(0, length(s_l[-1])))
      #Input fitted values where St-1 = s
      px[s_l[-1] == 1, ] = px1

      return(as.matrix(px))
    }, simplify = "array")

  } else{
    #If m1>2, we have more than two states, we will resort to a multinomial logistic regression
    estimate_condprobs <- sapply(combined_list, function(data) {
      #Define dependent variable
      y <- factor(data[, 1], levels = 1:max(data[, 1]))

      #Define lagged St
      s_l = Hmisc::Lag(data[, 3])

      #Estimate multinomial logistic regression
      res = nnet::multinom(y[s_l == 1] ~ data[, 'x'][s_l == 1, ])

      #Extract fitted values
      px1 = matrix(
        stats::fitted(res),
        ncol = ncol(stats::fitted(res)),
        nrow = nrow(stats::fitted(res))
      )
      #Create matrix to store fitted values
      px = matrix(rep(0, m1 * length(s_l[-1])),
                  nrow = length(s_l[-1]),
                  ncol = m1)
      #Input fitted values where St-1 = s
      px[s_l[-1] == 1, ] = px1

      return(as.matrix(px))
    })
  }

  # Create combinations of components index to obtain the probabilities P(Sjt | Sjt-1)
  combined_sum_probs <-
    matrix(
      data = 1:dim(estimate_condprobs)[3],
      ncol = max(s),
      byrow = TRUE
    )

  #Create list of all probabilities P(Sjt | Sjt-1)
  probs_all <- lapply(1:nrow(combined_sum_probs), function(index) {
    select = combined_sum_probs[index, ]
    apply(estimate_condprobs[, , select], c(1, 2), sum)
  })


  #Create combinations of components index to extract the realized probabilities P(Sjt | Sjt-1)
  combine_probs <-
    matrix(data = c(1:length(probs_all), rep(1:m1, each = ncol(s))), ncol =
             2)

  #Extract realized values of probabilities, according to Sjt-1 and Sjt
  probs <- sapply(1:nrow(combine_probs), function(i) {
    index = combine_probs[i, ]
    t_values = 2:(n - 1)
    col_values = s[t_values, index[[2]]]

    p <- sapply(t_values, function(t) {
      probs_all[[index[1]]][t - 1, col_values[t - 1]]
    })
    p
  })

  return(probs)

}
