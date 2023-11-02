#Calculate Frequency Transition Matrices (Ching (2002))

#' Title
#'
#' @param s numerical matrix with categorical data sequences
#'
#' @return numerical array with frequency transition matrices
#' @export
#'
#' @examples
#' ##############
ProbMatrixF <- function(s){
  n = nrow(s)
  m1 = max(s)
  m0 = min(s)
  num_cols = ncol(s)
  number_comb = num_cols*num_cols

  s_prev <- s[-nrow(s), ]
  s_next <- s[-1, ]

  column_combinations <- expand.grid(1:num_cols, 1:num_cols)

  all_list <- lapply(1:nrow(column_combinations), function(i) {
    selected_columns <- column_combinations[i, 2:1]
    combined_data <- data.frame(s_next[, selected_columns[,1]], s_prev[, selected_columns[,2]])
    return(combined_data)
  })

  f <- lapply(
    all_list,
    FUN = function(z){
      table(z[,1], z[,2])
    }
  )
  return(array(unlist(f), dim=c(m1, m1, number_comb)))
}
