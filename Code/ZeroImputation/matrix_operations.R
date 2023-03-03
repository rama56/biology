
drop_values_probabilistic <- function(sparse_matrix, zero_prob){
  dense_matrix = as.matrix(sparse_matrix)
  dimensions = dim(dense_matrix)
  n_rows = dimensions[1]
  n_cols = dimensions[2]
  
  zero_one_matrix = matrix(sample(c(0,1), size=n_rows*n_cols,
                                  replace=TRUE, prob=c(zero_prob,1-zero_prob)),
                           nrow = n_rows )
  
  artificial_matrix = zero_one_matrix * dense_matrix
  artificial_sparse_matrix = as(artificial_matrix, "dgCMatrix")
  return(artificial_sparse_matrix)
}

find_distance <- function(matrix_1, matrix_2){
  diff = matrix_1 - matrix_2
  squared_diff = diff * diff
  total_squared_distance = sum(squared_diff)
  
  dimensions = dim(squared_diff)
  n_rows = dimensions[1]
  n_cols = dimensions[2]
  
  mean_squared_distance = total_squared_distance/(n_rows*n_cols)
  return(mean_squared_distance)
}