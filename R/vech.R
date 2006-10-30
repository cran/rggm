"vech" <-
function(m, diag = TRUE)
{
  # Return the half vectororized elements of a given matrix.
  # m: a matrix
  
  return(as.matrix(m[lower.tri(m, diag)]))
}

