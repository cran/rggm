"rob.delta" <-
function(amat)
{
  # Return a matrix some columns are partialed out based on a given adajcent matrix
  # amat: an adjacent matrix
  # Value: a matrix of indices partilaed out

  p <- nrow(amat)

  amat <- diag(p) + amat
  del <- as.logical(amat[lower.tri(amat, TRUE)]) # 0: non-adjacent; 1: adjacent

  delta <- diag(0.5 * p * (p+1))
  return(delta[, del])
}

