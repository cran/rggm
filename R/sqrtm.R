"sqrtm" <-
function(m)
{
  # Get a square root matrix of a square matrix
  # m: a square matrix

  x1 <- eigen(m)
  if(any(x1$values < 0)) {
    stop("A matrix should be positive definite.")
  }
  x1e <- sqrt(diag(x1$values))
  u1 <- cbind(x1$vectors)
  u1i <- solve(u1)
  return(u1 %*% x1e %*% u1i)
}

