"duplication" <-
function(n)
{
  m <- n * (n+1) / 2
  a <- diag(n)
  a[lower.tri(a, TRUE)] <- 1:m
  a[upper.tri(a)] <- a[lower.tri(a)]
  
  a.vec <- as.vector(a)

  D <- matrix(rep(0, (n*n) * m), (n*n), m)
  for(r in 1:(n*n)) {
    D[r, a.vec[r]] <- 1
  }

  return(D)
}

