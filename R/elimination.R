"elimination" <-
function(n)
{
  d <- duplication(n)
  return(solve(t(d) %*% d) %*% t(d))
}

