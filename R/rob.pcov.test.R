"rob.pcov.test" <-
function(estcov, amat, w.par, smpl.frm)
{
  ################################################################
  nC2 <- function(vec) {
    c1 <- function(vec) {
      rslt <- rep(vec[1], length(vec))
      if (length(vec) > 1)
        rslt <- append(rslt, c1(vec[-1]))
      return(rslt)
    }
    
    c2 <- function(vec) {
      rslt <- vec
      if (length(vec) > 1)
        rslt <- append(rslt, c2(vec[-1]))
      
      return(rslt)
    }

    list.1 <- c1(vec)
    list.2 <- c2(vec)

    return(lapply((1:length(list.1)),
                  function(i) c(list.1[i], list.2[i])))
  }

  pvalue <- function(value, avalue, sample.size) {
    tstat <- value/sqrt(avalue/n)
    pval <- 2*(1 - pnorm(abs(tstat)))
    
    return(list(tstat = tstat,
                pvalue = pval))
  }

  ################################################################
  # Main routine
  p <- nrow(amat)

  sigma.inv.trad <-
    solve(rob.estimConGraph(amat = amat, w.par = 0.0, smpl.frm = smpl.frm)$Sigma)
  pcor.trad <- diag(2, p) - cov2cor(sigma.inv.trad)

  sigma.inv <- solve(estcov)
  pcor <- diag(2, p) - cov2cor(sigma.inv)

  asym.sgm <- rob.asymv.sigma(sigma.inv, amat = amat, w.par = w.par)

  n <- nrow(smpl.frm)
  vals <- as.vector(t(rob.delta(amat)) %*% vech(sigma.inv))
  avals <- diag(asym.sgm)
  vs <- nC2(vertices(amat))[as.logical(vech(amat + diag(p)))]
  trslts <- lapply(1:length(vals), function(i) {
    pval <- pvalue(vals[i], avals[i], n)
    lst <- list(vs = vs[[i]],
                val = vals[i],
                aval = avals[i],
                tstat = pval$tstat,
                pval = pval$pvalue)

    return(lst)
  })
  
  return(trslts)
}

