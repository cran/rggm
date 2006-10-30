"rob.browne.test" <-
function(estcov, amat, w.par, smpl.frm)
{

  complete.ug <- function(amat) {
    amat[] <- 1
    sat <- amat - diag(nrow(amat))
    vnames <- vertices(amat)
    dimnames(sat) <- list(vnames, vnames)

    return(sat)
  }

  sat <- complete.ug(amat)
  
  # coerce smpl.frm to be as data frame
  smpl.frm <- as.data.frame(smpl.frm)

  nost.inv   <-
    solve(rob.estimConGraph(amat = sat, w.par = w.par, smpl.frm = smpl.frm)$Sigma)
  struct.inv <- solve(estcov)

  diff <-vech(nost.inv - struct.inv)

  gamma <-
    rob.asymv.sigma(estcov = nost.inv, amat = sat, w.par = w.par)
  gamma.inv <- solve(gamma)

  dlt <- rob.delta(amat)

  V <- gamma.inv -
    gamma.inv %*% dlt %*%
      solve(t(dlt) %*% gamma.inv %*% dlt) %*%
        t(dlt) %*% gamma.inv

  chi.val <-  nrow(smpl.frm) * t(diff) %*% V %*% diff
  df      <- sum(diag(gamma %*% V))
  p.val   <- 1 - pchisq(chi.val, df)

  return(list(chi.val = chi.val,
              df = df,
              p.val = p.val))
}

