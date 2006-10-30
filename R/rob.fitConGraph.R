"rob.fitConGraph" <-
function(amat, w.par, smpl.frm, it.limit = 200, tol = 1e-06)
{
  cgm <- rob.estimConGraph(amat = amat, w.par = w.par, smpl.frm = smpl.frm,
                            it.limit = it.limit, tol = tol)

  rbt <- rob.browne.test(estcov = cgm$Sigma, amat = amat,
                         w.par = w.par, smpl.frm = smpl.frm)

  return(list(mhat = cgm$mu,
              Shat = cgm$Sigma,
              w.vec = cgm$weight,
              w.par = w.par,
              it = cgm$it,
              tstat = rbt$chi.val,
              df = rbt$df,
              p.val = rbt$p.val))
}

