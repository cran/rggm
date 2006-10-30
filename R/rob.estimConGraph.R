"rob.estimConGraph" <-
function(amat, w.par, smpl.frm, it.limit = 200, tol = 1e-06)
{
  ################################################################
  mu <- function(smpl.frm, w.par, w.vec) {
    dnm <- t(exp(w.par * w.vec)) %*% as.matrix(smpl.frm)
    nmr <- sum(exp(w.par * w.vec))
    if (is.infinite(nmr)) {
      print(paste("Denominator: ", dnm))
      print(paste("Numerator: ", nmr))
      stop("Convergence failed.\n")
    }

    return(dnm / nmr)
  }
  
  ################################################################
  sigma <- function(smpl.frm, mu, w.par, w.vec) {
    swped <- as.matrix(sweep(smpl.frm, 2, mu))
    exp.vec <- as.vector(exp(w.par * w.vec))

    p <- ncol(smpl.frm)
    n <- nrow(smpl.frm)
    sigma <- matrix(0, p, p)
    for (i in 1:n) {
      sigma <- sigma + exp.vec[i] * swped[i, ] %o% swped[i, ]
    }

    dnm <- mean(exp.vec) - (w.par / (w.par + 1)^(1 + 0.5 * p))
    sigma <- sigma / (n * dnm)

    return(sigma)
  }
  
  ################################################################
  weight <- function(smpl.frm, mu, Sigma) {
    S.inv <- solve(Sigma)
    swped <- sweep(smpl.frm, 2, mu)
    
    w.vec <- apply(swped, 1, function(swprow) {
      -0.5 * swprow %*% S.inv %*% swprow
    })
    
    return(w.vec)
  }
  
  ################################################################
  # Main routine
  p <- ncol(smpl.frm)
  n <- nrow(smpl.frm)

  tf <- TRUE

  # Initialze
  rM <- mean(smpl.frm)
  rS <- cov(smpl.frm)
  w.vec <- rep(1, n) # Conventional ML setting

  it <- 0
  repeat {
    it <- it + 1

    rS.old <- rS
    rM <-
      mu(smpl.frm = smpl.frm, w.par = w.par, w.vec = w.vec)
    rS <-
      sigma(smpl.frm = smpl.frm, mu = rM, w.par = w.par, w.vec = w.vec)
    rS <- fitConGraph(amat = amat, S = rS, n = n, pri = FALSE, alg = 2, tol = tol)$Shat

    w.vec <- weight(smpl.frm = smpl.frm, mu = rM, Sigma = rS)

    if (sum(abs(rS.old - rS)) < tol) {
      rslt <-
        list(mu = rM, Sigma = rS, weight = w.vec, w.par = w.par, it = it, success = tf)
      break
    }

    if (it > it.limit) {
      tf <- FALSE
      rslt <-
        list(mu = NULL, Sigma = NULL, weight = NULL, w.par = w.par, it = it, success = tf)
      break
    }
  }

  return(rslt)
}

