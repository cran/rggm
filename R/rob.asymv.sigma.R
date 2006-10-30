"rob.asymv.sigma" <-
function(estcov, amat, w.par) {

  ################################################################
  bindex <- function(vec)
  {
    list.1 <-
      as.vector(sapply(vec,
                       function(item) rep(item, length(vec)) ))
    list.2 <- rep(vec, length(vec))
    return(lapply((1:length(list.1)),
           function(i) c(list.1[i], list.2[i])))
  }
  
  ################################################################
  asymv.v <- function(w.par, p) {

    bicond <- function( bi ) {
      # return the list-value matrix indexed by TRUE/FALSE
      tf <- 
        sapply( bi, function( item.1 ) {
          lapply( bi, function( item.2 ) {
            tf.vec <- rep( FALSE, 3 )
            item.1 <- sort( item.1 ); item.2 <- sort( item.2 )
            if( item.1[1] == item.1[2] ) {
              if( (item.1[1] == item.2[1]) && (item.1[2] == item.2[2]) ) {
                tf.vec[1] <- TRUE
                return( tf.vec )
              }
              if( item.2[1] == item.2[2] ) {
                tf.vec[2] <- TRUE
                return( tf.vec )
              }
            }
            if((item.1[1] == item.2[1]) && (item.1[2] == item.2[2])) {
              tf.vec[3] <- TRUE
            }
            return( tf.vec )
          } )
        } )
      return( tf )
    }
    
    coef.a <- (4 * w.par^2 + 2)/( (2 * w.par + 1)^(0.5 * (p + 4)) ) -
      w.par^2 / (w.par + 1)^(p + 2)
    coef.b <- (4 * w.par^2)/( (2 * w.par + 1)^(0.5 * (p + 4)) ) -
      w.par^2 / (w.par + 1)^(p + 2)
    coef.c <- 1 / (2 * w.par + 1)^(0.5 * (p + 4))

    # Main routine
    idx.mtx <- as.matrix(bicond(bindex(1:p)))
    rslt <- 
      apply(idx.mtx, 1, function(vc) {
        sapply(vc, function(lst) {
          c(coef.a, coef.b, coef.c) %*% unlist(lst) }) })

    return(rslt)
  }
  
  ################################################################
  asymv.ve <- function(w.par, estcov) {

    ################################################################
    bicond <- function(bi) {
      tf <-
        sapply(bi, function(item.1) {
          lapply(bi, function(item.2) {
            tf.vec <- rep(FALSE, 2)
            item.1 <- sort(item.1); item.2 <- sort(item.2)
            if(item.1[1] == item.1[2] && item.2[1] == item.2[2]) {
              if(item.1[1] == item.2[1]) {
                tf.vec[1] <- TRUE
                return(tf.vec)
              } else {
                tf.vec[2] <- TRUE
                return(tf.vec)
              }
            }
            if((item.1[1] == item.2[1]) && (item.1[2] == item.2[2])) {
              tf.vec[2] <- TRUE
              return(tf.vec)
            }
            return(tf.vec) }) })
      return(tf)
    }
    
    ################################################################
    ve.core <- function(w.par, p) {
      coef.a1 <- -1.5 * w.par/(w.par + 1)^(0.5 * (p + 4))
      coef.a2 <- -0.5 * w.par/(w.par + 1)^(0.5 * (p + 4))
      coef.b <- 0.5 * w.par/(w.par + 1)^(0.5 * (p + 2))

      rslt <- 
        apply( as.matrix( bicond( bindex( 1:p ) ) ),
              1, function( vc ) {
                sapply( vc, function( lst ) {
                  c( coef.a1, coef.a2 ) %*% unlist( lst ) } ) } )

      vecp <- as.matrix(as.vector(diag(p)))
      rslt <- rslt + coef.b * (vecp %*% t(vecp))

      return( rslt )
    }

    ################################################################
    # Main routine
    p <- ncol(estcov)

    sqm.inv <- sqrtm(m = estcov)
    sqm <- solve(sqm.inv)
    
    first <- (sqm.inv %x% sqm.inv) %*%
      ve.core(w.par, p) %*% (sqm %x% sqm)

    coef.c <- 1/(w.par + 1)^(0.5 * (p + 2))
    second <- coef.c * diag(p * p)

    return(first + second)
  }
  
  ################################################################
  # Main routine
  p <- nrow(amat)
 
  sqm.inv <- sqrtm(m = estcov)

  dlt <- rob.delta(amat)
  kn <- sqm.inv %x% sqm.inv

  dp <- duplication(n = p)
  di <- duplication.inverse(n = p)

  coef.m <- t(dlt) %*% di %*% kn
  v <- coef.m %*% asymv.v(w.par = w.par, p = p) %*% t(coef.m)

  av.slv <-
    solve(t(dlt) %*% di %*% 
          asymv.ve(w.par = w.par, estcov) %*%
          dp %*% dlt)

  return(av.slv %*% v %*% t(av.slv))
}

