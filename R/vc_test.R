#' Perform variance component test
#' @param eta_prime A vector giving the transformed effect size
#' @param se A vector of the effect size standard errors
#' @param omega_inv_sq A matrix of square root inverse of genetic correlation
#' @param ce A matrix of the genetic correlation
#' @return none
#' @export
vc_test <- function(eta_prime, se, omega_inv_sq, ce, tol = 1e-8){
  n <- length(eta_prime)
  eig <- eigen(t(omega_inv_sq) %*% diag(se) %*% ce %*% diag(se) %*% omega_inv_sq,
               symmetric = T)
  eigval <- eig$values[n:1]; eigvec <- -eig$vectors[,n:1]
  delta2 <- (t(eigvec[, eigval > tol]) %*% eta_prime)^2
  Delta <- eigval[eigval > tol]
  return(vc_optim(n, delta2, Delta))
}

#' @keywords internal
vc_optim <- function(n, d1, d2){
  it <- 10^((-36:23)/4)
  init <- it[which.max(sapply(it, function(i) llike(i, n, d1, d2)))]
  tausq <- max(0, newton_optim(init, d1, d2))
  ll0 <- llike(0, n, d1, d2)
  ll1 <- llike(tausq, n, d1, d2)
  if(ll1 < ll0){
    tausq <- 0
    ll1 <- ll0
  }
  stat <- -2 * (ll0 - ll1)
  return(c(tausq, stat))  
}

#' @keywords internal
newton_optim <- function(x, delta2, Delta, maxit = 1e4, tol = 1e-6){
  i <- 0
  while(abs(dllike(x, delta2, Delta)) > tol & i < maxit){
    x <- x - dllike(x, delta2, Delta) / d2llike(x, delta2, Delta)
    i <- i + 1
  }
  return(x)
}

#' @keywords internal
llike <- function(x, n, delta2, Delta){
  return(-1/2 * (n * log(2 * pi) + sum(log(Delta + x)) + sum(delta2 / (Delta +  x))))
}

#' @keywords internal
dllike <- function(x, delta2, Delta){
  return(-0.5 * (sum(1 / (Delta + x)) - sum(delta2 / (Delta + x)^2)))
}

#' @keywords internal
d2llike <- function(x, delta2, Delta){
  return(-0.5 * (-sum(1 / (Delta + x)^2) + 2 * sum(delta2 / (Delta + x)^3)))
}

