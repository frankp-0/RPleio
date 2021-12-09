#' Perform importance sampling
#' @param N An integer, number of simulations to perform
#' @param se A vector of effect size standard errors
#' @param omega_inv_sq A matrix of square root inverse of genetic correlation
#' @param ce A matrix of the genetic correlation
#' @param out A string giving the output file
#' @param ncores An integer of the number of cores
#' @return none
#' @export
is_estim <- function(N, se, omega_inv_sq, ce, out, ncores = 1){
  n <- length(se)
  D <- diag(se) %*% ce %*% diag(se)
  D_0 <- diag(rep(1, n)) %*% ce %*% diag(rep(1, n))
  c_j <- c(1,1.1,1.2,1.3,1.4,1.7,2,2.5,3,4,5)
  n_c <- length(c_j)
  alpha <- rep(1 / n_c, n_c)
  Pf <- c_sample(N, alpha, D_0, c_j)
  eta_prime <- sweep(Pf, 2, se, `*`) %*% t(omega_inv_sq)
  data <- mclapply(1:nrow(eta_prime), function(i) {
    b <- unlist(eta_prime[i,])
    vc_test(b, se, omega_inv_sq, ce)[2]
  }, mc.cores = ncores) %>% unlist()
  
  qf <- mvtnorm::dmvnorm(Pf, rep(0, n), ce)
  pf <- lapply(1:n_c, function(i) {
    mvtnorm::dmvnorm(Pf, rep(0, n), D_0 * c_j[[i]]^2)
  })
  thetas <- c(40^(seq(-5,0,length=20))[1:19], 40^(seq(0,1,length=20)))
  p_alpha <- rowSums(t(do.call("rbind", pf)) * alpha)
  p <- sapply(thetas, function(theta) p_estim(data, theta, p_alpha, alpha, qf, pf, n_c, N))
  is <- data.frame(theta = thetas, p = p)
  is <- is[order(is$p, decreasing = T),]
  is <- rbind(c(0, 1), is)
  fwrite(is, file = out)
}

#' @keywords internal
p_estim <- function(data, theta, p_alpha, alpha, qf, pf, n_c, N){
  data <- as.numeric(data > theta)
  m <- data * qf / p_alpha
  vm <- sapply(1:length(alpha), function(i) cov(cbind(m, pf[[i]]))[1,2])
  v <- cov(t(do.call("rbind", pf))/p_alpha)
  s <- svd(v)
  inv_v <- s$v[,-ncol(s$v)] %*% diag(1 / s$d[-length(s$d)]) %*% t(s$u[,-ncol(s$u)])
  betas <- (inv_v %*% vm)[,1]
  control_variate <- rowSums(t(do.call("rbind", pf) * betas))
  p <- sum((data * qf - control_variate) / rowSums(t(do.call("rbind", pf) * alpha))) / N + sum(betas)
  return(p)
}

#' @keywords internal
c_sample <- function(N, alpha, D_0, c_j){
  counts <- rep(0, length(alpha))
  comp <- sample(1:length(alpha), N, replace = T, prob = alpha)
  counts[as.integer(names(table(comp)))] <- table(comp)
  Pf <- lapply(which(counts != 0), function(i) {
    MASS::mvrnorm(n = counts[i], mu = rep(0, nrow(D_0)), Sigma = D_0 * c_j[[i]]^2)
  })
  return(do.call("rbind", Pf))
}
