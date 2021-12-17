#' Perform RPleio
#' @param sum_file A string giving the summary statistic file
#' @param g_file A string giving the genetic correlation matrix file
#' @param e_file A string giving the environmental correlation matrix file
#' @param nis A integer giving the number of samples to be used for importance sampling
#' @param ncores An integer giing the number of cores to use for parallelization
#' @param isf A string giving the importance sampling file. Optional
#' @param out A string giving the output prefix
#' @return none
#' @import data.table
#' @import magrittr
#' @import stats
#' @import parallel
#' @export
pleio <- function(sum_file,
                  g_file,
                  e_file,
                  nis,
                  ncores,
                  isf = "",
                  out){
  
  # Read data ==============================
  sumstat <- data.table::fread(sum_file)[,-1]
  cg <- as.matrix(data.table::fread(g_file))
  ce <- as.matrix(data.table::fread(e_file))

  # Importance sampling ====================
  omega_inv_sq <- sqrt_ginv(cg)
  n <- ncol(sumstat) / 2
  if (isf == "") {
    ind <- seq(2, 2 * n, 2)
    se <- 1 / sqrt(colMeans(1 / sumstat[, ..ind]^2))
    isf <- paste0(out, ".is")
    print("Performing importance sampling...")
    is_estim(nis, se, omega_inv_sq, ce, isf, ncores)
  }

  # Variance Component Test ==============================
  ind <- seq(1, 2 * n, 2)
  print("Performing variance component test...")
  df <- mclapply(1:nrow(sumstat), function(i) {
    se <- unlist(sumstat[i, .SD, .SDcols = ind + 1])
    eta_prime <- t(omega_inv_sq) %*% unlist(sumstat[i, ..ind])
    vc_test(eta_prime, se, omega_inv_sq, ce)},
    mc.cores = ncores
    )
  df <- data.table(do.call("rbind", df))
  names(df) <- c("tausq", "pleio_stat")

  # P-values =============================================
  is <- fread(isf)
  ss <- smooth.spline(x = is$theta, y = log(is$p))

  df[, pleio_p := exp(predict(ss, pleio_stat)$y)]
  b <- cov(is$theta, log(is$p)) / var(is$theta)
  a <- log(is$p)[nrow(is)] - b * is$theta[nrow(is)]
  df[pleio_stat <= is$theta[2], pleio_p := 1]
  df[pleio_stat > is$theta[40], pleio_p := exp(a + b * pleio_stat)]

  # Write output ========================================
  print("Writing output file...")
  fwrite(df, paste0(out, ".txt.gz"), compress = "gzip")
}

#' @keywords internal
sqrt_ginv <- function(X, tol = sqrt(1e-8)){
  Xsvd <- svd(X);
  u <- Xsvd$u; d <- Xsvd$d; v <- t(Xsvd$v)    
  pos <- d > max(tol * d[1], 0)
  res <- t(v)[,pos] %*% sqrt(diag(1 / d[pos])) %*% t(u[, pos])
  return(res)
}
