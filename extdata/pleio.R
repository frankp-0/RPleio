#!/usr/bin/env Rscript

# Setup ==============================
library(optparse)
library(data.table)
library(parallel)
library(magrittr)
library(RPleio)

option_list <- list(
  make_option(c("-s", "--summary"),
              type = "character",
              help = "formatted summary statistics file"
              ),
  make_option(c("-g", "--cg"),
              type = "character",
              help = "genetic covariance matrix file"
              ),
  make_option(c("-e", "--ce"),
              type = "character",
              help = "environmental covariance matrix file"
              ),
  make_option(c("-n", "--nis"),
              type = "integer",
              help = "numer of simulations for importance sampling"
              ),
  make_option(c("-c", "--cores"),
              type = "integer",
              default = detectCores() - 2,
              help = "number of cores to use for parallelization"
              ),
  make_option(c("-f", "--isf"),
              default = "",
              help = "importance sampling output file"
              ),
  make_option(c("-o", "--out"),
              type = "character",
              help = "output file prefix"
              )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Read data
sumstat <- fread(opt$summary)[,-1]
cg <- as.matrix(fread(opt$cg))
ce <- as.matrix(fread(opt$ce))
nis <- opt$nis
out_pre <- opt$out
ncores <- opt$cores
isf <- opt$isf

# Importance Sampling ==============================
n <- ncol(sumstat) / 2
#' @keywords internal
sqrt_ginv <- function(X, tol = sqrt(1e-8)){
  Xsvd <- svd(X);
  u <- Xsvd$u; d <- Xsvd$d; v <- t(Xsvd$v)    
  pos <- d > max(tol * d[1], 0)
  res <- t(v)[,pos] %*% sqrt(diag(1 / d[pos])) %*% t(u[, pos])
  return(res)
}

omega_inv_sq <- sqrt_ginv(cg)
if (isf == "") {
  ind <- seq(2, 2 * n, 2)
  se <- 1 / sqrt(colMeans(1 / sumstat[, ..ind]^2))
  isf <- paste0(out_pre, ".is")
  is_estim(nis, se, omega_inv_sq, ce, isf, ncores)
}

# Variance Component Test ==============================
ind <- seq(1, 2 * n, 2)

df <- mclapply(1:nrow(sumstat), function(i) {
  se <- unlist(sumstat[i, .SD, .SDcols = ind + 1])
  eta_prime <- t(omega_inv_sq) %*% unlist(sumstat[i, ..ind])
  vc_test(eta_prime, se, omega_inv_sq, ce)},
  mc.cores = ncores
  )
df <- data.table(do.call("rbind", df))
names(df) <- c("tausq", "pleio_stat")

# P-values
is <- fread(isf)
ss <- smooth.spline(x = is$theta, y = log(is$p))

df[, pleio_p := exp(predict(ss, pleio_stat)$y)]
b <- cov(is$theta, log(is$p)) / var(is$theta)
a <- log(is$p)[nrow(is)] - b * is$theta[nrow(is)]
df[pleio_stat <= is$theta[2], pleio_p := 1]
df[pleio_stat > is$theta[40], pleio_p := exp(a + b * pleio_stat)]


fwrite(df, paste0(out_pre, ".vc"))
