% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vc_test.R
\name{vc_test}
\alias{vc_test}
\title{Perform variance component test}
\usage{
vc_test(eta_prime, se, omega_inv_sq, ce, tol = 1e-08)
}
\arguments{
\item{eta_prime}{A vector giving the transformed effect size}

\item{se}{A vector of the effect size standard errors}

\item{omega_inv_sq}{A matrix of square root inverse of genetic correlation}

\item{ce}{A matrix of the genetic correlation}
}
\value{
none
}
\description{
Perform variance component test
}
