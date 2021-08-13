#' Functions for computing the TMTI tests
#'
#' @param pvals A vector of pvalues
#' @param n A positive number (or Inf) indicating which type of local minimum
#' to consider. Defaults to Infm, corresponding to the global minimum.
#' @param tau Number between 0 and 1 or NULL, describing the truncation level.
#' @param K Integer between >1 and m describing the truncation index.
#' @param gamma Function; function to be used as the gamma approximation. If NULL, then
#' the gamma function will be bootstrapped assuming independence. Defaults
#' to NULL
#' @param m_max Integer; the highest number of test for which the analytical
#' computation of the TMTI CDF is used. When m is above m_max it will be
#' bootstrapped or user supplied instead.
#' @param log.p Logical; indicating whether to compute Y's on log-scale.
#' Defaults to TRUE
#' @param B Numeric; number of bootstrap replicates to be used when estimating
#' the gamma function. If a gamma is supplied, this argument is ignored.
#' Defaults to 1e3.
#' @param ... Additional parameters
#'
#' @return A p-value from the TMTI test
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'     rbeta(10, 1, 20),  ## Mean value of .05
#'     runif(10)
#' )
#' TMTI(pvals)
TMTI <- function (
  pvals,
  n     = Inf,
  tau   = NULL,
  K     = NULL,
  gamma = NULL,
  log.p = TRUE,
  B = 1e3,
  m_max = 100,
  ...
) {
  m <- length(pvals)

  if (m == 1) {
    return (pvals)
  }

  Y <- make_Y(pvals, tau = tau, K = K, log.p = log.p)
  Z <- Y[.GetMinima(Y, n)]

  if (m <= m_max & log.p) {
    if (!is.null(K) & !is.null(tau)) stop("Please supply only one of tau and K")
    else if (!is.null(K)) gamma <- function (x) rtTMTI_CDF(exp(x), m, K)
    else if (!is.null(tau)) gamma <- function (x) tTMTI_CDF(exp(x), m, tau)
    else gamma <- function (x) TMTI_CDF(exp(x), m)
  } else if (m <= 100 & !log.p) {
    if (!is.null(K) & !is.null(tau)) stop("Please supply only one of tau and K")
    else if (!is.null(K)) gamma <- function (x) rtTMTI_CDF(x, m, K)
    else if (!is.null(tau)) gamma <- function (x) tTMTI_CDF(x, m, tau)
    else gamma <- function (x) TMTI_CDF(x, m)
  } else if(is.null(gamma)) {
    gamma <- gamma_bootstrapper(m, log.p = log.p, B = B, tau = tau, K = K, ...)
    # gamma_bootstrapper(Z = Z, m = m, n = n, log.p = log.p, ...)
  }
  # else {
  #   gamma(Y[.GetMinima(Y, n)])
  # }

  gamma(Z)
}

tTMTI <- function (
  pvals,
  tau,
  n = Inf,
  gamma = NULL,
  ...
) {
  TMTI(pvals, n = n, tau = tau, K = NULL, gamma = gamma, ...)
}

rtTMTI <- function (
  pvals,
  K,
  n = Inf,
  gamma = NULL,
  ...
) {
  TMTI(pvals, n = n, tau = NULL, K = K, gamma = gamma, ...)
}
