#' Computes the TMTI test for a joint hypothesis given input p-values.
#'
#' @param pvals A vector of pvalues.
#' @param n A positive number (or Inf) indicating which type of local minimum
#' to consider. Defaults to Inf, corresponding to the global minimum.
#' @param tau Number between 0 and 1 or NULL, describing the truncation level.
#' @param K Integer between >1 and m describing the truncation index.
#' @param gamma Function; function to be used as the gamma approximation. If NULL, then
#' the gamma function will be bootstrapped assuming independence. Defaults
#' to NULL.
#' @param m_max Integer; the highest number of test for which the analytical
#' computation of the TMTI CDF is used. When m is above m_max it will be
#' bootstrapped or user supplied instead.
#' @param B Numeric; number of bootstrap replicates to be used when estimating
#' the gamma function. If a gamma is supplied, this argument is ignored.
#' Defaults to 1e3.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param ... Additional parameters.
#'
#' @return A p-value from the TMTI test
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' TMTI(pvals)

TMTI = function(pvals,
                 n = Inf,
                 tau = NULL,
                 K = NULL,
                 gamma = NULL,
                 B = 1e3,
                 m_max = 100,
                 is.sorted = FALSE,
                 ...) {
  if (!is.null(tau) & !is.null(K)) {
    stop("At most one of tau and K can be non NULL")
  }

  m = length(pvals)

  if (m == 1) {
    return(pvals)
  }

  if (is.sorted) {
    if (!is.null(tau)) {
      pvals = if (sum(pvals <= tau) > 0) pvals[pvals <= tau] else pvals[1]
    } else if (!is.null(K)) {
      pvals = pvals[1:K]
    }
  } else {
    if (!is.null(tau)) {
      pvals = if (sum(pvals <= tau) > 0) sort(pvals[pvals <= tau]) else min(pvals)
    } else if (!is.null(K)) {
      pvals = sort(pvals)[1:K]
    } else {
      pvals = pvals[order(pvals)]
    }
  }

  if (n >= m) {
    Z = MakeZ_C(pvals, m)
  } else {
    # Y = TMTI::MakeY_C(pvals = pvals, m)
    # Z = Y[.GetMinima(Y, n)]
    Z = MakeZ_C_nsmall(pvals, n, m)
  }



  if (!is.null(gamma)) {
    return(gamma(Z))
  } else if ((m <= m_max) & (n >= m)) {
    if (!is.null(K) & !is.null(tau)) {
      stop("Please supply only one of tau and K")
    } else if (!is.null(K)) {
      gamma = function(x) rtTMTI_CDF(x, m, K)
    } else if (!is.null(tau)) {
      gamma = function(x) tTMTI_CDF(x, m, tau)
    } else {
      gamma = function(x) TMTI_CDF(x, m)
    }
  } else if (is.null(gamma)) {
    gamma = gamma_bootstrapper(m, B = B, tau = tau, K = K, n = n, ...)
  }

  gamma(Z)
}
