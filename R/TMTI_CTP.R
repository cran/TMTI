#' A Closed Testing Procedure for the TMTI using an O(n^2) shortcut
#'
#' @param pvals A vector of p-values
#' @param alpha Level to perform each intersection test at. Defaults to 0.05
#' @param B Number of bootstrap replications if gamma needs to be approximated.
#' Not used if specifying a list of functions using the gammaList argument
#' or if length(pvals) <= 100. Defaults to 1000
#' @param gammaList A list of pre-specified gamma functions. If NULL, gamma
#' functions will be approximated via bootstrap, assuming independence. Defaults
#' to NULL.
#' @param log.p Logical, indicating whether to compute Y's on log-scale.
#' Defaults to TRUE
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param ... Additional arguments
#'
#' @return A data.frame containing:
#' * i: The sorted index of each p-value.
#' * p_adjust: The CTP adjusted p-value, controlling the FWER strongly.
#' * FirstAccept: The first level of the test tree at which the hypothesis could
#' not be rejected. NA if it is never rejected.
#' * Index: The original index of the unsorted p-value inputs.
#' @export TMTI_CTP
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' TMTI_CTP(pvals, earlyStop = TRUE)


TMTI_CTP <- function (
  pvals, alpha = 0.05, B = 1e3,
  gammaList = NULL,
  log.p = TRUE,
  tau = NULL, K = NULL,
  ...
) {
  ord <- order(pvals)
  p2  <- pvals
  p   <- sort(pvals)
  m   <- length(pvals)

  if (!is.null(K)) {
    if (length(K) < length(pvals))
      K <- rep(K, length.out = length(pvals))
  }

  Q <- matrix(0, m, m)

  Q[m, m] <- p[m]
  for (i in 1:(m - 1)) {
    counter <- m

    Q[counter, i] <- p[i]

    for (j in m:(i + 1)) {
      counter <- counter - 1

      subp <- p[c(i, m:j)]
      m2   <- length(subp)

      Q[counter, i] <- TMTI (
        pvals = subp,
        tau   = tau,
        K     = K[m2],
        log.p = log.p,
        gamma = gammaList[[m2]]
      )
    }
  }
  for (i in 1:(m - 1)) {
    Q[i, (i + 1):m] <- diag(Q)[i]
  }
  adjusted_p <- apply(Q, 2, max)

  data.frame (
    "p_adjusted" = adjusted_p,
    "index"      = ord
  )
}
