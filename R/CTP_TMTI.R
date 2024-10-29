#' A Closed Testing Procedure for the TMTI using an O(n^2) shortcut
#'
#' @name CTP_TMTI
#' @aliases TMTI_CTP
#' @param pvals A vector of p-values.
#' @param alpha Level to perform each intersection test at. Defaults to 0.05.
#' @param B Number of bootstrap replications if gamma needs to be approximated.
#' Not used if specifying a list of functions using the gammaList argument
#' or if length(pvals) <= 100. Defaults to 1000.
#' @param gammaList A list of pre-specified gamma functions. If NULL, gamma
#' functions will be approximated via bootstrap, assuming independence. Defaults
#' to NULL.
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param is.sorted Logical, indicating the p-values are pre-sorted. Defaults
#' to FALSE.
#' @param EarlyStop Logical indicating whether to exit as soon as a non-significant
#' p-value is found. Defaults to FALSE.
#' @param ... Additional arguments.
#'
#' @return A data.frame containing adjusted p-values and the original index of
#' the p-values.
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' CTP_TMTI(pvals)
#'
CTP_TMTI = function(pvals, alpha = 0.05, B = 1e3,
                    gammaList = NULL,
                    tau = NULL, K = NULL,
                    is.sorted = FALSE,
                    EarlyStop = FALSE,
                    ...) {
  LocalTest = function (x) {
    TMTI::TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }

  CTP_LocalTest (
    LocalTest = LocalTest,
    pvals = pvals,
    alpha = alpha,
    is.sorted = is.sorted,
    EarlyStop = EarlyStop,
    ...
  )
}

#'
#' @rdname CTP_TMTI
#' @export
TMTI_CTP = function(pvals, alpha = 0.05, B = 1e3,
                    gammaList = NULL,
                    tau = NULL, K = NULL,
                    is.sorted = FALSE,
                    ...) {
  .Deprecated(new = "CTP_TMTI")

  CTP_TMTI (
    pvals = pvals,
    alpha = alpha,
    B = B,
    gammaList = gammaList,
    tau = tau,
    K = K,
    is.sorted = is.sorted,
    ...
  )
}
