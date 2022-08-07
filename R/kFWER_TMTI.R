#' kFWER_TMTI. Computes the largest rejection set possible with kFWER control.
#'
#' @param pvals A vector p-values.
#' @param k An integer denoting the desired k at which to control the kFWER.
#' @param alpha Significance level.
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
#' @param verbose Logical, indicating whether or not to print progress.
#'
#' @return The number of marginal hypotheses that can be rejected with kFWER control.
#' @export
#'
#' @examples
#' nfalse = 50
#' m = 100
#' pvals = c (
#'   sort(runif(nfalse, 0, 0.05 / m)),
#'   sort(runif(m - nfalse, 0.1, 1))
#' )
#' kFWER_TMTI (
#'   pvals = pvals,
#'   k = 5,
#'   alpha = 0.05,
#'   verbose = FALSE
#' )
kFWER_TMTI = function (
    pvals,
    k,
    alpha = 0.05,
    B = 1e3,
    gammaList = NULL,
    tau = NULL,
    K = NULL,
    verbose = FALSE
) {
  if(is.unsorted(pvals))
    pvals = sort(pvals)
  LocalTest = function (x) {
    TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }
  if (k <= 1L) {
    return (
      FWER_set_C(
        LocalTest = LocalTest,
        pvals = pvals,
        alpha = alpha,
        low = 0,
        high = length(pvals) - 1,
        verbose = verbose
      )
    )
  } else {
    return (
      kFWER_set_C(
        LocalTest = LocalTest,
        pvals = pvals,
        k = k + 1,
        alpha = alpha,
        low = 0,
        high = length(pvals) - 1,
        verbose = verbose
      )
    )
  }
}
