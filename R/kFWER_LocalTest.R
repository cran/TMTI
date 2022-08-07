#' kFWER_LocalTest. Computes the largest rejection set possible with kFWER control.
#'
#' @param LocalTest A function that returns a p-value for a joint hypothesis test.
#' @param pvals A vector p-values.
#' @param k An integer denoting the desired k at which to control the kFWER.
#' @param alpha Significance level.
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
#' kFWER_LocalTest (
#'   LocalTest = function (x) min(x) * length(x),
#'   pvals = pvals,
#'   k = 5,
#'   alpha = 0.05,
#'   verbose = FALSE
#' )
kFWER_LocalTest = function (
    LocalTest,
    pvals,
    k,
    alpha = 0.05,
    verbose = FALSE
) {
  if(is.unsorted(pvals))
    pvals = sort(pvals)
  if (k <= 0L) {
    stop("k must be greater than or equal to 1")
  } else if (k <= 1L) {
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
        # low = k - 1,
        high = length(pvals),# - 1,
        verbose = verbose
      )
    )
  }


}
