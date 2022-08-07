#' A Closed Testing Procedure for any local test satisfying the conditions of Mogensen and Markussen (2021) using an O(n^2) shortcut.
#'
#' @name CTP_LocalTest
#' @aliases localTest_CTP
#' @param LocalTest A function which defines the choice of local test to use.
#' @param pvals A vector of p-values.
#' @param alpha Level to perform each intersection test at. Defaults to 0.05.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
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
#' ## Perform the CTP using a local Bonferroni test
#' CTP_LocalTest(function(x) {
#'   min(c(length(x) * min(x), 1))
#' }, pvals)
#'

CTP_LocalTest = function(
    LocalTest,
    pvals,
    alpha = 0.05,
    is.sorted = FALSE,
    EarlyStop = FALSE,
    ...
) {
  if (is.sorted) {
    ord = 1:length(pvals)
  } else {
    ord = order(pvals)
    pvals = sort(pvals)
  }
  f = function (x, y) {
    TestSet_C (
      LocalTest = LocalTest,
      pSub = x,
      pRest = y,
      alpha = 0.05,
      is_subset_sequence = TRUE,
      EarlyStop = EarlyStop,
      verbose = FALSE
    )
  }

  p_adjusted = FullCTP_C (
    LocalTest = LocalTest,
    f = f,
    pvals = pvals,
    EarlyStop = EarlyStop,
    alpha = alpha
  )
  data.frame (
    "p_adjusted" = p_adjusted,
    "index"      = ord[1:length(p_adjusted)]
  )
}

#'
#' @rdname CTP_LocalTest
#' @param localTest A function specifying a local test (deprecated).
#' @export

localTest_CTP = function(localTest, pvals, alpha = 0.05, is.sorted = FALSE, ...) {
  if (is.sorted) {
    ord = 1:length(pvals)
  } else {
    ord = order(pvals)
    pvals = sort(pvals)
  }
  .Deprecated(new = "CTP_LocalTest")

  CTP_LocalTest (
    LocalTest = localTest,
    pvals = pvals,
    alpha = alpha,
    is.sorted = is.sorted,
    ...
  )
}
