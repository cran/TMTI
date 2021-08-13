#' Test a subset of hypotheses in its closure using a user-specified local test
#'
#' @param localTest Function which defines a combination test.
#' @param pvals Numeric vector of p-values
#' @param subset Numeric vector; the subset to be tested
#' @param alpha Numeric; the level to test at, if stopping early. Defaults
#' to 0.05
#' @param earlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' lower bounds on the p-values for the global test.
#' @param verbose Logical; set to TRUE to print progress.
#' @param ... Additional arguments to be passed onto TMTI()
#'
#' @return The adjusted p-value for the test of the hypothesis that there are
#' no false hypotheses among the selected subset.
#' @export
#'
#' @examples
#' ## Simulate p-values; 10 from false hypotheses, 10 from true
#' pvals <- sort(c (
#'     rbeta(10, 1, 20),  # Mean value of .1
#'     runif(10)
#' ))
#' ## Test whether the highest 10 contain any false hypotheses using a Bonferroni test
#' TestSet_localTest(function(x) {min(c(1, length(x) * min(x)))}, pvals, subset = 11:20)
TestSet_localTest <- function (
  localTest,
  pvals,
  subset,
  alpha = 0.05,
  earlyStop = FALSE,
  verbose   = FALSE,
  ...
) {
  m     <- length(pvals)
  m2    <- length(subset)
  pSub  <- sort(pvals[subset])
  pRest <- sort(pvals[-subset])

  out <- list()

  pfirst <- localTest(pSub)
  out[[1]] <- c(
    "p"     = pfirst,
    "layer" = 0,
    "Accept" = (pfirst >= alpha)
  )

  if (earlyStop & out[[1]][3]) {
    return(out[[1]][1])
  }

  stepCounter <- 0

  for (i in length(pRest):1) {
    stepCounter <- stepCounter + 1
    if (verbose) {
      cat("\rStep", stepCounter, " of ", length(pRest))
    }
    ptilde <- c(pSub, pRest[length(pRest):i])
    pp <- localTest (
      ptilde
    )
    out[[stepCounter + 1]] <- c (
      "p" = pp,
      "layer" = stepCounter,
      "Accept" = (pp >= alpha)
    )

    if (earlyStop & out[[stepCounter + 1]][3])
      break
  }

  out <- do.call("rbind", out)
  max(out[, 1])
}
