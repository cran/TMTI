#' Test a subset of hypotheses in its closure using the TMTI
#'
#' @param pvals Numeric vector of p-values
#' @param subset Numeric vector; the subset to be tested
#' @param alpha Numeric; the level to test at, if stopping early. Defaults
#' to 0.05
#' @param tau Numeric; the treshhold to use if using rTMTI. Set to NULL for TMTI
#' or rtTMTI. Defaults to NULL
#' @param K Integer; The number of p-values to use if using rtTMTI. Set to NULL
#' for TMTI or tTMTI. Defaults to NULL.
#' @param earlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' lower bounds on the p-values for the global test.
#' @param verbose Logical; set to TRUE to print progress.
#' @param gammalist List of functions. Must be such that the i'th element
#' is the gamma function for sets of size i. Set to NULL to bootstrap the
#' functions assuming independence. Defaults to NULL.
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
#' ## Test whether the highest 10 contain any false hypotheses
#' TestSet_TMTI(pvals, subset = 11:20)
TestSet_TMTI <- function (
  pvals,
  subset,
  alpha = 0.05,
  tau   = NULL,
  K     = NULL,
  earlyStop = FALSE,
  verbose   = FALSE,
  gammalist = NULL,
  ...
) {
  m     <- length(pvals)
  m2    <- length(subset)
  pSub  <- sort(pvals[subset])
  pRest <- sort(pvals[-subset])
  # pRest <- sort(pvals[pvals > max(pSub)])

  out <- list()

  if (!is.null(K) & length(K) < m) {
    K <- rep(K, length.out = m)
  }

  pfirst <- TMTI(pSub, tau = tau, K = K[m2], gamma = gammalist[[m2]], ...)
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
    pp <- TMTI (
      ptilde,
      tau = tau,
      K   = K[length(ptilde)],
      gamma = gammalist[[length(ptilde)]],
      ...
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
