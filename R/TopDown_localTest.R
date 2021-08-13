#' TopDown localTest algorithm for estimating a 1-alpha confidence set for the number
#' of false hypotheses among a set.
#'
#' @param localTest A function specifying a local test.
#' @param pvals A vector of p-values
#' @param subset Numeric vector specifying a subset a p-values to estimate a
#' confidence set for the number of false hypotheses for. Defaults to NULL
#' corresponding to estimating a confidence set for the number of false
#' hypotheses in the entire set.
#' @param alpha Level in [0,1] at which to generate confidence set. Defaults
#' to 0.05
#' @param verbose Logical, indicating whether or not to write out the progress.
#' Defaults to TRUE
#' @param ... Additional parameters
#'
#' @return A lower 1-alpha bound for the number of false hypotheses among the
#' set of supplied p-values
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals <- c (
#'   rbeta(10, 1, 20),  ## Mean value of .05
#'   runif(10)
#' )
#' ## Estimate the confidence set using a local Bonferroni test
#' TopDown_localTest(function(x) {min(c(1, length(x) * min(x)))}, pvals)

TopDown_localTest <- function (
  localTest,
  pvals,
  subset = NULL,
  alpha = 0.05,
  verbose = TRUE,
  ...
) {
  ord     <- order(pvals)
  pvals   <- sort(pvals)
  m       <- length(pvals)
  t_alpha <- 0

  if (!is.null(subset) & length(subset) < length(pvals)) {
    counter <- 0
    top <- length(subset)
    for (i in subset) {
      counter <- counter + 1
      if (verbose) {
        cat("\rStep", counter, " of ", length(subset))
      }
      subset2 <- subset[length(subset)]:i
      p_loc <- TestSet_localTest (
        localTest,
        pvals,
        subset2,
        alpha = alpha,
        earlyStop = TRUE,
        ...
      )
      accept <- (p_loc >= alpha)
      if (accept) {
        t_alpha <- length(subset2)
        break
      }
    }
  } else {
    top <- length(pvals)
    for (i in m:1) {
      if (verbose) cat("\rStep", i)
      pvals_tilde <- pvals[(m - i + 1):m]
      p_loc <- localTest(pvals_tilde)
      accept <- (p_loc >= alpha)
      if (accept) {
        t_alpha <- i
        break
      }
    }
  }

  cat (
    paste0 (
      "Confidence set for the number of false hypotheses is {",
      top - t_alpha,
      ",..., ",
      top,
      "}\n"
    )
  )
  return(top - t_alpha)
}
