#' TopDown TMTI algorithm for estimating a 1-alpha confidence set for the number
#' of false hypotheses among a set.
#'
#' @param pvals A vector of p-values.
#' @param subset Numeric vector specifying a subset a p-values to estimate a
#' confidence set for the number of false hypotheses for. Defaults to NULL
#' corresponding to estimating a confidence set for the number of false
#' hypotheses in the entire set.
#' @param alpha Level in [0,1] at which to generate confidence set. Defaults
#' to 0.05.
#' @param gammaList List of pre-specified gamma functions. If NULL, the
#' functions will be approximated by bootstrap assuming independence. Defaults
#' to NULL.
#' @param verbose Logical, indicating whether or not to write out the progress.
#' Defaults to TRUE.
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' is.sorted. Defaults to FALSE.
#' @param mc.cores Number of cores to parallelize onto.
#' @param chunksize Integer indicating the size of chunks to parallelize. E.g.,
#' if setting chunksize = mc.cores, each time a parallel computation is set up,
#' each worker will perform only a single task. If mc.cores > chunksize, some
#' threads will be inactive.
#' @param direction A string indicating whether to perform a binary search ('binary'/'b')
#' or decreasing ('decreasing'/'d') search. Defaults to 'binary', which has better
#' computational complexity.
#' @param ... Additional parameters.
#'
#' @return A 1-alpha lower bound for the number of false hypotheses among the
#' set of supplied p-values
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' TopDown_TMTI(pvals)
#'

TopDown_TMTI = function(pvals,
                        subset = NULL,
                        alpha = 0.05,
                        gammaList = NULL,
                        verbose = TRUE,
                        tau = NULL, K = NULL,
                        is.sorted = FALSE,
                        mc.cores = 1L,
                        chunksize = 4 * mc.cores,
                        direction = 'binary',
                        ...) {
  LocalTest = function (x) {
    TMTI::TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }
  TMTI::TopDown_LocalTest (
    LocalTest = LocalTest,
    pvals = pvals,
    subset = subset,
    alpha = alpha,
    is.sorted = is.sorted,
    verbose = verbose,
    mc.cores = mc.cores,
    chunksize = chunksize,
    direction = direction,
    ...
  )
}

