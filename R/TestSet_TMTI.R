#' Test a subset of hypotheses in its closure using the TMTI
#'
#' @param pvals Numeric vector of p-values.
#' @param subset Numeric vector; the subset to be tested.
#' @param alpha Numeric; the level to test at, if stopping early. Defaults
#' to 0.05.
#' @param tau Numeric; the treshold to use if using rTMTI. Set to NULL for TMTI
#' or rtTMTI. Defaults to NULL.
#' @param K Integer; The number of p-values to use if using rtTMTI. Set to NULL
#' for TMTI or tTMTI. Defaults to NULL.
#' @param EarlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' lower bounds on the p-values for the global test.
#' @param verbose Logical; set to TRUE to print progress.
#' @param gammaList List of functions. Must be such that the ith element
#' is the gamma function for sets of size i. Set to NULL to bootstrap the
#' functions assuming independence. Defaults to NULL.
#' @param mc.cores Number of cores to parallelize onto.
#' @param chunksize Integer indicating the size of chunks to parallelize. E.g.,
#' if setting chunksize = mc.cores, each time a parallel computation is set up,
#' each worker will perform only a single task. If mc.cores > chunksize, some
#' threads will be inactive.
#' @param is.sorted Logical, indicating the p-values are pre-sorted. Defaults
#' to FALSE.
#' @param ... Additional arguments.
#'
#' @return The adjusted p-value for the test of the hypothesis that there are
#' no false hypotheses among the selected subset.
#' @export
#'
#' @examples
#' ## Simulate p-values; 10 from false hypotheses, 10 from true
#' pvals = sort(c(
#'   rbeta(10, 1, 20), # Mean value of .1
#'   runif(10)
#' ))
#' ## Test whether the highest 10 contain any false hypotheses
#' TestSet_TMTI(pvals, subset = 11:20)

TestSet_TMTI = function(pvals,
                         subset,
                         alpha = 0.05,
                         tau = NULL,
                         K = NULL,
                         EarlyStop = FALSE,
                         verbose = FALSE,
                         gammaList = NULL,
                         mc.cores = 1L,
                         chunksize = 4 * mc.cores,
                         is.sorted = FALSE,
                         ...) {
  LocalTest = function (x) {
    TMTI::TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }
  TMTI::TestSet_LocalTest (
    LocalTest = LocalTest,
    pvals = pvals,
    subset = subset,
    alpha = alpha,
    is.sorted = is.sorted,
    EarlyStop = EarlyStop,
    verbose = verbose,
    mc.cores = mc.cores,
    chunksize = chunksize,
    ...
  )
}