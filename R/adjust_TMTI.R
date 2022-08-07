#' Adjust all p-values using a Closed Testing Procedeure and the
#' TMTI family of tests.
#'
#' @param pvals vector of p-values.
#' @param alpha significance level. Defaults to 0.05.
#' @param B Number of bootstrap replications. Only relevant if length(pvals) > 100
#' and no gammaList is supplied.
#' @param gammaList A list of functions. These functions should be the CDFs of
#' the chosen TMTI test for different m.
#' @param tau Number between 0 and 1 or NULL, describing the truncation level.
#' @param K Integer between >1 and m describing the truncation index.
#' @param is.sorted Logical, indicating whether the supplied p-values are already
#' sorted. Defaults to FALSE.
#' @param EarlyStop Logical; set to TRUE to stop as soon as a hypothesis can be
#' accepted at level alpha. This speeds up the procedure, but now only provides
#' upper bounds on the adjusted p-values that are below alpha.
#' @param verbose Logical; set to TRUE to print progress. Defaults to FALSE.
#' @param mc.cores Number of cores to parallelize onto.
#' @param chunksize Integer indicating the size of chunks to parallelize. E.g.,
#' if setting chunksize = mc.cores, each time a parallel computation is set up,
#' each worker will perform only a single task. If mc.cores > chunksize, some
#' threads will be inactive.
#' @param direction String that is equal to either "increasing"/"i", "decreasing"/"d" or "binary"/"b".
#' Determines the search direction. When set to"increasing", the function computes the exact adjusted p-value
#' for all those hypotheses that can be rejected (while controlling the FWER),
#' but is potentially slower than "decreasing". "decreasing"identifies all hypotheses that can
#' be rejected with FWER control, but does not compute the actual adjusted p-values.
#' "binary" performs a binary search for the number of hypotheses
#' that can be rejected with FWER control.  Defaults to "increasing". Note that
#' 'binary' does not work with parallel.direction == 'breadth'.
#' @param parallel.direction A string that is either "breadth" or "depth"
#' (or abbreviated to "b" or "d), indicating in which direction to parallelize.
#' Breadth-first parallelization uses a more efficient C++ implementation to
#' adjust each p-value, but depth-first parallelization potentially finishes
#' faster if using early stopping (EarlyStop = TRUE) and very few hypotheses
#' can be rejected.
#' @param AdjustAll Logical, indicating whether to adjust all p-values (TRUE)
#' or only those that are marginally significant (FALSE). Defaults to FALSE.
#' @param ... Additional arguments.
#'
#' @return a data.frame containing adjusted p-values and their respective indices.
#' If direction == 'decreasing' or 'binary', an integer describing the number of
#' hypotheses that can be rejected with FWER control is returned.
#' @export
#'
#' @examples
#' p = sort(runif(100)) # Simulate and sort p-values
#' p[1:10] = p[1:10]**3 # Make the bottom 10 smaller, such that they correspond to false hypotheses
#' adjust_TMTI(p, alpha = 0.05, is.sorted = TRUE)
adjust_TMTI = function(pvals, alpha = 0.05, B = 1e3,
                       gammaList = NULL,
                       tau = NULL, K = NULL,
                       is.sorted = FALSE,
                       EarlyStop = FALSE,
                       verbose = FALSE,
                       mc.cores = 1L,
                       chunksize = 4 * mc.cores,
                       direction = "increasing",
                       parallel.direction = "breadth",
                       AdjustAll = FALSE,
                       ...) {
  LocalTest = function (x) {
    TMTI::TMTI(x, tau = tau, K = K, gamma = gammaList[[length(x)]])
  }
  TMTI::adjust_LocalTest (
    LocalTest = LocalTest,
    pvals = pvals,
    alpha = alpha,
    is.sorted = is.sorted,
    EarlyStop = EarlyStop,
    verbose = verbose,
    mc.cores = mc.cores,
    chunksize = chunksize,
    direction = direction,
    parallel.direction = parallel.direction,
    AdjustAll = AdjustAll,
    ...
  )

}
