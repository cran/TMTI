#' Function to bootstrap the Cumulative Distribution Functions (CDFs) of the TMTI statistics.
#'
#' @param m Number of tests.
#' @param n Number (or Inf) indicating what kind of minimum to consider.
#' Defaults to Inf, corresponding to the global minimum.
#' @param B Number of bootstrap replicates. Rule of thumb is to use at least
#' 10 * m.
#' @param mc.cores Integer denoting the number of cores to use when using
#' parallelization, Defaults to 1, corresponding to single-threaded computations.
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#'
#' @return An approximation of the function \eqn{\gamma^m(x)} under the
#' assumption that all p-values are independent and exactly uniform.
#' @export
#'
#' @examples
#' ## Get an approximation of gamma
#' gamma_function = gamma_bootstrapper(10)
#' ## Evaluate it in a number, say .2
#' gamma_function(.2)
gamma_bootstrapper = function(m,
                               n = Inf,
                               B = 1e3,
                               mc.cores = 1L,
                               tau = NULL,
                               K = NULL) {
  if (!is.finite(n)) .GetMinima = function(y, n) which.min(y)

  if (m == 1) {
    return(identity)
  } else if (m < 1) {
    stop("Please supply a positive integer for m")
  } else {
    forCDF = parallel::mclapply(
      1:B,
      function(i) {
        pvals = sort(stats::runif(m))

        if (!is.null(tau) & !is.null(K)) {
          stop("At most one of tau and K can be non NULL")
        } else if (!is.null(tau)) {
          pvals = if (sum(pvals <= tau) > 0) sort(pvals[pvals <= tau]) else min(pvals)
        } else if (!is.null(K)) {
          pvals = sort(pvals)[1:K]
        } else {
          pvals = pvals[order(pvals)]
        }

        if (n >= m) {
          Z = MakeZ_C(pvals, m)
        } else {
          Z = MakeZ_C_nsmall(pvals, n, m)
        }
      },
      mc.cores = mc.cores
    )

    return(
      function(x) mean(forCDF < x, na.rm = T)
    )
  }
}
