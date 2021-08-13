#' Function to bootstrap the Cumulative Distribution Functions (CDFs) of the TMTI statistics.
#'
#' @param m Number of tests
#' @param n Number (or Inf) indicating what kind of minimum to consider.
#' Defaults to Inf, corresponding to the global minimum.
#' @param B Number of bootstrap replicates. Rule of thumb is to use at least
#' 10 * m
#' @param log.p Logical indicating whether to calculate p-values on log-scale.
#' Defaults to TRUE
#' @param mc.cores Integer denoting the number of cores to use when using
#' parallelization, Defaults to 1, corresponding to single-threaded computations
#' @param tau Numerical (in (0,1)); threshold to use in tTMTI. If set to NULL,
#' then either TMTI (default) or rtTMTI is used.
#' @param K Integer; Number of smallest p-values to use in rtTMTI. If se to NULL,
#' then either TMTI (default) or tTMTI is used.
#'
#' @return An approximation of the function \eqn{\gamma^m(x)} under the
#' assumption that all p-values are independent and exactly uniform
#' @export
#'
#' @examples
#' ## Get an approximation of gamma
#' gamma_function <- gamma_bootstrapper(10)
#' ## Evaluate it in a number, say .2
#' gamma_function(.2)
gamma_bootstrapper <- function (
  m,
  n = Inf,
  B = 1e3,
  log.p = TRUE,
  mc.cores = 1L,
  tau = NULL,
  K = NULL
) {
  if (!is.finite(n)) .GetMinima <- function(y, n) which.min(y)

  if (m == 1)
    return(identity)
  else if (m < 1)
    stop("Please supply a positive integer for m")
  else {
    # if (log.p)
    #   xtrans <- log
    # else
    #   xtrans <- identity

    forCDF <- parallel::mclapply (
      1:B,
      function (i) {
        pvals <- sort(stats::runif(m))
        # Y <- pbeta(pvals, 1:m, m + 1 - 1:m, log.p = log.p)
        Y <- make_Y(pvals, tau = tau, K = K, log.p = log.p)

        Y[.GetMinima(Y, n)]
      },
      mc.cores = mc.cores
    )

    return (
      function (x) mean(forCDF < x, na.rm = T)
    )
  }
}
