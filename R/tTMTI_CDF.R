#' Computes the analytical version of the tTMTI_infty CDF. When m>100, this should
#' not be used.
#'
#' @param x Point in which to evaluate the CDF
#' @param m Number of independent tests to combine
#' @param tau The truncation point of the tTMTI procedure
#'
#' @return The probability that the test statistic is at most x assuming
#' independence under the global null hypothesis.
#' @export
#'
#' @examples
#' tTMTI_CDF(0.05, 100, 0.05)

tTMTI_CDF <- function (x, m, tau) {
  if (m == 1) {
    x
  }

  P <- function(x, a) { ## Constructs the necessary polynomials
    m  <- length(a)

    sum(1 / factorial(m:1) * x^(m:1) * a)
  }

  xs <- numeric(m)
  xs <- stats::qbeta(x, 1:m, m + 1 - 1:m)

  PP <- list()
  PP[[1]] <- xs[1]
  for (i in 2:m) {
    PP[[i]] <- P(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  hh_terms <- c (
    1 - (1 - PP[[1]] / tau) * (xs[1] <= tau),
    # ifelse(xs[1] <= tau, PP[[1]] / tau, 1),
    sapply (
      2:m,
      function (K) {
        Ptau <- P(tau, c(1, -do.call("c", PP[1:(K - 1)])))

        1 - (factorial(K) / tau**K) * (
          Ptau -
            PP[[K]]
        ) * (xs[K] <= tau)
      }
    )
  )
  tau_terms <- sapply (
    1:m,
    function (i) {
      choose(m, i) * tau^i * (1 - tau)^(m - i)
    }
  )

  pbeta_tau <- stats::pbeta(tau, 1, m)

  (1 - tau)**m * (x - pbeta_tau) / (1 - pbeta_tau) * (xs[1] > tau) + sum(hh_terms * tau_terms)
}
