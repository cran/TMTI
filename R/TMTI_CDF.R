#' Computes the analytical version of the TMTI_infty CDF. When m>100, this should
#' not be used.
#'
#' @param x Point in which to evaluate the CDF.
#' @param m Number of independent tests to combine.
#'
#' @return The probability that the test statistic is at most x assuming
#' independence under the global null hypothesis.
#' @export
#'
#' @examples
#' TMTI_CDF(0.05, 100)
TMTI_CDF = function(x, m) { ## This is the explicit form of gamma
  P = function(x, a) { ## Constructs the necessary polynomials
    m = length(a)

    sum(1 / factorial(m:1) * x^(m:1) * a)
  }

  xs = stats::qbeta(x, 1:m, m + 1 - 1:m)

  PP = list()
  PP[1] = xs[1]
  for (i in 2:m) {
    PP[[i]] = P(xs[i], c(1, -do.call("c", PP[1:(i - 1)])))
  }

  1 - factorial(m) * (
    P(1, c(1, -do.call("c", PP[1:(m - 1)]))) -
      PP[[m]]
  )
}
