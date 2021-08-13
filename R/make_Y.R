#' Returns the transformed p-values to be used in the TMTI, tTMTI or rtTMTI. Internal function.
#'
#' @param pvals A vector of p-values
#' @param tau A numeric between 0 and 1 indicating the truncation level.
#' Defaults to NULL
#' @param K An integer > 1 indicating the the truncation index. Defaults
#' to NULL
#' @param log.p Logical, indicating whether transformations are on log scale.
#' Defaults to TRUE.
#'
#' @return A vector, Y, of transformed p-values
#' @export
#'
#' @examples
#' ## Simulate p-values
#' p <- runif(10)
#' make_Y(p)  # Normal Y transformation
#' make_Y(p, tau = 0.5) # Using only the p-values below .5
#' make_Y(p, K = 5) # Using only the five smallest p-values.
make_Y <- function (
  pvals,
  tau    = NULL,
  K      = NULL,
  log.p  = TRUE
) {
  m <- length(pvals)

  if (!is.null(tau) & !is.null(K))
    stop("Only one of the arguments tau and K can be different from NULL")
  else if (!is.null(tau)) {
    p_subset <- 1:sum(pvals < tau)
    if (identical(p_subset, numeric(0))) {
      warning("There are no p-values below tau, using the smallest p-value instead")
      p_subset <- 1
    }
  } else if (!is.null(K))
    p_subset <- 1:K
  else
    p_subset <- 1:m

  stats::pbeta (
    sort(pvals),
    1:m,
    m + 1 - 1:m,
    log.p = log.p
  )[p_subset]

}
