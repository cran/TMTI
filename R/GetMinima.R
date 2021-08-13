#' Helper function to return the first local minimum that is strictly smaller than the following n entries in a vector.
#'
#' @param Y A vector of transformed p-values
#' @param n The number describing the kind of local minimum to return; n=1
#' returns the first local minimum, n=Inf returns the global minimum.
#'
#' @return The index of the n'th local minimum
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' p <- runif(10)
#' ## Generate the transformations
#' Y <- pbeta(sort(p), 1:10, 10 + 1 - 1:10)
#' ## Get the first local minimum
#' .GetMinima(Y, 1)
#' ## Get the global minimum:
#' .GetMinima(Y, Inf)  # Note, equivalent to which.min(Y)

.GetMinima <- function (Y, n) {
  if (!is.finite(n)) {
    return (
      which.min(Y)
    )
  } else if (n == 1) {
    z <- c(1, Y, 1)

    tf <- (z < dplyr::lag(z)) & (z < dplyr::lead(z))
    tf <- tf[2:(length(Y) + 1)]

    return(which(tf)[1])
  } else  {
    m <- length(Y)

    K <- list()
    K[[1]] <- 1
    for (i in 2:m) {
      vec <- K[[i - 1]] + 1
      checker <- (which (
        Y[K[[i - 1]]] > Y[(vec:m)]
      ) + vec - 1)[1]

      if (identical(checker, numeric(0)) | is.na(checker)) {
        break
      } else {
        K[[i]] <- checker
      }
    }

    K <- do.call("c", K)

    kmax <- length(K)
    K2 <- K[1:(kmax - 1)]
    check1 <- which (
      (K2 + n) < dplyr::lead(K)[-kmax]
    )[1]

    if (identical(check1, integer(0)) | is.na(check1))

      kn <- kmax
    else
      kn <- min(c(kmax, check1))

    return (K[kn])
  }
}
