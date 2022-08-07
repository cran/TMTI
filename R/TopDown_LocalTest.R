#' TopDown LocalTest algorithm for estimating a 1-alpha confidence set for the number
#' of false hypotheses among a set.
#'
#' @name TopDown_LocalTest
#' @aliases TopDown_localTest
#' @param LocalTest A function specifying a local test.
#' @param pvals A vector of p-values.
#' @param subset Numeric vector specifying a subset a p-values to estimate a
#' confidence set for the number of false hypotheses for. Defaults to NULL
#' corresponding to estimating a confidence set for the number of false
#' hypotheses in the entire set.
#' @param alpha Level in [0,1] at which to generate confidence set. Defaults
#' to 0.05.
#' @param verbose Logical, indicating whether or not to write out the progress.
#' Defaults to TRUE.
#' @param mc.cores Integer specifying the number of cores to parallelize onto.
#' @param chunksize Integer indicating the size of chunks to parallelize. E.g.,
#' if setting chunksize = mc.cores, each time a parallel computation is set up,
#' each worker will perform only a single task. If mc.cores > chunksize, some
#' threads will be inactive.
#' @param direction A string indicating whether to perform a binary search ('binary'/'b')
#' or decreasing ('decreasing'/'d') search. Defaults to 'binary', which has better
#' computational complexity.
#' @param ... Additional parameters.
#'
#' @return A 1-alpha bound lower for the number of false hypotheses among the
#' specified subset of the supplied p-values
#' @export
#'
#' @examples
#' ## Simulate some p-values
#' ## The first 10 are from false hypotheses, the next 10 are from true
#' pvals = c(
#'   rbeta(10, 1, 20), ## Mean value of .05
#'   runif(10)
#' )
#' ## Estimate the confidence set using a local Bonferroni test
#' TopDown_LocalTest(function(x) {
#'   min(c(1, length(x) * min(x)))
#' }, pvals)
#'
TopDown_LocalTest = function(LocalTest,
                             pvals,
                             subset = NULL,
                             alpha = 0.05,
                             verbose = FALSE,
                             mc.cores = 1L,
                             chunksize = 4 * mc.cores,
                             direction = 'binary',
                              ...) {
  if (is.unsorted(subset))
    subset = sort(subset)
  ord = order(pvals)
  pvals = sort(pvals)
  m = length(pvals)
  t_alpha = 0

  if (!is.null(subset) & length(subset) < length(pvals)) {
    if (mc.cores <= 1L) {
      if (any(tolower(direction) == c('d', 'decreasing'))) {
        counter = 0
        top = length(subset)
        for (i in seq_along(subset)) {
          counter = counter + 1
          if (verbose) {
            cat(
              sprintf(
                "\rOuter step %i of %i\n", counter, length(subset)
              )
            )
          }
          subset2 = subset[i:length(subset)]
          p_loc = TestSet_LocalTest(
            LocalTest,
            pvals,
            subset2,
            alpha = alpha,
            EarlyStop = TRUE,
            verbose = verbose,
            mc.cores = mc.cores,
            ...
          )
          accept = (p_loc >= alpha)
          if (accept) {
            t_alpha = length(subset2)
            break
          }
        }
        out = top - t_alpha
      } else if (any(tolower(direction) == c('b', 'binary'))) {
        return(TopDown_C_binary_subset (
          LocalTest = LocalTest,
          pSub = pvals[subset],
          pRest = pvals[-subset],
          alpha = alpha,
          low = 0,
          high = length(subset) - 1,
          verbose = verbose
        ))
      }
    } else {
      if (any(tolower(direction) == c('d', 'decreasing'))) {
        l_sub = length(subset)
        chunks = split(seq_along(subset), ceiling(seq_along(subset) / chunksize))
        results = list()
        counter = 1
        for (x in chunks) {
          results[[counter]] = unlist(parallel::mclapply (
            x,
            function (i) {
              TestSet_C (
                LocalTest = LocalTest,
                pSub = pvals[subset[i:l_sub]],
                pRest = pvals[-subset[i:l_sub]],
                alpha = alpha,
                is_subset_sequence = FALSE,
                EarlyStop = TRUE,
                verbose = FALSE
              )
            },
            mc.cores = mc.cores
          ))
          if (any(unlist(results) > alpha)) {
            break
          }
          counter = counter + 1
        }
        t_alpha = length(subset) + 1 - which(unlist(results) > alpha)[1]
        out = length(subset) - t_alpha
      } else if (any(tolower(direction) == c('b', 'binary'))) {
        low  = 1
        high = length(subset)
        while (TRUE) {
          mid = floor((low + high) / 2)
          p = TMTI::TestSet_LocalTest (
            LocalTest = LocalTest,
            pvals = pvals,
            subset = subset[mid:length(subset)],
            alpha = alpha,
            EarlyStop = TRUE,
            verbose = FALSE,
            mc.cores = mc.cores,
            is.sorted = TRUE
          )
          if (p < alpha) {
            low = mid + 1
          } else {
            high = mid
          }
          if (verbose)
            cat(sprintf("\nlow = %i, mid = %i, high = %i, p = %.2f%%          ",
                        low, mid, high, p * 100))
          if (low == length(subset) & mid < low) {
            pnew = TMTI::TestSet_LocalTest (
              LocalTest = LocalTest,
              pvals = pvals,
              subset = subset[length(subset)],
              alpha = alpha,
              EarlyStop = TRUE,
              verbose = FALSE,
              mc.cores = mc.cores,
              is.sorted = TRUE
            )
            if (pnew < alpha)
              R = length(subset) + 1
          } else if (high == low) {
            if (p < alpha)
              R = mid + 1
            else
              R = mid
            break
          } else if (low > high) {
            R = low
            break
          }
        }
        return (R - 1)
      }
    }
  } else {
    if (mc.cores <= 1L) {
      if (any(tolower(direction) == c('d', 'decreasing')))  {
        out = TopDown_C (
          LocalTest = LocalTest,
          pvals = pvals,
          alpha = alpha
        )
      } else if (any(tolower(direction) == c('b', 'binary')))  {
        return (
          TopDown_C_binary (
            LocalTest = LocalTest,
            pvals = pvals,
            alpha = alpha,
            low = 0,
            high = length(pvals) - 1,
            verbose = verbose
          )
        )
      }
    } else {
      if (any(tolower(direction) == c('b', 'binary')))
        message(
          "Binary search not supported for mc.cores > 1 when computing CSs among all hypotheses. Reverting to direction = decreasing"
        )
      .f = function(i) {
        LocalTest(pvals[i:m])
      }
      chunks = split(seq(m), ceiling(seq(m) / chunksize))
      results = list()
      counter = 1
      for (x in chunks) {
        if (verbose) {
          cat(sprintf("\rProcessing chunk %i of %i", counter, length(chunks)))
        }
        results_ = parallel::mclapply(
          x,
          .f,
          mc.cores = mc.cores
        )
        results[[counter]] = unlist(results_)
        if (any(unlist(results) > alpha)) {
          break
        }
        counter = counter + 1
      }
      t_alpha = m + 1 - which(unlist(results) > alpha)[1]
      out = length(pvals) - t_alpha
    }
  }

  if (verbose) {
    cat(
      paste0(
        "Confidence set for the number of false hypotheses is {",
        out,
        ",..., ",
        if (is.null(subset)) length(pvals) else length(subset),
        "}\n"
      )
    )
  }
  return(out)
}


#'
#' @rdname TopDown_LocalTest
#' @param localTest A function specifying a local test (deprecated).
#' @export
TopDown_localTest = function(localTest,
                             pvals,
                             subset = NULL,
                             alpha = 0.05,
                             verbose = TRUE,
                             mc.cores = 1L,
                             chunksize = 4 * mc.cores,
                             ...) {
  .Deprecated(new = "TopDown_LocalTest")
  TopDown_LocalTest (
    LocalTest = localTest,
    pvals = pvals,
    subset = NULL,
    alpha = 0.05,
    verbose = TRUE,
    mc.cores = 1L,
    chunksize = 4 * mc.cores,
    ...
  )
}
