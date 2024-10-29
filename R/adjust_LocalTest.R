#' Adjust all p-values using a Closed Testing Procedure and a
#' user-defined local test which satisfies the quadratic shortcut given in Mogensen and Markussen (2021)
#'
#' @param LocalTest A function specifying a local test.
#' @param pvals vector of p-values.
#' @param alpha significance level. Defaults to 0.05.
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
#' adjust_LocalTest(
#'   LocalTest = function(x) {
#'     min(c(1, length(x) * min(x)))
#'   },
#'   p, alpha = 0.05, is.sorted = TRUE
#' )
adjust_LocalTest = function(LocalTest,
                             pvals, alpha = 0.05,
                             is.sorted = FALSE,
                             EarlyStop = FALSE,
                             verbose = FALSE,
                             mc.cores = 1L,
                             chunksize = 4 * mc.cores,
                             direction = "increasing",
                             parallel.direction = "breadth",
                             AdjustAll = FALSE,
                             ...) {
  if (AdjustAll & mc.cores <= 1) {
    return (
      CTP_LocalTest(
        LocalTest = LocalTest,
        pvals = pvals,
        alpha = alpha,
        is.sorted = !is.unsorted(pvals),
        EarlyStop = EarlyStop,
        ...
      )
    )
  }

  if (AdjustAll) {
    m2 = length(pvals)
  } else {
    m2 = sum(pvals <= alpha)

    if (m2 <= 0) {
      stop("There are no p-values that are marginally significant at level alpha")
    }
    if (verbose) {
      cat(sprintf("\rThere are %i marginally significant p-values to adjust", m2))
    }
  }

  if (is.sorted) {
    ord = seq_along(pvals)
  } else {
    ord = order(pvals)
  }

  if (mc.cores == 1) {
    .f = function(i, es = EarlyStop) {
      if (verbose) {
        cat(sprintf("\nAdjusting p-value %i of %i.\n", i, m2))
      }

      out = TMTI::TestSet_LocalTest(
        LocalTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        EarlyStop = es,
        verbose = verbose,
        is.sorted = is.sorted,
        ...
      )
      out
    }
    results = list()
    if (tolower(direction) == "increasing" | tolower(direction) == "i") {
      for (i in seq(m2)) {
        results[[i]] = .f(i)
        if (results[[i]] >= alpha & EarlyStop) {
          # message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
          break
        }
      }
      return(
        data.frame(
          "p"     = unlist(results),
          "index" = ord[1:length(results)]
        )
      )
    } else if (tolower(direction) == "decreasing" | tolower(direction) == "d") {
      for (i in m2:1) {
        results[[i]] = .f(i, TRUE)
        if ((results[[i]] < alpha) & (i > 1)) {
          # message(paste0("Adjusted p-value ", i, " is below alpha, implying that the remaining are also below alpha. Exiting"))
          brokeEarly = TRUE
          break
        }
      }
      if (!exists("brokeEarly") | i == 1) {
        return(
          data.frame(
            "p"     = unlist(results),
            "index" = ord[length(results):1]
          )
        )
      } else {
        nonnull = length(unlist(results))

        return(
          data.frame(
            "p" = c(
              rep(paste0("p > ", alpha), nonnull - 1),
              results[[i]],
              rep(paste0("p <= ", results[[i]]), m2 - nonnull)
            ),
            "index" = ord[m2:1]
          )
        )
      }
    } else if (tolower(direction) == "binary" | tolower(direction) == "b") {
      return(FWER_set_C(
        LocalTest = LocalTest,
        pvals = pvals,
        alpha = alpha,
        low = 0,
        high = length(pvals) - 1,
        verbose = verbose
      ))
    }
  } else if (any(tolower(parallel.direction) == c("d", "depth"))) {
    .f = function(i, es = EarlyStop) {
      if (verbose) {
        cat(sprintf("\rAdjusting p-value %i of %i.", i, m2))
      }

      out = TMTI::TestSet_LocalTest(
        LocalTest,
        pvals = pvals,
        subset = ord[i],
        alpha = alpha,
        EarlyStop = es,
        verbose = verbose,
        is.sorted = is.sorted,
        mc.cores = mc.cores,
        chunksize = chunksize,
        ...
      )
      out
    }
    if (any(tolower(direction) == c("i", "increasing"))) {
      results = list()
      for (i in seq(m2)) {
        results[[i]] = .f(i)
        if (results[[i]] >= alpha & EarlyStop) {
          # message(paste0("Adjusted p-value ", i, " was above alpha, implying that the remaining are also above alpha. Exiting"))
          break
        }
      }
      return(
        data.frame(
          "p"     = unlist(results),
          "index" = ord[1:length(results)]
        )
      )
    } else if (any(tolower(direction) == c('d', 'decreasing'))) {
      results = list()
      for (i in rev(seq(m2))) {
        results[[i]] = .f(i)
        if (results[[i]] < alpha) {
          out = i
          if(verbose) {
            cat(sprintf("There are %i p-values that can be rejected with FWER control.\nAn upper bound for these is %.4f%%",
                      out, 100*results[[i]] ))
          }
          return(i)
        }
      }
      return(
        data.frame(
          "p"     = unlist(results),
          "index" = ord[1:length(results)]
        )
      )
    } else if (any(tolower(direction) == c('b', 'binary'))) {
      low  = 1
      high = if(AdjustAll) m2 else m2 + 1
      while (TRUE) {
        mid = floor((low + high) / 2)
        p = .f(mid, TRUE)
        if (p < alpha) {
          low = mid + 1
        } else {
          high = mid
        }
        if (verbose)
          cat(sprintf("\nlow = %i, mid = %i, high = %i", low, mid, high))
        if (high == low) {
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
      if (verbose) {
        cat(sprintf("\rThere are %i hypotheses that can be rejected with FWER control", R - 1))
      }
      return (R - 1)
    }

  } else if (any(tolower(parallel.direction) == c("b", "breadth"))) {
    if (any(tolower(direction) == c('b', 'binary')))
      stop("Binary search only works with parallel.direction == 'depth'")

    chunks  = if (any(tolower(direction) == c("i", "increasing")))
      split(seq(m2), ceiling(seq(m2) / chunksize))
    else
      split(rev(seq(m2)), ceiling(seq(m2) / chunksize))
    results = list()
    counter = 1
    for (x in chunks) {
      results[[counter]] = unlist(parallel::mclapply (
        x,
        function (j) {
          TMTI::TestSet_LocalTest (
            LocalTest = LocalTest,
            pvals = pvals,
            subset = ord[j],
            alpha = alpha,
            EarlyStop = EarlyStop,
            verbose = verbose,
            mc.cores = 1,
            is.sorted = is.sorted
          )
        },
        mc.cores = mc.cores
      ))
      counter = counter + 1
      if (any(tolower(direction) == c("i", "increasing"))) {
        if ((any(unlist(results) > alpha)) & EarlyStop)
          break
      } else {
        if (any(unlist(results) < alpha)) {
          out = m2 + 1 - which(unlist(results) < alpha)[1]
          if (verbose) {
            cat(sprintf("There are %i p-values that can be rejected with FWER control.\nAn upper bound for these is %.4f%%",
                      out, 100*unlist(results)[which(unlist(results) < alpha)[1]] ))
          }
          return (out)
        }
      }

    }
    return(
      data.frame(
        "p"     = unlist(results),
        "index" = ord[1:length(unlist(results))]
      )
    )
  }
}
