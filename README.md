# TMTI
This repository contains the R package 'TMTI', which implements the methods of Mogensen and Markussen (2021). A preprint of the paper can be found here: https://arxiv.org/abs/2108.04731

# Installation
The package can be installed using e.g. devtools or the remotes package:
```
devtools::install_github("phillipmogensen/TMTI")
# or
remotes::install_github("phillipmogensen/TMTI")
```

# Usage
## Testing global hypotheses
The workhorse function for global hypothesis testing is called `TMTI`. It takes as it's main input a vector of $p$-values, which (unless a $\gamma$-function is explicitly specified) are assumed to be independent. The $p$-value calculation is exact when the length of the input vector is below $100$, but a bootstrap approximation is used for longer vectors.

Additional arguments that can be supplied are (among others) `K` (an integer which indicates a rank truncation index) and `tau` (a float in $(0,1)$ which indicates a truncation threshold).

Example:
```
pvalues <- runif(50)  # 50 independent p-values

TMTI::TMTI(pvalues)
TMTI::TMTI(pvalues, K = 5)  # Uses only the smallest 5 p-values to test the global null
TMTI::TMTI(pvalues, alpha = 0.1)  # Uses only p-values that are marginally significant at level 0.1
```

## Multiple testing
The package includes two functions for conduction Closed Testing Procedures (CTPs) in quadratic time. One function, `TMTI_CTP`, which uses the `TMTI` as local test in every layer of the CTP. The other function, `localTest_CTP`, employs a user-input local test which must satisfy the monotonicity conditions of Mogensen and Markussen (2021). Both functions return the CTP adjusted $p$-values for all elementary hypotheses and control the FWER rate strongly.

Example:
```
pvalues <- c(runif(10, 0, 0.01), runif(10))  # p-values from 10 false hypotheses and 10 true

TMTI::TMTI_CTP (  # CTP using the tTMTI
  pvalues, 
  tau = 0.05
)

TMTI::localTest_CTP (  # CTP using the Fisher Combination Test
  localTest = function (x) {
    if (length(x) == 1L)
      return(x)
    else {
      T_FCT <- -2 * sum(log(x))
      
      return(pchisq(T_FCT, df = 2 * length(x), lower.tail = F))
    }
  },
  pvalues
)
```
