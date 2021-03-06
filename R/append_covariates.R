#' @export
append_covariates <- function(model_func,
                              data, zyg, ...,
                              covs1 = character(0), covs2 = character(0)) {
  # Based on Schwabe et al. (2016): https://link.springer.com/article/10.1007/s10519-015-9771-1
  # covs1 - covariates that can mismatch within the pair (expected variables 'X1', 'X2' in the table)
  # covs2 - covariates that matches within the pair (covs1 = 'X' - expected variable 'X' in the table)
  if (!all(c('data', 'zyg') %in% formalArgs(model_func)))
    stop('This model_func is not compatible with this constructor, arguments data and zyg must be present')
  model <- model_func(data = data, zyg = zyg, ...)
  #  model <- univariate_ace(twinData2, 'zyg2', 'bmi')

  nc1 <- length(covs1)
  nc2 <- length(covs2)
  if (nc1 + nc2 == 0) {
    return(model)
  }

  sep <- list(...)$sep
  if (length(sep) == 0)
    sep <- ''

  covNames <- c(paste(covs1, rep(1:2, each = nc1), sep = sep), covs2)
  if (!all(covNames %in% names(data)))
    stop("Some covariates are not in the data table")

  selvars <- c(model$MZ$expectation$dims, covNames)

  mzdata <- data[data[[zyg]] == 'MZ', selvars]
  dzdata <- data[data[[zyg]] == 'DZ', selvars]

  # TODO: Make dummy variables for cathegorical covariates

  nv <- mxEval(ncol(expMeans), model, compute = TRUE) / 2 # number of phenotypic
  nv <- as.vector(nv)

  pheno_names <- unique(
    gsub(
      paste0(sep, '[1|2]$'),
      '',
      model$MZ$expectation$dims)
  )

  phenoCovMZ <- model$expCovMZ
  phenoCovMZ$name <- 'phenoCovMZ'

  phenoCovDZ <- model$expCovDZ
  phenoCovDZ$name <- 'phenoCovDZ'

  phenoMeans <- model$expMeans
  phenoMeans$name <- 'phenoMeans'

  expMZ <- model$MZ$expectation
  expMZ$dims <- selvars
  expMZ$threshnames <- selvars

  expDZ <- model$DZ$expectation
  expDZ$dims <- selvars
  expDZ$threshnames <- selvars


  model <-
    mxModel(model,
            # Covariate means
            mxMatrix(type = 'Full',
                     nrow = 1, ncol = nc1,
                     value = 0, free = TRUE,
                     name = 'covMeans1'),
            mxMatrix(type = 'Full',
                     nrow = 1, ncol = nc2,
                     value = 0, free = TRUE,
                     name = 'covMeans2'),

            # Updating expMeans
            phenoMeans,
            mxAlgebra(cbind(phenoMeans,
                            covMeans1, covMeans1,
                            covMeans2),
                      name = 'expMeans'),

            # betas
            mxMatrix(type = 'Full',
                     nrow = nc1, ncol = nv,
                     values = 0, free = TRUE,
                     name = 'beta1'),
            mxMatrix(type = 'Full',
                     nrow = nc2, ncol = nv,
                     values = 0, free = TRUE,
                     name = 'beta2'),
            mxAlgebra(rbind(beta1, beta2),
                      dimnames = list(c(covs1, covs2), pheno_names),
                      name = 'beta'),

            # Covariate covariances within twin (CovW, CovWprime)
            mxMatrix(type = 'Lower',
                     nrow = nc1, ncol = nc1,
                     values = 0.5, free = TRUE,
                     name = 'CholCovW'),
            mxAlgebra(CholCovW %*% t(CholCovW),
                      name = 'CovW'),
            mxMatrix(type = 'Zero',
                     nrow = nc1, ncol = nc2,
                     name = 'Append1'),
            mxMatrix(type = 'Zero',
                     nrow = nc2, ncol = nc1 + nc2,
                     name = 'Append2'),
            mxAlgebra(rbind(cbind(CovW, Append1), Append2),
                      name = 'CovWprime'),

            mxMatrix(type = 'Lower',
                     nrow = nc1 + nc2, ncol = nc1 + nc2,
                     values = 5, free = TRUE,
                     name = 'CholCovB'),
            mxAlgebra(CholCovB %*% t(CholCovB),
                      name = 'CovB'),

            # Sector 1: covariation across phenotypic variables explained by covariates
            mxAlgebra(rbind(cbind(t(beta) %&% (CovWprime + CovB), t(beta) %&% CovB),
                            cbind(t(beta) %&% CovB, t(beta) %&% (CovWprime + CovB))),
                      name = 'Sect1'),

            # Sector 2: covariation between ph variables and covariates
            nc1 <- mxMatrix(type = 'Full',
                            ncol = nc1, nrow = 1,
                            values = 1:nc1, free = FALSE,
                            name = 'nc1mat'),
            mxAlgebra(rbind(cbind((t(beta) %*% (CovWprime + CovB))[, nc1mat], t(beta) %*% CovB),
                            cbind((t(beta) %*% CovB)[, nc1mat], t(beta) %*% (CovWprime + CovB))),
                      name = 'Sect2'),

            # Sector3: covariation across covariates
            mxAlgebra(rbind(cbind(CovB[nc1mat, nc1mat] + CovW, CovB[nc1mat, ]),
                            cbind(CovB[, nc1mat], CovB + CovWprime)),
                      name = 'Sect3'),

            # Updating expCovMZ
            phenoCovMZ,
            mxAlgebra(rbind(cbind(phenoCovMZ + Sect1, Sect2),
                            cbind(t(Sect2), Sect3)),
                      name = 'expCovMZ',
                      dimnames = list(selvars, selvars)),

            # Updating expCovDZ
            phenoCovDZ,
            mxAlgebra(rbind(cbind(phenoCovDZ + Sect1, Sect2),
                            cbind(t(Sect2), Sect3)),
                      name = 'expCovDZ',
                      dimnames = list(selvars, selvars)),

            # Replace the data
            mxModel(model$MZ,
                    expMZ,
                    mxData(observed = mzdata, type = 'raw')),
            mxModel(model$DZ,
                    expDZ,
                    mxData(observed = dzdata, type = 'raw'))
    )
  return(model)
}
