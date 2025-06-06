.onAttach <- function(libname, pkgname) {
  if (length(getOption('TwinAnalysis.sep')) == 0)
    options(TwinAnalysis.sep = '')
}

#' @import OpenMx mlth.data.frame glue

#' @rdname twin_ref_models
#' @title Define twin reference models
#'
#' @description
#' Define and fit saturated and indepencence models that follow the assumptions of twin model.
#' The twin model has to be fit as the function takes its estimates
#' of covariances and means to compute starting values. The twin model has
#' to include MZ and DZ sub-models.
#'
#'
#' @param model a twin model, \code{MxModel} object
#' @param run whether the reference models are to be fit, default is TRUE
#' @param ... the parameters passed to \code{ref_models}
#'
#' @return
#'
#' @export
twin_ref_models <- function(
  model,
  run = FALSE,
  ...
) {
  selvars <- model$MZ$expectation$dims
  nv <- length(selvars) / 2

  # FIXME: Doesn't work when the model hasn't been   fit
  startmz <- 0.9 * attr(model$MZ$fitfunction$result,'expCov')
  startmz[upper.tri(startmz)] <- t(startmz)[upper.tri(startmz)]
  startdz <- 0.9 * attr(model$DZ$fitfunction$result,'expCov')
  startdz[upper.tri(startdz)] <- t(startdz)[upper.tri(startdz)]
  startmean <- 0.9 * attr(model$MZ$fitfunction$result,'expMean')[1:nv]

  labs <- matrix(NA, nrow = nv, ncol = nv)
  labs[lower.tri(labs, diag = TRUE)] <- paste0('V', 1: (nv * (nv + 1) / 2))
  labs[upper.tri(labs)] <- t(labs)[upper.tri(labs)]
  labs <- rbind(
    cbind(labs, matrix(NA, ncol = nv, nrow = nv)),
    cbind(matrix(NA, ncol = nv, nrow = nv), labs)
  )

  sat_model <- mxModel(
    name = 'Saturated',
    mxMatrix(
      type = 'Symm',
      ncol = 2 * nv,
      nrow = 2 * nv,
      free = TRUE,
      values = t(startmz),
      labels = labs,
      name = 'expCovMZ'
    ),
    mxMatrix(
      type = 'Symm',
      ncol = 2 * nv,
      nrow = 2 * nv,
      free = TRUE,
      values = t(startdz),
      labels = labs,
      name = 'expCovDZ'
    ),
    mxMatrix(
      type = 'Full',
      nrow = 1,
      ncol = nv,
      free = TRUE,
      values = startmean,
      name = 'mean'
    ),
    mxAlgebra(
      cbind(mean, mean),
      name = 'expMeans'
    ),
    mxModel(
      name = 'MZ',
      model$MZ$data,
      mxExpectationNormal(
        covariance = 'Saturated.expCovMZ',
        means = 'Saturated.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),
    mxModel(
      name = 'DZ',
      model$DZ$data,
      mxExpectationNormal(
        covariance = 'Saturated.expCovDZ',
        means = 'Saturated.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),
    mxFitFunctionMultigroup(
      c('MZ','DZ')
    )
  )

  ind_model <- mxModel(
    name = 'Independence',
    mxMatrix(
      type = 'Diag',
      ncol = nv,
      nrow = nv,
      free = TRUE,
      values = diag2vec(startmz)[1:nv],
      name = 'Chol'
    ),
    mxMatrix(
      type = 'Zero',
      ncol = nv,
      nrow = nv,
      name = 'Z'
    ),
    mxAlgebra(
      Chol %*% t(Chol),
      name = 'V'
    ),
    mxAlgebra(
      rbind(
        cbind(V, Z),
        cbind(Z, V)
      ),
      name = 'expCov'
    ),
    mxMatrix(
      type = 'Full',
      nrow = 1,
      ncol = nv,
      values = startmean,
      free = TRUE,
      name = 'mean'
    ),
    mxAlgebra(
      cbind(mean, mean),
      name = 'expMeans'
    ),
    mxModel(
      name = 'MZ',
      model$MZ$data,
      mxExpectationNormal(
        covariance = 'Independence.expCov',
        means = 'Independence.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),
    mxModel(
      name = 'DZ',
      model$DZ$data,
      mxExpectationNormal(
        covariance = 'Independence.expCov',
        means = 'Independence.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),
    mxFitFunctionMultigroup(
      c('MZ','DZ')
    )
  )

  model <- ref_models(
    model,
    ref_models = list(
      Saturated = sat_model,
      Independence = ind_model
    ),
    run = run,
    ...
  )
  return(model)
}

#' @title Process twin data into MxData objects
#'
#' @description Process twin data into MxData objects to add to MxModel.
#' The function accepts raw data and covariance matrices.
#'
#' @param data \code{data.frame} or \code{list}
#' @param data_type 'raw', 'cov' or 'cor', as required by \code{type} parameter
#' of \code{mxData}
#' @param zyg factor used to split the data frame when data_type = 'raw'
#'
#' @details
#' When data_type = 'raw' (default), the data must be a data frame or matrix
#' with raw data that will be split upon factor \code{zyg} and passed to
#' \code{mxData}. When data_type = 'cov' or 'cor', the data  must be a (named)
#' list of arguments passed to \code{mxData}, including \code{observed} -
#' covariation or correlation matrix, \code{numObs} - number of observations,
#' \code{means} - vector of means. The \code{type} argument is not required as
#' it is added from \code{data_type}. Refer to \code{?mxData} for details on
#' \code{mxData} arguments.
#'
#' @return A list of \code{MxData} objects
#'
#' @export
process_twin_data <- function(data,
                              selvars,
                              data_type = 'raw',
                              zyg = character(0)) {

  if (data_type == 'raw') {
    # Checking that zygosity is indeed a factor
    zygosity <- data[, zyg]
    if (!is.factor(zygosity)) {
      if (all(zygosity %in% 1:2)) {
        zygosity <- factor(
          zygosity,
          levels = 1:2,
          labels = c('MZ', 'DZ')
        )
      } else if (all(zygosity %in% c('MZ', 'DZ'))) {
        zygosity <- factor(
          zygosity,
          levels = c('MZ', 'DZ')
        )
      } else {
        stop('zygosity must be labeled as either "MZ" and "DZ" or 1 and 2')
      }
    } else if (!all(levels(zygosity) %in% c('MZ', 'DZ'))) {
      stop('zygosity levels must be "MZ" and "DZ"')
    }

    data <- split(
      as.data.frame(data[, selvars]),
      zygosity
    )

    lapply(data, function(x) {
      mxData(observed = x,
             type = 'raw')
    })
  } else {
    data <- lapply(data, function(x) {
      x$observed <- as.matrix(x$observed[selvars, selvars])
      x
    })
    data <- lapply(data, function(x) c(x, list(type = data_type)))
    means <- lapply(data, `[[`, 'means')
    if (any(sapply(means, length) == 0) || any(sapply(means, is.na))) {
      data <- lapply(data,
                     function(x)
                       c(x, list(means = setNames(rep(0, length(selvars)),
                                                  selvars))))
    }
    lapply(data, do.call, what = mxData)
  }
}

#' @export
compute_starting_values <- function(model, data,
                                    data_type = 'raw',
                                    vars = character(0),
                                    sep = '',
                                    ...) {
  # This function computes starting values for a given model
  # It assumes that the variables are coded var1 and var2 for two twins
  if (data_type == 'raw') {
    longdata <- reshape(data,
                        varying = lapply(vars, paste, 1:2, sep = sep),
                        v.names = vars,
                        direction = 'long',
                        sep = sep)
    startmean <- 0.9 * colMeans(longdata[, vars, drop = FALSE], na.rm = TRUE)
    freemean <- TRUE
    covmat <- cov(longdata[, vars, drop = FALSE], use = 'pairwise')
  } else {
    means <- lapply(data, `[[`, 'means')
    if (all(sapply(means, length) > 0)) {
      startmean <- 0.9 * Reduce(`+`, means) / length(means)
      freemean <- TRUE
    } else {
      startmean <- 0
      freemean <- FALSE
    }

    covs <- lapply(data, `[[`, 'observed')
    covs <- lapply(covs, function(x) {
      (x[paste(vars, 1, sep = sep), paste(vars, 1, sep = sep), drop = FALSE] +
         x[paste(vars, 2, sep = sep), paste(vars, 2, sep = sep), drop = FALSE]) / 2
    })

    covmat <- Reduce(`+`, covs) / length(covs)
    covmat <- as.matrix(covmat)
    dimnames(covmat) <- list(vars, vars)
  }

  outp <- list(means = startmean,
               freemean = freemean)

  dots <- list(...)

  if (model == 'univariate_ace') {
    outp$var <- 0.9 * covmat
  }

  if (model == 'multivariate_ace') {
    outp$chol <- t(chol(0.9 * covmat))
  }

  if (model == 'cross_lag_ace') {
    pheno <- cross_lag(list(observed = covmat / 3,
                            numObs = 1000),
                       data_type = 'cov',
                       definition = dots$definition)
    pheno <- mxRun(pheno)
    outp$pheno <- pheno
  }

  return(outp)
}

adjust_twin_data <- function(DV, IV, interaction = TRUE) {
  # Not used
  result <- DV
  for (i in deps){
    LM <- lm(as.formula(
      paste0(i, '~', paste(indeps,
                           collapse = ifelse(interaction, '*', '+')))),
      data = data)
    outp[, i] <- LM$residuals[rownames(outp)]
  }
}

