#' @rdname twin_ref_models
#' @title Define twin reference models
#'
#' @description
#' Define and fit saturated and indepencence models corresponding to a twin
#' model. The twin model has to be fit as the function takes its estimates
#' of covariances and means to compute starting values. The twin model has
#' to include MZ and DZ sub-models.
#'
#'
#' @param model a twin model, MxModel object
#' @param run whether the reference models are to be fit, default is TRUE
#' @param ... the parameters passed to mxRun
#'
#' @return
#'
#' @examples
#'
#' @export
twin_ref_models <- function(model, run = FALSE, ...) {
  selvars <- model$MZ$expectation$dims
  nv <- length(selvars) / 2

  startmz <- 0.9 * attr(model$MZ$fitfunction$result,'expCov')
  startmz[upper.tri(startmz)] <- t(startmz)[upper.tri(startmz)]
  startdz <- 0.9 * attr(model$DZ$fitfunction$result,'expCov')
  startdz[upper.tri(startdz)] <- t(startdz)[upper.tri(startdz)]
  startmean <- 0.9 * attr(model$MZ$fitfunction$result,'expMean')[1:nv]

  labs <- matrix(NA, nrow = nv, ncol = nv)
  labs[lower.tri(labs, diag = TRUE)] <- paste0('V', 1: (nv * (nv + 1) / 2))
  labs[upper.tri(labs)] <- t(labs)[upper.tri(labs)]
  labs <- rbind(cbind(labs, matrix(NA, ncol = nv, nrow = nv)),
                cbind(matrix(NA, ncol = nv, nrow = nv), labs))

  sat_model <- mxModel(name = 'Saturated',
                       mxMatrix(type = 'Symm', ncol = 2 * nv, nrow = 2 * nv,
                                free = TRUE, values = t(startmz), labels = labs,
                                name = 'expCovMZ'),
                       mxMatrix(type = 'Symm', ncol = 2 * nv, nrow = 2 * nv,
                                free = TRUE, values = t(startdz), labels = labs,
                                name = 'expCovDZ'),
                       mxMatrix(type = 'Full', nrow = 1, ncol = nv,
                                free = TRUE, values = startmean, name = 'mean'),
                       mxAlgebra(cbind(mean, mean),
                                 name = 'expMeans'),
                       mxModel(name = 'MZ',
                               model$MZ$data,
                               mxExpectationNormal(covariance = 'Saturated.expCovMZ',
                                                   means = 'Saturated.expMeans',
                                                   dimnames = selvars),
                               mxFitFunctionML()),
                       mxModel(name = 'DZ',
                               model$DZ$data,
                               mxExpectationNormal(covariance = 'Saturated.expCovDZ',
                                                   means = 'Saturated.expMeans',
                                                   dimnames = selvars),
                               mxFitFunctionML()),
                       mxFitFunctionMultigroup(c('MZ','DZ')))

  ind_model <- mxModel(name = 'Independence',
                       mxMatrix(type = 'Diag', ncol = nv, nrow = nv,
                                free = TRUE, values = diag2vec(startmz)[1:nv],
                                name = 'Chol'),
                       mxMatrix(type = 'Zero', ncol = nv, nrow = nv,
                                name = 'Z'),
                       mxAlgebra(Chol %*% t(Chol), name = 'V'),
                       mxAlgebra(rbind(cbind(V, Z),
                                       cbind(Z, V)),
                                 name = 'expCov'),
                       mxMatrix(type = 'Full', nrow = 1, ncol = nv,
                                values = startmean,
                                free = TRUE, name = 'mean'),
                       mxAlgebra(cbind(mean, mean),
                                 name = 'expMeans'),
                       mxModel(name = 'MZ',
                               model$MZ$data,
                               mxExpectationNormal(covariance = 'Independence.expCov',
                                                   means = 'Independence.expMeans',
                                                   dimnames = selvars),
                               mxFitFunctionML()),
                       mxModel(name = 'DZ',
                               model$DZ$data,
                               mxExpectationNormal(covariance = 'Independence.expCov',
                                                   means = 'Independence.expMeans',
                                                   dimnames = selvars),
                               mxFitFunctionML()),
                       mxFitFunctionMultigroup(c('MZ','DZ')))

  if (run) {
    return(list(Saturated = mxRun(sat_model, ...),
                Independence = mxRun(ind_model, ...)))
  } else {
    return(list(Saturated = sat_model,
                Independence = ind_model))
  }
}

adjust_twin_data <- function(DV, IV, interaction = TRUE) {
  result <- DV
  for (i in deps){
    LM <- lm(as.formula(
      paste0(i, '~', paste(indeps,
                           collapse = ifelse(interaction, '*', '+')))),
      data = data)
    outp[, i] <- LM$residuals[rownames(outp)]
  }
}
