#' @export
univariate_ace <- function(data, zyg, vars) { # vars is length 1
  require(OpenMx)
  if (length(vars) > 1) {
    warning('vars has to be length 1')
    vars <- vars[1]
  }
  selvars <- paste0(vars, 1:2)
  nv <- 1

  # Starting values
  longdata <-reshape(data[, selvars],
                     varying = selvars,
                     v.names = vars,
                     direction = 'long')
  startmean <- 0.9 * mean(longdata[, vars], na.rm = TRUE)
  startvar <- sqrt(0.3 * var(longdata[, vars], na.rm = TRUE))

  mzdata <- data[data[[zyg]] == 'MZ', selvars]
  dzdata <- data[data[[zyg]] == 'DZ', selvars]

  twmodel <- mxModel(name = 'ACE',
                     mxMatrix(type = 'Full',
                              nrow = 1, ncol = nv,
                              free = TRUE, values = startmean,
                              name = 'mean'),
                     mxAlgebra(cbind(mean, mean),
                               name = 'expMeans'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = startvar,
                              name = 'a'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = startvar,
                              name = 'c'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = startvar,
                              name = 'e'),
                     mxAlgebra(t(a) %*% a, name = 'A'),
                     mxAlgebra(t(c) %*% c, name = 'C'),
                     mxAlgebra(t(e) %*% e, name = 'E'),
                     mxAlgebra(rbind(cbind(A + C + E, A + C    ),
                                     cbind(A + C,     A + C + E)),
                               name = 'expCovMZ'),
                     mxAlgebra(rbind(cbind(A + C + E, .5%x%A + C),
                                     cbind(.5%x%A + C, A + C + E)),
                               name = 'expCovDZ'),

                     mxModel(name = 'MZ',
                             mxData(observed = mzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),
                     mxModel(name = 'DZ',
                             mxData(observed = dzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),

                     mxFitFunctionMultigroup(c('MZ','DZ')),

                     # Output tables
                     mxAlgebra(A + C + E, name = 'V'),
                     mxAlgebra(A/V, name = 'Aprop'),
                     mxAlgebra(C/V, name = 'Cprop'),
                     mxAlgebra(E/V, name = 'Eprop'),
                     mxAlgebra(cbind(diag2vec(Aprop),
                                     diag2vec(Cprop),
                                     diag2vec(Eprop)),
                               dimnames = list(vars, c('A', 'C', 'E')),
                               name = 'Variance_components')
                     )
  output_tables(twmodel) <- c('Variance_components')
  twmodel <- def_ci(twmodel, output_tables(twmodel))
  return(twmodel)
}

#' @export
multivariate_ace <- function(data, zyg, vars) {
  require(OpenMx)

  nv <- length(vars)
  selvars <- paste0(vars, rep(1:2, each = nv))

  indnames <- outer(vars, vars,
                    function(x, y) paste0(y, '<->', x))
  indnames <- indnames[row(indnames) > col(indnames)]

  # Starting values
  longdata <-reshape(data[, selvars],
                     varying = selvars,
                     v.names = vars,
                     direction = 'long')
  startmean <- 0.9 * sapply(longdata[, vars], mean, na.rm = TRUE)
  startchol <- chol(0.3 * var(longdata[, vars], na.rm = TRUE))

  mzdata <- data[data[[zyg]] == 'MZ', selvars]
  dzdata <- data[data[[zyg]] == 'DZ', selvars]

  twmodel <- mxModel(name = 'ACE',
                     mxMatrix(type = 'Full',
                              nrow = 1, ncol = nv,
                              free = TRUE, values = startmean,
                              name = 'mean'),
                     mxAlgebra(cbind(mean, mean),
                               name = 'expMeans'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = t(startchol),
                              name = 'a'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = t(startchol),
                              name = 'c'),
                     mxMatrix(type = 'Lower',
                              nrow = nv, ncol = nv,
                              free = TRUE, values = t(startchol),
                              name = 'e'),
                     mxAlgebra(t(a) %*% a, name = 'A'),
                     mxAlgebra(t(c) %*% c, name = 'C'),
                     mxAlgebra(t(e) %*% e, name = 'E'),

                     mxAlgebra(rbind(cbind(A + C + E, A + C    ),
                                     cbind(A + C,     A + C + E)),
                               name = 'expCovMZ'),
                     mxAlgebra(rbind(cbind(A + C + E, .5%x%A + C),
                                     cbind(.5%x%A + C, A + C + E)),
                               name = 'expCovDZ'),

                     mxModel(name = 'MZ',
                             mxData(observed = mzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),
                     mxModel(name = 'DZ',
                             mxData(observed = dzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),
                     mxFitFunctionMultigroup(c('MZ','DZ')),

                     # Output tables
                     mxAlgebra(A + C + E, name = 'V'),
                     mxMatrix(type = "Iden",
                              nrow = nv, ncol = nv,
                              name = "I"),
                     mxAlgebra(solve(sqrt(I * V)), name = "iSD"),
                     mxAlgebra(solve(sqrt(I * A)), name = 'aSD'),
                     mxAlgebra(solve(sqrt(I * C)), name = 'cSD'),
                     mxAlgebra(solve(sqrt(I * E)), name = 'eSD'),

                     mxAlgebra(iSD %&% V, name='Cor'),
                     mxAlgebra(Cor,
                               dimnames = list(vars, vars),
                               name='Total_correlations'),

                     mxAlgebra(A / V, name = 'Aprop'),
                     mxAlgebra(C / V, name = 'Cprop'),
                     mxAlgebra(E / V, name = 'Eprop'),
                     mxAlgebra(cbind(diag2vec(Aprop),
                                     diag2vec(Cprop),
                                     diag2vec(Eprop)),
                               dimnames=list(vars, c('A', 'C', 'E')),
                               name='Variance_components'),

                     mxAlgebra(aSD %&% A,
                               dimnames = list(vars, vars),
                               name = 'Cor_A'),
                     mxAlgebra(cSD %&% C,
                               dimnames = list(vars, vars),
                               name = 'Cor_C'),
                     mxAlgebra(eSD %&% E,
                               dimnames = list(vars, vars),
                               name = 'Cor_E'),
                     mxAlgebra(Cor_A,
                               dimnames = list(vars, vars),
                               name='A_correlations'),
                     mxAlgebra(Cor_C,
                               dimnames = list(vars, vars),
                               name='C_correlations'),
                     mxAlgebra(Cor_E,
                               dimnames = list(vars, vars),
                               name='E_correlations'),

                     mxAlgebra(sqrt(I * Aprop) %&% Cor_A, name = 'Cov_A'),
                     mxAlgebra(sqrt(I * Cprop) %&% Cor_C, name = 'Cov_C'),
                     mxAlgebra(sqrt(I * Eprop) %&% Cor_E, name = 'Cov_E'),
                     mxAlgebra(abs(Cov_A) + abs(Cov_C) + abs(Cov_E), name='absV'),
                     mxAlgebra(abs(Cov_A) / absV, name = 'CorPr_A'),
                     mxAlgebra(abs(Cov_C) / absV, name = 'CorPr_C'),
                     mxAlgebra(abs(Cov_E) / absV, name = 'CorPr_E'),

                     mxAlgebra(cbind(vechs(Cor_A), vechs(Cor_C), vechs(Cor_E),
                                     vechs(Cor),
                                     vechs(CorPr_A), vechs(CorPr_C), vechs(CorPr_E)),
                               dimnames=list(indnames,
                                             c('A', 'C', 'E', 'Total', 'A(%)', 'C(%)', 'E(%)')),
                               name='All_correlations'))

  output_tables(twmodel) <- c('Variance_components', 'All_correlations',
                              'A_correlations', 'C_correlations', 'E_correlations',
                              'Total_correlations')
  twmodel <- def_ci(twmodel, output_tables(twmodel))
  return(twmodel)
}

#' @export
cross_lag_ace <- function(data, zyg, definition = list()) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste0(vars, rep(1:2, each = length(vars)))

  # Starting values
  longdata <-reshape(data[, selvars],
                     varying = selvars,
                     v.names = vars,
                     direction = 'long')
  phmodel <- cross_lag(longdata, definition)
  paths <- rownames(phmodel$Parameter_table)

  phmodel <- mxRun(phmodel)
  # Sterile RAM template
  # The phenotypic model can is used as a template
  template <- mxModel(name = 'template',
                      phmodel$A,
                      phmodel$S,
                      phmodel$I,
                      phmodel$V,
                      phmodel$SDs,
                      phmodel$Sst,
                      phmodel$Ast,
                      phmodel$Parameter_table)

  template$A$values <- phmodel$A$values
  template$S$values <- 0.9 * phmodel$S$values / 3

  mzdata <- data[data[[zyg]] == 'MZ', selvars]
  dzdata <- data[data[[zyg]] == 'DZ', selvars]

  twmodel <- mxModel(name = 'ACE',
                     mxMatrix(type = 'Full',
                              nrow = 1,
                              ncol = nv,
                              free = TRUE,
                              values = phmodel$M$values,
                              name = 'mean'),
                     mxAlgebra(cbind(mean, mean),
                               name = 'expMeans'),
                     mxModel(template,
                             name = 'A'),
                     mxModel(template,
                             name = 'C'),
                     mxModel(template,
                             name = 'E'),
                     mxAlgebra(A.V + C.V + E.V,
                               name = 'V'),

                     mxAlgebra(rbind(cbind(V, A.V + C.V),
                                     cbind(A.V + C.V, V)),
                               name = 'expCovMZ'),
                     mxAlgebra(rbind(cbind(V, 0.5 %x% A.V + C.V),
                                     cbind(0.5 %x% A.V + C.V, V)),
                               name = 'expCovDZ'),

                     mxModel(name = 'MZ',
                             mxData(observed = mzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),
                     mxModel(name = 'DZ',
                             mxData(observed = dzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                                 means = 'ACE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),

                     mxFitFunctionMultigroup(c('MZ','DZ')),

                     # Output tables
                     mxAlgebra(cbind(diag2vec(A.V) / diag2vec(V),
                                     diag2vec(C.V) / diag2vec(V),
                                     diag2vec(E.V) / diag2vec(V),
                                     (diag2vec(A.V) - diag2vec(A.S)) / diag2vec(V),
                                     diag2vec(A.S) / diag2vec(V),
                                     (diag2vec(C.V) - diag2vec(C.S)) / diag2vec(V),
                                     diag2vec(C.S) / diag2vec(V),
                                     (diag2vec(E.V) - diag2vec(E.S)) / diag2vec(V),
                                     diag2vec(E.S) / diag2vec(V)),
                               dimnames = list(vars, c('Atot', 'Ctot', 'Etot',
                                                       'Acom', 'Ares',
                                                       'Ccom', 'Cres',
                                                       'Ecom', 'Eres')),
                               name = 'Variance_components'),
                     mxAlgebra(cbind(A.Parameter_table[,'st estimate'],
                                     C.Parameter_table[,'st estimate'],
                                     E.Parameter_table[,'st estimate']),
                               dimnames = list(paths, c('A', 'C', 'E')),
                               name = 'Path_estimates')
                     )
  output_tables(twmodel) <- c('Variance_components', 'Path_estimates')
  twmodel <- def_ci(twmodel, output_tables(twmodel))
  return(twmodel)
}

#' @export
cross_lag_ade <- function(data, zyg, definition = list()) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste0(vars, rep(1:2, each = length(vars)))

  # Starting values
  longdata <-reshape(data[,selvars],
                     varying = selvars,
                     v.names = vars,
                     direction = 'long')
  phmodel <- cross_lag(longdata, definition)
  paths <- rownames(phmodel$Parameter_table)

  phmodel <- mxRun(phmodel)
  # Sterile RAM template
  # The phenotypic model can is used as a template
  template <- mxModel(name = 'template',
                      phmodel$A,
                      phmodel$S,
                      phmodel$I,
                      phmodel$V,
                      phmodel$SDs,
                      phmodel$Sst,
                      phmodel$Ast,
                      phmodel$Parameter_table)

  template$A$values <- phmodel$A$values
  template$S$values <- 0.9 * phmodel$S$values / 3

  mzdata <- data[data[[zyg]] == 'MZ', selvars]
  dzdata <- data[data[[zyg]] == 'DZ', selvars]

  twmodel <- mxModel(name = 'ADE',
                     mxMatrix(type = 'Full',
                              nrow = 1,
                              ncol = nv,
                              free = TRUE,
                              values = phmodel$M$values,
                              name = 'mean'),
                     mxAlgebra(cbind(mean, mean),
                               name = 'expMeans'),
                     mxModel(template,
                             name = 'A'),
                     mxModel(template,
                             name = 'D'),
                     mxModel(template,
                             name = 'E'),
                     mxAlgebra(A.V + D.V + E.V,
                               name = 'V'),

                     mxAlgebra(rbind(cbind(V, A.V + D.V),
                                     cbind(A.V + D.V, V)),
                               name = 'expCovMZ'),
                     mxAlgebra(rbind(cbind(V, 0.5 %x% A.V + 0.25 %x% D.V),
                                     cbind(0.5 %x% A.V + 0.25 %x% D.V, V)),
                               name = 'expCovDZ'),

                     mxModel(name = 'MZ',
                             mxData(observed = mzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ADE.expCovMZ',
                                                 means = 'ADE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),
                     mxModel(name = 'DZ',
                             mxData(observed = dzdata, type = 'raw'),
                             mxExpectationNormal(covariance = 'ADE.expCovDZ',
                                                 means = 'ADE.expMeans',
                                                 dimnames = selvars),
                             mxFitFunctionML()),

                     mxFitFunctionMultigroup(c('MZ','DZ')),

                     # Output tables
                     mxAlgebra(cbind(diag2vec(A.V) / diag2vec(V),
                                     diag2vec(D.V) / diag2vec(V),
                                     diag2vec(E.V) / diag2vec(V),
                                     (diag2vec(A.V) - diag2vec(A.S)) / diag2vec(V),
                                     diag2vec(A.S) / diag2vec(V),
                                     (diag2vec(D.V) - diag2vec(D.S)) / diag2vec(V),
                                     diag2vec(D.S) / diag2vec(V),
                                     (diag2vec(E.V) - diag2vec(E.S)) / diag2vec(V),
                                     diag2vec(E.S) / diag2vec(V)),
                               dimnames = list(vars, c('Atot', 'Dtot', 'Etot',
                                                       'Acom', 'Ares',
                                                       'Dcom', 'Dres',
                                                       'Ecom', 'Eres')),
                               name = 'Variance_components'),
                     mxAlgebra(cbind(A.Parameter_table[,'st estimate'],
                                     D.Parameter_table[,'st estimate'],
                                     E.Parameter_table[,'st estimate']),
                               dimnames = list(paths, c('A', 'D', 'E')),
                               name = 'Path_estimates')
  )
  output_tables(twmodel) <- c('Variance_components', 'Path_estimates')
  twmodel <- def_ci(twmodel, output_tables(twmodel))
  return(twmodel)
}
