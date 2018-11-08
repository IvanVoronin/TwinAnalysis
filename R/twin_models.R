#' @export
cross_lag_ace <- function(data, zyg, definition = list()) {
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
  return(twmodel)
}
