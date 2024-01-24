#' @export
#' @name univ_multiv_ace_ade
#' @title Define univariate and multivariate ACE and ADE models
#'
#' @description
#' The functions define univariate and multivariate (correlated factors) ACE
#' and ADE models.
#'
#' @param data either \code{data.frame} (for raw data) or \code{list} (for
#'        covariance/correlation input). See Note.
#' @param zyg name of the variable that labels zygosity in the
#'        \code{data.frame}.
#' @param ph name of the phenotype, a single character value for univariate
#'        ACE/ADE or a character vector for multivariate ACE/ADE.
#' @param data_type type of the data (raw, cov or cor). See Note.
#' @param sep separator between the name of the phenotype and the label of a
#'        twin. Default is "" (empty string).
#'
#' @return One unfitted \code{mxModel}.
#'
#' @note
#' The function accepts two forms of the data: \code{data.frame}
#' (\code{data_type = 'raw'}) or \code{list} of covariance/correlation matrices
#' (\code{data_type = 'cov'}, \code{data_type = 'cor'}).
#'
#' When \code{data} is a \code{data.frame}, \code{zyg} is expected to point at
#' the variable in \code{data} that defines zygosity groups . Zygosity
#' variable MUST be a factor with two labels: 'MZ' and 'DZ'.
#'
#' When \code{data} is a list of covariance/correlation matrices, it must
#' include two named elements, 'MZ' and 'DZ'. These elements must be the lists
#' with following elements: \code{observed} (covariance/correlation matrix),
#' \code{means} (numeric vector of observed means, optional) and \code{numObs}
#' (number of observations).
#'
#' By default, it is expected that phenotypic trait X is labeled as 'X1' in twin
#' 1 and 'X2' in twin 2.
#'
#' Output tables (univariate ACE/ADE): "Variance components": proportion of
#' total variance explained by A/C/D/E; "Raw variance": total variance,
#' variance of A/C/D/E.
univariate_ace <- function(data,
                           zyg = character(0),
                           ph,
                           data_type = 'raw',
                           sep = getOption('TwinAnalysis.sep')) { # vars is length 1
  if (length(ph) != 1) {
    warning('phenotype has to be of length 1, using only the first phenotypic variable')
    ph <- ph[1]
  }

  nv <- 1
  selvars <- paste(ph, 1:2, sep = sep)

  # Starting values
  startval <- compute_starting_values(
    model = 'univariate_ace',
    data = data,
    data_type = data_type,
    vars = ph,
    sep = sep
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  twmodel <- mxModel(
    name = 'ACE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = startval$freemean,
             values = startval$means,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'a'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'c'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'e'),
    mxAlgebra(t(a) %*% a, name = 'A'),
    mxAlgebra(t(c) %*% c, name = 'C'),
    mxAlgebra(t(e) %*% e, name = 'E'),

    # Covariance
    mxAlgebra(rbind(cbind(A + C + E, A + C    ),
                    cbind(A + C,     A + C + E)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(A + C + E, .5%x%A + C),
                    cbind(.5%x%A + C, A + C + E)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ','DZ')),

    # Output tables
    mxAlgebra(A + C + E, name = 'V'),
    mxAlgebra(cbind(V, A, C, E),
              dimnames = list(ph, c('Total', 'A', 'C', 'E')),
              name = 'Raw_variance'),
    mxAlgebra(cbind(A / V, C / V, E / V),
              dimnames = list(ph, c('A(%)', 'C(%)', 'E(%)')),
              name = 'Variance_components')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components', 'Raw_variance')
  return(twmodel)
}
#' @export
#' @rdname univ_multiv_ace_ade
univariate_ade <- function(data,
                           zyg = character(0),
                           ph,
                           data_type = 'raw',
                           sep = getOption('TwinAnalysis.sep')) { # vars is length 1
  if (length(ph) != 1) {
    warning('phenotype has to be of length 1, using only the first phenotypic variable')
    ph <- ph[1]
  }

  nv <- 1
  selvars <- paste(ph, 1:2, sep = sep)

  # Starting values
  startval <- compute_starting_values(
    model = 'univariate_ace',
    data = data,
    data_type = data_type,
    vars = ph,
    sep = sep
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  twmodel <- mxModel(
    name = 'ADE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = startval$freemean,
             values = startval$means,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'a'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'd'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE, values = sqrt(0.3 * startval$var),
             name = 'e'),
    mxAlgebra(t(a) %*% a, name = 'A'),
    mxAlgebra(t(d) %*% d, name = 'D'),
    mxAlgebra(t(e) %*% e, name = 'E'),

    # Covariance
    mxAlgebra(rbind(cbind(A + D + E, A + D    ),
                    cbind(A + D,     A + D + E)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(A + D + E, .5%x%A + .25%x%D),
                    cbind(.5%x%A + .25%x%D, A + D + E)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ','DZ')),

    # Output tables
    mxAlgebra(A + D + E, name = 'V'),
    mxAlgebra(cbind(V, A, D, E),
              dimnames = list(ph, c('Total', 'A', 'D', 'E')),
              name = 'Raw_variance'),
    mxAlgebra(cbind(A / V, D / V, E / V),
              dimnames = list(ph, c('A(%)', 'D(%)', 'E(%)')),
              name = 'Variance_components')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components', 'Raw_variance')
  return(twmodel)
}

#' @export
#' @rdname univ_multiv_ace_ade
#' @note
#' Output tables (multivariate ACE/ADE): "Variance components": proportion of
#' total variance explained by A/C/D/E, per phenotype; "Covariation
#' components": total covariance, covariance of A/C/D/E, proportion of
#' covariation explained by A/C/D/E; "All correlations": correlations
#' for A/C/D/E and total (phenotypic) correlation; "A/C/D/E/Total
#' correlations": correlation table for A/C/D/E or phenotypic correlation
#' (not returned by \code{get_output_tables}).
multivariate_ace <- function(data,
                             zyg = character(0),
                             ph,
                             data_type = 'raw',
                             sep = getOption('TwinAnalysis.sep')) {
  nv <- length(ph)
  selvars <- paste(ph, rep(1:2, each = nv), sep = sep)

  indnames <- outer(ph, ph,
                    function(x, y) paste0(y, ' <-> ', x))
  indnames <- indnames[row(indnames) >= col(indnames)]

  # Starting values
  startval <- compute_starting_values(
    model = 'multivariate_ace',
    data = data,
    data_type = data_type,
    vars = ph,
    sep = sep
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  twmodel <- mxModel(
    name = 'ACE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = startval$freemean,
             values = startval$means,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'a'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'c'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'e'),
    mxAlgebra(t(a) %*% a, name = 'A'),
    mxAlgebra(t(c) %*% c, name = 'C'),
    mxAlgebra(t(e) %*% e, name = 'E'),

    # Covariance
    mxAlgebra(rbind(cbind(A + C + E, A + C    ),
                    cbind(A + C,     A + C + E)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(A + C + E, .5%x%A + C),
                    cbind(.5%x%A + C, A + C + E)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ','DZ')),

    # Output tables
    mxAlgebra(A + C + E, name = 'V'),
    mxAlgebra(abs(A) + abs(C) + abs(E), name = 'Vabs'),
    mxAlgebra(abs(A) / Vabs, name = 'Aperc'),
    mxAlgebra(abs(C) / Vabs, name = 'Cperc'),
    mxAlgebra(abs(E) / Vabs, name = 'Eperc'),
    mxAlgebra(cbind(diag2vec(Aperc),
                    diag2vec(Cperc),
                    diag2vec(Eperc)),
              dimnames = list(ph, c('A(%)', 'C(%)', 'E(%)')),
              name = 'Variance_components'),

    mxAlgebra(cbind(vech(V),
                    vech(A), vech(C), vech(E),
                    vech(Aperc), vech(Cperc), vech(Eperc)),
              dimnames = list(indnames,
                              c('Total', 'A', 'C', 'E', 'A(%)', 'C(%)', 'E(%)')),
              name = 'Covariation_components'),

    mxMatrix(type = "Iden",
             nrow = nv, ncol = nv,
             name = "I"),
    mxAlgebra(sqrt(I * V), name = "SD_total"),
    mxAlgebra(sqrt(I * A), name = 'SD_A'),
    mxAlgebra(sqrt(I * C), name = 'SD_C'),
    mxAlgebra(sqrt(I * E), name = 'SD_E'),

    mxAlgebra(solve(SD_total) %&% V, name='r_total'),
    mxAlgebra(solve(SD_A) %&% A, name = 'r_A'),
    mxAlgebra(solve(SD_C) %&% C, name = 'r_C'),
    mxAlgebra(solve(SD_E) %&% E, name = 'r_E'),

    mxAlgebra(cbind(vech(r_A), vech(r_C), vech(r_E),
                    vech(r_total)),
              dimnames = list(indnames,
                              c('r_A', 'r_C', 'r_E', 'Total')),
              name='All_correlations'),

    mxAlgebra(r_total,
              dimnames = list(ph, ph),
              name = 'Total_correlations'),
    mxAlgebra(r_A,
              dimnames = list(ph, ph),
              name = 'A_correlations'),
    mxAlgebra(r_C,
              dimnames = list(ph, ph),
              name = 'C_correlations'),
    mxAlgebra(r_E,
              dimnames = list(ph, ph),
              name = 'E_correlations'),

    mxAlgebra(V,
              dimnames = list(ph, ph),
              name = 'Total_covariations'),
    mxAlgebra(A,
              dimnames = list(ph, ph),
              name = 'A_covariations'),
    mxAlgebra(C,
              dimnames = list(ph, ph),
              name = 'C_covariations'),
    mxAlgebra(E,
              dimnames = list(ph, ph),
              name = 'E_covariations')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components', 'Covariation_components', 'All_correlations')
  return(twmodel)
}

#' @export
#' @rdname univ_multiv_ace_ade
multivariate_ade <- function(data,
                             zyg = character(0),
                             ph,
                             data_type = 'raw',
                             sep = getOption('TwinAnalysis.sep')) {
  nv <- length(ph)
  selvars <- paste(ph, rep(1:2, each = nv), sep = sep)

  indnames <- outer(ph, ph,
                    function(x, y) paste0(y, ' <-> ', x))
  indnames <- indnames[row(indnames) >= col(indnames)]

  # Starting values
  startval <- compute_starting_values(
    model = 'multivariate_ace',
    data = data,
    data_type = data_type,
    vars = ph,
    sep = sep
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  twmodel <- mxModel(
    name = 'ADE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = startval$freemean,
             values = startval$means,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'a'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'd'),
    mxMatrix(type = 'Lower',
             nrow = nv, ncol = nv,
             free = TRUE,
             values = 0.6 * startval$chol,
             name = 'e'),
    mxAlgebra(t(a) %*% a, name = 'A'),
    mxAlgebra(t(d) %*% d, name = 'D'),
    mxAlgebra(t(e) %*% e, name = 'E'),

    # Covariance
    mxAlgebra(rbind(cbind(A + D + E, A + D    ),
                    cbind(A + D,     A + D + E)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(A + C + E, .5%x%A + .25%x%D),
                    cbind(.5%x%A + .25%x%D, A + D + E)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ','DZ')),

    # Output tables
    mxAlgebra(A + D + E, name = 'V'),
    mxAlgebra(abs(A) + abs(D) + abs(E), name = 'Vabs'),
    mxAlgebra(abs(A) / Vabs, name = 'Aperc'),
    mxAlgebra(abs(D) / Vabs, name = 'Dperc'),
    mxAlgebra(abs(E) / Vabs, name = 'Eperc'),
    mxAlgebra(cbind(diag2vec(Aperc),
                    diag2vec(Dperc),
                    diag2vec(Eperc)),
              dimnames = list(ph, c('A(%)', 'D(%)', 'E(%)')),
              name = 'Variance_components'),

    mxAlgebra(cbind(vech(V),
                    vech(A), vech(D), vech(E),
                    vech(Aperc), vech(Dperc), vech(Eperc)),
              dimnames = list(indnames,
                              c('Total', 'A', 'D', 'E', 'A(%)', 'D(%)', 'E(%)')),
              name = 'Covariation_components'),

    mxMatrix(type = "Iden",
             nrow = nv, ncol = nv,
             name = "I"),
    mxAlgebra(sqrt(I * V), name = "SD_total"),
    mxAlgebra(sqrt(I * A), name = 'SD_A'),
    mxAlgebra(sqrt(I * D), name = 'SD_D'),
    mxAlgebra(sqrt(I * E), name = 'SD_E'),

    mxAlgebra(solve(SD_total) %&% V, name='r_total'),
    mxAlgebra(solve(SD_A) %&% A, name = 'r_A'),
    mxAlgebra(solve(SD_D) %&% D, name = 'r_D'),
    mxAlgebra(solve(SD_E) %&% E, name = 'r_E'),

    mxAlgebra(cbind(vech(r_A), vech(r_D), vech(r_E),
                    vech(r_total)),
              dimnames = list(indnames,
                              c('r_A', 'r_D', 'r_E', 'Total')),
              name='All_correlations'),

    mxAlgebra(r_total,
              dimnames = list(ph, ph),
              name = 'Total_correlations'),
    mxAlgebra(r_A,
              dimnames = list(ph, ph),
              name = 'A_correlations'),
    mxAlgebra(r_D,
              dimnames = list(ph, ph),
              name = 'C_correlations'),
    mxAlgebra(r_E,
              dimnames = list(ph, ph),
              name = 'E_correlations'),

    mxAlgebra(V,
              dimnames = list(ph, ph),
              name = 'Total_covariations'),
    mxAlgebra(A,
              dimnames = list(ph, ph),
              name = 'A_covariations'),
    mxAlgebra(D,
              dimnames = list(ph, ph),
              name = 'D_covariations'),
    mxAlgebra(E,
              dimnames = list(ph, ph),
              name = 'E_covariations')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components', 'Covariation_components', 'All_correlations')
  return(twmodel)
}


#' @export
#' @name crosslag_ace_ade
#' @title Cross-laged ACE model
#'
#' @description
#' A cross-lagged ACE model is used to decompose phenotypic cross-lagged
#' relationships into ACE. This is the improved version of the model where
#' the algebra that computes phenotypic paths and their decomposition into ACE
#' was improved. ll other components of the model, including core parameters,
#' remained unchanged. The original model is kept as \code{crosslag_ace_old}.
#'
#' @param data either \code{data.frame} (for raw data) or \code{list} for
#'        covariation/correlation input. See Note.
#' @param zyg name of the variable that labels zygosity in the
#'        \code{data.frame}.
#' @param definition a \code{list} that describes measurement occasions.
#'        See Note.
#' @param data_type type of the data (raw, cov or cor). See Note.
#' @param sep separator between the name of the phenotype and the label of a
#'        twin. Default is ''.
#'
#' @return One unfitted \code{mxModel}.
#'
#' @note
#' The function accepts two forms of the data: \code{data.frame}
#' (\code{data_type = 'raw'}) or \code{list} of covariance/correlation matrices
#' (\code{data_type = 'cov'}, \code{data_type = 'cor'}).
#'
#' When \code{data} is a \code{data.frame}, \code{zyg} is expected to point at
#' the variable in \code{data} that defines zygosity groups . Zygosity
#' variable MUST be a factor with two labels: 'MZ' and 'DZ'.
#'
#' When \code{data} is a list of covariance/correlation matrices, it must
#' include two named elements, 'MZ' and 'DZ'. These elements must be the lists
#' with following elements: \code{observed} (covariance/correlation matrix),
#' \code{means} (numeric vector of observed means, optional) and \code{numObs}
#' (number of observations).
#'
#' By default, it is expected that phenotypic trait X is labeled as 'X1' in twin
#' 1 and 'X2' in twin 2.
#'
#' \code{definition} is a list of character vectors. The first character vector
#' includes the variables from the first measurement occasion, the second - the
#' variables from the second measurement occasion, etc. In the model, the
#' variables from the same measurement occasion are assumed to correlate.
#' The variables from one measurement occasion are assumed to predict the
#' variables from the next measurement occasion. Refer to the vignette on
#' genetic cross-lag for more information.
#'
#' Output tables: "Variance components": proportion of variance explained by
#' A/C/D/E, per variable, proportion of variance unaccounted for by preceding
#' measurements (specific variance); "Raw variance": total variance and variance
#' of A/C/D/E, per variable; "Raw paths": unstandardized regression and
#' covariance paths for A/C/D/E; "Standardized paths": standardized regression
#' and covariance paths for A/C/D/E; "Phenotypic paths": total (phenotypic)
#' regression and covariance paths, unstandardized and standardized, and
#' proportion of each path explained by A/C/D/E.
cross_lag_ace <- function(
    data,
    zyg = character(0),
    definition,
    data_type = 'raw',
    sep = getOption('TwinAnalysis.sep')
) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste(vars, rep(1:2, each = nv), sep = sep)

  ## Filter matrices
  niv <- nv - length(rev(definition)[[1]])
  ndv <- nv - length(definition[[1]])
  F1vals <- diag(1, nrow = nv, ncol = niv)
  F2vals <- diag(1, nrow = nv, ncol = ndv)[c((ndv + 1):nv, 1:ndv), ]

  VbyT <- sapply(vars, \(y) sapply(definition, \(x) y %in% x))
  F3vals <- t(VbyT) %*% VbyT
  F4vals <- t(VbyT[-1, ]) %*% VbyT[-nrow(VbyT), ]

  # Starting values
  startval <- compute_starting_values(
    'cross_lag_ace',
    data = data,
    data_type = data_type,
    vars = vars,
    sep = sep,
    definition = definition
  )

  S <- startval$pheno$S
  cholS <- t(chol(S$values))
  s <- mxMatrix(
    type = 'Lower',
    free = S$free[lower.tri(S$free, diag = TRUE)],
    values = cholS[lower.tri(cholS, diag = TRUE)],
    nrow = nrow(cholS),
    ncol = ncol(cholS)
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  # Locate parameters
  params <- omxLocateParameters(startval$pheno)
  params <- params[params$matrix %in% c('A', 'S') &
                     params$row >= params$col, ]
  param_names <- glue_data(params,
                           '{vars[col]}',
                           ' {ifelse(matrix == "A", "-->", "<->")} ',
                           '{vars[row]}')

  twmodel <- mxModel(
    name = 'ACE',

    # Means
    mxMatrix(
      type = 'Full',
      nrow = 1, ncol = nv,
      free = TRUE, values = 0,
      name = 'mean'
    ),
    mxAlgebra(
      cbind(mean, mean),
      name = 'expMeans'
    ),

    # Variance components
    mxMatrix(
      type = 'Iden',
      nrow = nv, ncol = nv,
      name = 'I'
    ),
    lapply(
      c('A', 'C', 'E'),
      function(x) {
        list(
          `$<-`(startval$pheno$A, 'name', glue('A_{x}')),
          `$<-`(s, 'name', glue('s_{x}')),
          mxAlgebraFromString(
            glue('t(s_{x}) %*% s_{x}'),
            name = glue('S_{x}')),
          mxAlgebraFromString(
            glue('solve(I - A_{x}) %&% S_{x}'),
            name = glue('V_{x}'))
        )
      }
    ),

    # Covariance
    mxAlgebra(
      V_A + V_C + V_E,
      name = 'V'
    ),
    mxAlgebra(
      rbind(cbind(V, V_A + V_C),
            cbind(V_A + V_C, V)),
      name = 'expCovMZ'
    ),
    mxAlgebra(
      rbind(cbind(V, 0.5 %x% V_A + V_C),
            cbind(0.5 %x% V_A + V_C, V)),
      name = 'expCovDZ'
    ),

    # MZ submodel
    mxModel(
      name = 'MZ',
      datalist$MZ,
      mxExpectationNormal(
        covariance = 'ACE.expCovMZ',
        means = 'ACE.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),

    # DZ submodel
    mxModel(
      name = 'DZ',
      datalist$DZ,
      mxExpectationNormal(
        covariance = 'ACE.expCovDZ',
        means = 'ACE.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),

    mxFitFunctionMultigroup(
      c('MZ', 'DZ')
    ),

    ## Output tables
    ## Variance components
    mxAlgebra(
      cbind(diag2vec(V_A) / diag2vec(V),
            diag2vec(V_C) / diag2vec(V),
            diag2vec(V_E) / diag2vec(V),
            diag2vec(S_A) / diag2vec(V),
            diag2vec(S_C) / diag2vec(V),
            diag2vec(S_E) / diag2vec(V)
      ),
      dimnames = list(
        vars,
        c('A(%)_tot', 'C(%)_tot', 'E(%)_tot',
          'A(%)_spec', 'C(%)_spec', 'E(%)_spec')
      ),
      name = 'Variance_components'
    ),

    ## Unstandardized variance
    mxAlgebra(
      cbind(diag2vec(V),
            diag2vec(V_A),
            diag2vec(V_C),
            diag2vec(V_E)),
      dimnames = list(vars, c('Total', 'A', 'C', 'E')),
      name = 'Raw_variance'
    ),

    ## Standardized cross-lagged ACE estimates
    lapply(
      c('A', 'C', 'E'),
      function(x) {
        list(
          mxAlgebraFromString(
            glue('sqrt(I * V_{x})'),
            name = glue('SD_{x}')
          ),
          mxAlgebraFromString(
            glue('solve(SD_{x}) %&% S_{x}'),
            name = glue('Sst_{x}')
          ),
          mxAlgebraFromString(
            glue('solve(SD_{x}) %*% A_{x} %*% SD_{x}'),
            name = glue('Ast_{x}')
          )
        )
      }
    ),

    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}st_A[{row}, {col}],
                   {matrix}st_C[{row}, {col}],
                   {matrix}st_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(param_names, c('A', 'C', 'E')),
      name = 'Standardized_paths'
    ),

    ## Unstandardized cross-lagged ACE estimates
    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}_A[{row}, {col}],
                   {matrix}_C[{row}, {col}],
                   {matrix}_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(param_names, c('A', 'C', 'E')),
      name = 'Raw_paths'
    ),

    ## Phenotypic paths

    ## Filter matrices
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = niv,
      free = FALSE,
      values = F1vals,
      name = 'F1'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = ndv,
      free = FALSE,
      values = F2vals,
      name = 'F2'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = nv,
      free = FALSE,
      values = F3vals,
      name = 'F3'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = nv,
      free = FALSE,
      values = F4vals,
      name = 'F4'
    ),

    ## Covariance of independent variables
    mxAlgebra(
      t(F1) %&% (F3 * V),
      name = 'Vx'
    ),

    ## Covariance of dependent variables
    mxAlgebra(
      t(F2) %&% (F3 * V),
      name = 'Vy'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_A),
      name = 'Vy_A'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_C),
      name = 'Vy_C'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_E),
      name = 'Vy_E'
    ),

    ## Covariance between dependent and independent variables
    mxAlgebra(
      t(F2) %*% (F4 * V) %*% F1,
      name = 'Vxy'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_A) %*% F1,
      name = 'Vxy_A'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_C) %*% F1,
      name = 'Vxy_C'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_E) %*% F1,
      name = 'Vxy_E'
    ),

    ## beta
    mxAlgebra(
      Vxy %*% solve(Vx),
      name = 'beta'
    ),
    mxAlgebra(
      Vxy_A %*% solve(Vx),
      name = 'beta_A'
    ),
    mxAlgebra(
      Vxy_C %*% solve(Vx),
      name = 'beta_C'
    ),
    mxAlgebra(
      Vxy_E %*% solve(Vx),
      name = 'beta_E'
    ),

    ## Residual covariance of dependent variables
    mxAlgebra(
      Vy - Vxy %*% t(beta),
      name = 'E'
    ),
    mxAlgebra(
      Vy_A - Vxy_A %*% t(beta),
      name = 'E_A'
    ),
    mxAlgebra(
      Vy_C - Vxy_C %*% t(beta),
      name = 'E_C'
    ),
    mxAlgebra(
      Vy_E - Vxy_E %*% t(beta),
      name = 'E_E'
    ),

    ## Phenotypic A-matrix
    mxAlgebra(
      F2 %*% beta %*% t(F1),
      name = 'A'
    ),
    mxAlgebra(
      abs(F2 %*% beta_A %*% t(F1)) +
        abs(F2 %*% beta_C %*% t(F1)) +
        abs(F2 %*% beta_E %*% t(F1)),
      name = 'absA'
    ),
    mxAlgebra(
      abs(F2 %*% beta_A %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_A'
    ),
    mxAlgebra(
      abs(F2 %*% beta_C %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_C'
    ),
    mxAlgebra(
      abs(F2 %*% beta_E %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_E'
    ),

    ## Phenotypic S-matrix
    mxAlgebra(
      (F3 - F2 %*% t(F2) %*% F3) * V + ## First measurement
        F2 %&% E,                      ## Other measurements
      name = 'S'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_A + F2 %&% E_A) +
        abs((F3 - F2 %*% t(F2) %*% F3) * V_C + F2 %&% E_C) +
        abs((F3 - F2 %*% t(F2) %*% F3) * V_E + F2 %&% E_E),
      name = 'absS'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_A + F2 %&% E_A)
      / (absS + omxNot(absS)),
      name = 'Sperc_A'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_C + F2 %&% E_C)
      / (absS + omxNot(absS)),
      name = 'Sperc_C'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_E + F2 %&% E_E)
      / (absS + omxNot(absS)),
      name = 'Sperc_E'
    ),

    ## Standardized A- and S- matrices
    mxAlgebra(
      sqrt(I * V),
      name = 'SD'
    ),
    mxAlgebra(
      solve(SD) %&% S,
      name = 'Sst'
    ),
    mxAlgebra(
      solve(SD) %*% A %*% SD,
      name = 'Ast'
    ),

    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}[{row}, {col}],
                   {matrix}st[{row}, {col}],
                   {matrix}perc_A[{row}, {col}],
                   {matrix}perc_C[{row}, {col}],
                   {matrix}perc_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(
        param_names,
        c('raw', 'st', 'A(%)', 'C(%)', 'E(%)')
      ),
      name = 'Phenotypic_paths'
    )
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c(
    'Variance_components',
    'Raw_variance',
    'Raw_paths',
    'Standardized_paths',
    'Phenotypic_paths'
  )
  return(twmodel)
}

#' @export
#' @rdname crosslag_ace_ade
cross_lag_ade <- function(
    data,
    zyg = character(0),
    definition,
    data_type = 'raw',
    sep = getOption('TwinAnalysis.sep')
) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste(vars, rep(1:2, each = nv), sep = sep)

  ## Filter matrices
  niv <- nv - length(rev(definition)[[1]])
  ndv <- nv - length(definition[[1]])
  F1vals <- diag(1, nrow = nv, ncol = niv)
  F2vals <- diag(1, nrow = nv, ncol = ndv)[c((ndv + 1):nv, 1:ndv), ]

  VbyT <- sapply(vars, \(y) sapply(definition, \(x) y %in% x))
  F3vals <- t(VbyT) %*% VbyT
  F4vals <- t(VbyT[-1, ]) %*% VbyT[-nrow(VbyT), ]

  # Starting values
  startval <- compute_starting_values(
    'cross_lag_ace',
    data = data,
    data_type = data_type,
    vars = vars,
    sep = sep,
    definition = definition
  )

  S <- startval$pheno$S
  cholS <- t(chol(S$values))
  s <- mxMatrix(
    type = 'Lower',
    free = S$free[lower.tri(S$free, diag = TRUE)],
    values = cholS[lower.tri(cholS, diag = TRUE)],
    nrow = nrow(cholS),
    ncol = ncol(cholS)
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  # Locate parameters
  params <- omxLocateParameters(startval$pheno)
  params <- params[params$matrix %in% c('A', 'S') &
                     params$row >= params$col, ]
  param_names <- glue_data(params,
                           '{vars[col]}',
                           ' {ifelse(matrix == "A", "-->", "<->")} ',
                           '{vars[row]}')

  twmodel <- mxModel(
    name = 'ADE',

    # Means
    mxMatrix(
      type = 'Full',
      nrow = 1, ncol = nv,
      free = TRUE, values = 0,
      name = 'mean'
    ),
    mxAlgebra(
      cbind(mean, mean),
      name = 'expMeans'
    ),

    # Variance components
    mxMatrix(
      type = 'Iden',
      nrow = nv, ncol = nv,
      name = 'I'
    ),
    lapply(
      c('A', 'D', 'E'),
      function(x) {
        list(
          `$<-`(startval$pheno$A, 'name', glue('A_{x}')),
          `$<-`(s, 'name', glue('s_{x}')),
          mxAlgebraFromString(
            glue('t(s_{x}) %*% s_{x}'),
            name = glue('S_{x}')),
          mxAlgebraFromString(
            glue('solve(I - A_{x}) %&% S_{x}'),
            name = glue('V_{x}'))
        )
      }
    ),

    # Covariance
    mxAlgebra(
      V_A + V_D + V_E,
      name = 'V'
    ),
    mxAlgebra(
      rbind(cbind(V, V_A + V_D),
            cbind(V_A + V_D, V)),
      name = 'expCovMZ'
    ),
    mxAlgebra(
      rbind(cbind(V, 0.5 %x% V_A + 0.25 %x% V_D),
            cbind(0.5 %x% V_A + 0.25 %x% V_D, V)),
      name = 'expCovDZ'
    ),

    # MZ submodel
    mxModel(
      name = 'MZ',
      datalist$MZ,
      mxExpectationNormal(
        covariance = 'ADE.expCovMZ',
        means = 'ADE.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),

    # DZ submodel
    mxModel(
      name = 'DZ',
      datalist$DZ,
      mxExpectationNormal(
        covariance = 'ADE.expCovDZ',
        means = 'ADE.expMeans',
        dimnames = selvars
      ),
      mxFitFunctionML()
    ),

    mxFitFunctionMultigroup(
      c('MZ', 'DZ')
    ),

    ## Output tables
    ## Variance components
    mxAlgebra(
      cbind(diag2vec(V_A) / diag2vec(V),
            diag2vec(V_D) / diag2vec(V),
            diag2vec(V_E) / diag2vec(V),
            diag2vec(S_A) / diag2vec(V),
            diag2vec(S_D) / diag2vec(V),
            diag2vec(S_E) / diag2vec(V)
      ),
      dimnames = list(
        vars,
        c('A(%)_tot', 'D(%)_tot', 'E(%)_tot',
          'A(%)_spec', 'D(%)_spec', 'E(%)_spec')
      ),
      name = 'Variance_components'
    ),

    ## Unstandardized variance
    mxAlgebra(
      cbind(diag2vec(V),
            diag2vec(V_A),
            diag2vec(V_D),
            diag2vec(V_E)),
      dimnames = list(vars, c('Total', 'A', 'D', 'E')),
      name = 'Raw_variance'
    ),

    ## Standardized cross-lagged ADE estimates
    lapply(
      c('A', 'D', 'E'),
      function(x) {
        list(
          mxAlgebraFromString(
            glue('sqrt(I * V_{x})'),
            name = glue('SD_{x}')
          ),
          mxAlgebraFromString(
            glue('solve(SD_{x}) %&% S_{x}'),
            name = glue('Sst_{x}')
          ),
          mxAlgebraFromString(
            glue('solve(SD_{x}) %*% A_{x} %*% SD_{x}'),
            name = glue('Ast_{x}')
          )
        )
      }
    ),

    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}st_A[{row}, {col}],
                   {matrix}st_D[{row}, {col}],
                   {matrix}st_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(param_names, c('A', 'D', 'E')),
      name = 'Standardized_paths'
    ),

    ## Unstandardized cross-lagged ACE estimates
    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}_A[{row}, {col}],
                   {matrix}_D[{row}, {col}],
                   {matrix}_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(param_names, c('A', 'D', 'E')),
      name = 'Raw_paths'
    ),

    ## Phenotypic paths

    ## Filter matrices
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = niv,
      free = FALSE,
      values = F1vals,
      name = 'F1'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = ndv,
      free = FALSE,
      values = F2vals,
      name = 'F2'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = nv,
      free = FALSE,
      values = F3vals,
      name = 'F3'
    ),
    mxMatrix(
      type = 'Full',
      nrow = nv,
      ncol = nv,
      free = FALSE,
      values = F4vals,
      name = 'F4'
    ),

    ## Covariance of independent variables
    mxAlgebra(
      t(F1) %&% (F3 * V),
      name = 'Vx'
    ),

    ## Covariance of dependent variables
    mxAlgebra(
      t(F2) %&% (F3 * V),
      name = 'Vy'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_A),
      name = 'Vy_A'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_D),
      name = 'Vy_D'
    ),
    mxAlgebra(
      t(F2) %&% (F3 * V_E),
      name = 'Vy_E'
    ),

    ## Covariance between dependent and independent variables
    mxAlgebra(
      t(F2) %*% (F4 * V) %*% F1,
      name = 'Vxy'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_A) %*% F1,
      name = 'Vxy_A'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_D) %*% F1,
      name = 'Vxy_D'
    ),
    mxAlgebra(
      t(F2) %*% (F4 * V_E) %*% F1,
      name = 'Vxy_E'
    ),

    ## beta
    mxAlgebra(
      Vxy %*% solve(Vx),
      name = 'beta'
    ),
    mxAlgebra(
      Vxy_A %*% solve(Vx),
      name = 'beta_A'
    ),
    mxAlgebra(
      Vxy_D %*% solve(Vx),
      name = 'beta_D'
    ),
    mxAlgebra(
      Vxy_E %*% solve(Vx),
      name = 'beta_E'
    ),

    ## Residual covariance of dependent variables
    mxAlgebra(
      Vy - Vxy %*% t(beta),
      name = 'E'
    ),
    mxAlgebra(
      Vy_A - Vxy_A %*% t(beta),
      name = 'E_A'
    ),
    mxAlgebra(
      Vy_D - Vxy_D %*% t(beta),
      name = 'E_D'
    ),
    mxAlgebra(
      Vy_E - Vxy_E %*% t(beta),
      name = 'E_E'
    ),

    ## Phenotypic A-matrix
    mxAlgebra(
      F2 %*% beta %*% t(F1),
      name = 'A'
    ),
    mxAlgebra(
      abs(F2 %*% beta_A %*% t(F1)) +
        abs(F2 %*% beta_D %*% t(F1)) +
        abs(F2 %*% beta_E %*% t(F1)),
      name = 'absA'
    ),
    mxAlgebra(
      abs(F2 %*% beta_A %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_A'
    ),
    mxAlgebra(
      abs(F2 %*% beta_D %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_D'
    ),
    mxAlgebra(
      abs(F2 %*% beta_E %*% t(F1)) / (absA + omxNot(absA)),
      name = 'Aperc_E'
    ),

    ## Phenotypic S-matrix
    mxAlgebra(
      (F3 - F2 %*% t(F2) %*% F3) * V + ## First measurement
        F2 %&% E,                      ## Other measurements
      name = 'S'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_A + F2 %&% E_A) +
        abs((F3 - F2 %*% t(F2) %*% F3) * V_D + F2 %&% E_D) +
        abs((F3 - F2 %*% t(F2) %*% F3) * V_E + F2 %&% E_E),
      name = 'absS'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_A + F2 %&% E_A)
      / (absS + omxNot(absS)),
      name = 'Sperc_A'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_D + F2 %&% E_D)
      / (absS + omxNot(absS)),
      name = 'Sperc_C'
    ),
    mxAlgebra(
      abs((F3 - F2 %*% t(F2) %*% F3) * V_E + F2 %&% E_E)
      / (absS + omxNot(absS)),
      name = 'Sperc_E'
    ),

    ## Standardized A- and S- matrices
    mxAlgebra(
      sqrt(I * V),
      name = 'SD'
    ),
    mxAlgebra(
      solve(SD) %&% S,
      name = 'Sst'
    ),
    mxAlgebra(
      solve(SD) %*% A %*% SD,
      name = 'Ast'
    ),

    mxAlgebraFromString(
      glue(
        'rbind(',
        paste(
          glue_data(
            params,
            'cbind({matrix}[{row}, {col}],
                   {matrix}st[{row}, {col}],
                   {matrix}perc_A[{row}, {col}],
                   {matrix}perc_D[{row}, {col}],
                   {matrix}perc_E[{row}, {col}])'
          ),
          collapse = ', '
        ),
        ')'
      ),
      dimnames = list(
        param_names,
        c('raw', 'st', 'A(%)', 'D(%)', 'E(%)')
      ),
      name = 'Phenotypic_paths'
    )
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c(
    'Variance_components',
    'Raw_variance',
    'Raw_paths',
    'Standardized_paths',
    'Phenotypic_paths'
  )
  return(twmodel)
}

#' @export
#' @name crosslag_ace_ade_old
#' @title The cross-lagged ACE model used before 2024
#'
#' @description
#' This is the cross-lagged ACE model used before 2024. The new model which is
#' now under the name \code{crosslag_ace} is has an improved algebra for
#' computation of phenotypic paths and their decomposition into ACE. All other
#' components of the model, including the core parameters, remained unchanged.
#'
#' The onld model is kept for reproducibility.
#'
#' @param data either \code{data.frame} (for raw data) or \code{list} for
#'        covariation/correlation input. See Note.
#' @param zyg name of the variable that labels zygosity in the
#'        \code{data.frame}.
#' @param definition a \code{list} that describes measurement occasions.
#'        See Note.
#' @param data_type type of the data (raw, cov or cor). See Note.
#' @param sep separator between the name of the phenotype and the label of a
#'        twin. Default is ''.
#'
#' @return One unfitted \code{mxModel}.
#'
#' @note
#' The function accepts two forms of the data: \code{data.frame}
#' (\code{data_type = 'raw'}) or \code{list} of covariance/correlation matrices
#' (\code{data_type = 'cov'}, \code{data_type = 'cor'}).
#'
#' When \code{data} is a \code{data.frame}, \code{zyg} is expected to point at
#' the variable in \code{data} that defines zygosity groups . Zygosity
#' variable MUST be a factor with two labels: 'MZ' and 'DZ'.
#'
#' When \code{data} is a list of covariance/correlation matrices, it must
#' include two named elements, 'MZ' and 'DZ'. These elements must be the lists
#' with following elements: \code{observed} (covariance/correlation matrix),
#' \code{means} (numeric vector of observed means, optional) and \code{numObs}
#' (number of observations).
#'
#' By default, it is expected that phenotypic trait X is labeled as 'X1' in twin
#' 1 and 'X2' in twin 2.
#'
#' \code{definition} is a list of character vectors. The first character vector
#' includes the variables from the first measurement occasion, the second - the
#' variables from the second measurement occasion, etc. In the model, the
#' variables from the same measurement occasion are assumed to correlate.
#' The variables from one measurement occasion are assumed to predict the
#' variables from the next measurement occasion. Refer to the vignette on
#' genetic cross-lag for more information.
#'
#' Output tables: "Variance components": proportion of variance explained by
#' A/C/D/E, per variable, proportion of variance unaccounted for by preceding
#' measurements (specific variance); "Raw variance": total variance and variance
#' of A/C/D/E, per variable; "Raw paths": unstandardized regression and
#' covariance paths for A/C/D/E; "Standardized paths": standardized regression
#' and covariance paths for A/C/D/E; "Phenotypic paths": total (phenotypic)
#' regression and covariance paths, unstandardized and standardized, and
#' proportion of each path explained by A/C/D/E.
cross_lag_ace_old <- function(
    data,
    zyg = character(0),
    definition,
    data_type = 'raw',
    sep = getOption('TwinAnalysis.sep')
) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste(vars, rep(1:2, each = nv), sep = sep)

  # Starting values
  startval <- compute_starting_values('cross_lag_ace',
                                      data = data,
                                      data_type = data_type,
                                      vars = vars,
                                      sep = sep,
                                      definition = definition)

  S <- startval$pheno$S
  cholS <- t(chol(S$values))
  s <- mxMatrix(
    type = 'Lower',
    free = S$free[lower.tri(S$free, diag = TRUE)],
    values = cholS[lower.tri(cholS, diag = TRUE)],
    nrow = nrow(cholS),
    ncol = ncol(cholS)
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  # Locate parameters
  params <- omxLocateParameters(startval$pheno)
  params <- params[params$matrix %in% c('A', 'S') &
                     params$row >= params$col, ]
  param_names <- glue_data(params,
                           '{vars[col]}',
                           ' {ifelse(matrix == "A", "-->", "<->")} ',
                           '{vars[row]}')

  twmodel <- mxModel(
    name = 'ACE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = TRUE, values = 0,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Iden',
             nrow = nv, ncol = nv,
             name = 'I'),
    lapply(c('A', 'C', 'E'),
           function(x) {
             list(
               `$<-`(startval$pheno$A, 'name', glue('A_{x}')),
               `$<-`(s, 'name', glue('s_{x}')),
               mxAlgebraFromString(
                 glue('t(s_{x}) %*% s_{x}'),
                 name = glue('S_{x}')),
               mxAlgebraFromString(
                 glue('solve(I - A_{x}) %&% S_{x}'),
                 name = glue('V_{x}'))
             )
           }),

    # Covariance
    mxAlgebra(V_A + V_C + V_E, name = 'V_tot'),
    mxAlgebra(rbind(cbind(V_tot, V_A + V_C),
                    cbind(V_A + V_C, V_tot)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(V_tot, 0.5 %x% V_A + V_C),
                    cbind(0.5 %x% V_A + V_C, V_tot)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ', 'DZ')),

    # Output tables
    mxAlgebra(cbind(diag2vec(V_A) / diag2vec(V_tot),
                    diag2vec(V_C) / diag2vec(V_tot),
                    diag2vec(V_E) / diag2vec(V_tot),
                    diag2vec(S_A) / diag2vec(V_tot),
                    diag2vec(S_C) / diag2vec(V_tot),
                    diag2vec(S_E) / diag2vec(V_tot)),
              dimnames = list(vars, c('A(%)_tot', 'C(%)_tot', 'E(%)_tot',
                                      'A(%)_spec', 'C(%)_spec', 'E(%)_spec')),
              name = 'Variance_components'),

    mxAlgebra(cbind(diag2vec(V_tot),
                    diag2vec(V_A),
                    diag2vec(V_C),
                    diag2vec(V_E)),
              dimnames = list(vars, c('Total', 'A', 'C', 'E')),
              name = 'Raw_variance'),

    mxAlgebra(sqrt(I * V_tot),
              name = 'SD'),
    lapply(c('A', 'C', 'E'),
           function(x) {
             list(
               mxAlgebraFromString(glue('sqrt(I * V_{x})'),
                                   name = glue('SD_{x}')),
               mxAlgebraFromString(glue('solve(SD_{x}) %&% S_{x}'),
                                   name = glue('Sst_{x}')),
               mxAlgebraFromString(glue('solve(SD_{x}) %*% A_{x} %*% SD_{x}'),
                                   name = glue('Ast_{x}')),
               mxAlgebraFromString(glue('solve(SD) %&% S_{x}'),
                                   name = glue('Sst2_{x}')),
               mxAlgebraFromString(glue('solve(SD) %*% A_{x} %*% SD'),
                                   name = glue('Ast2_{x}')),
               mxAlgebraFromString(glue('abs(S_{x})'),
                                   name = glue('Sabs_{x}')),
               mxAlgebraFromString(glue('abs(A_{x} %*% (I * V_{x}))'),
                                   name = glue('Aabs_{x}'))
             )
           }),

    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}st_A[{row}, {col}],
                             {matrix}st_C[{row}, {col}],
                             {matrix}st_E[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('A', 'C', 'E')),
      name = 'Standardized_paths'),

    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}_A[{row}, {col}],
                             {matrix}_C[{row}, {col}],
                             {matrix}_E[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('A', 'C', 'E')),
      name = 'Raw_paths'),

    mxAlgebra(A_A %*% (I * V_A) + A_C %*% (I * V_C) + A_E %*% (I * V_E),
              name = 'A_ph'),
    mxAlgebra(Aabs_A + Aabs_C + Aabs_E,
              name = 'Aabs_ph'),
    mxAlgebra(S_A + S_C + S_E,
              name = 'S_ph'),
    mxAlgebra(Sabs_A + Sabs_C + Sabs_E,
              name = 'Sabs_ph'),
    mxAlgebra(solve(SD) %*% A_ph %*% SD,
              name = 'Ast_ph'),
    mxAlgebra(solve(SD) %&% S_ph,
              name = 'Sst_ph'),
    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}_ph[{row}, {col}],
                             {matrix}st_ph[{row}, {col}],
                            ({matrix}abs_A/{matrix}abs_ph)[{row}, {col}],
                            ({matrix}abs_C/{matrix}abs_ph)[{row}, {col}],
                            ({matrix}abs_E/{matrix}abs_ph)[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('raw', 'st', 'A(%)', 'C(%)', 'E(%)')),
      name = 'Phenotypic_paths')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components',
                              'Raw_variance',
                              'Raw_paths',
                              'Standardized_paths',
                              'Phenotypic_paths')
  return(twmodel)
}

#' @export
#' @rdname crosslag_ace_ade_old
cross_lag_ade_old <- function(data,
                          zyg = character(0),
                          definition,
                          data_type = 'raw',
                          sep = getOption('TwinAnalysis.sep')) {
  vars <- unlist(definition)
  nv <- length(vars)
  selvars <- paste(vars, rep(1:2, each = nv), sep = sep)

  # Starting values
  startval <- compute_starting_values('cross_lag_ace',
                                      data = data,
                                      data_type = data_type,
                                      vars = vars,
                                      sep = sep,
                                      definition = definition)

  S <- startval$pheno$S
  cholS <- t(chol(S$values))
  s <- mxMatrix(
    type = 'Lower',
    free = S$free[lower.tri(S$free, diag = TRUE)],
    values = cholS[lower.tri(cholS, diag = TRUE)],
    nrow = nrow(cholS),
    ncol = ncol(cholS)
  )

  # Process twin data into MxData
  datalist <- process_twin_data(
    data = data,
    selvars = selvars,
    data_type = data_type,
    zyg = zyg
  )

  # Locate parameters
  params <- omxLocateParameters(startval$pheno)
  params <- params[params$matrix %in% c('A', 'S') &
                     params$row >= params$col, ]
  param_names <- glue_data(params,
                           '{vars[col]}',
                           ' {ifelse(matrix == "A", "-->", "<->")} ',
                           '{vars[row]}')

  twmodel <- mxModel(
    name = 'ADE',

    # Means
    mxMatrix(type = 'Full',
             nrow = 1, ncol = nv,
             free = TRUE, values = 0,
             name = 'mean'),
    mxAlgebra(cbind(mean, mean),
              name = 'expMeans'),

    # Variance components
    mxMatrix(type = 'Iden',
             nrow = nv, ncol = nv,
             name = 'I'),
    lapply(c('A', 'D', 'E'),
           function(x) {
             list(
               `$<-`(startval$pheno$A, 'name', glue('A_{x}')),
               `$<-`(s, 'name', glue('s_{x}')),
               mxAlgebraFromString(
                 glue('t(s_{x}) %*% s_{x}'),
                 name = glue('S_{x}')),
               mxAlgebraFromString(
                 glue('solve(I - A_{x}) %&% S_{x}'),
                 name = glue('V_{x}'))
             )
           }),

    # Covariance
    mxAlgebra(V_A + V_D + V_E, name = 'V_tot'),
    mxAlgebra(rbind(cbind(V_tot, V_A + V_D),
                    cbind(V_A + V_D, V_tot)),
              name = 'expCovMZ'),
    mxAlgebra(rbind(cbind(V_tot, 0.5 %x% V_A + 0.25 %x% V_D),
                    cbind(0.5 %x% V_A + 0.25 %x% V_D, V_tot)),
              name = 'expCovDZ'),

    # MZ submodel
    mxModel(name = 'MZ',
            datalist$MZ,
            mxExpectationNormal(covariance = 'ACE.expCovMZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    # DZ submodel
    mxModel(name = 'DZ',
            datalist$DZ,
            mxExpectationNormal(covariance = 'ACE.expCovDZ',
                                means = 'ACE.expMeans',
                                dimnames = selvars),
            mxFitFunctionML()),

    mxFitFunctionMultigroup(c('MZ', 'DZ')),

    # Output tables
    mxAlgebra(cbind(diag2vec(V_A) / diag2vec(V_tot),
                    diag2vec(V_D) / diag2vec(V_tot),
                    diag2vec(V_E) / diag2vec(V_tot),
                    diag2vec(S_A) / diag2vec(V_tot),
                    diag2vec(S_D) / diag2vec(V_tot),
                    diag2vec(S_E) / diag2vec(V_tot)),
              dimnames = list(vars, c('A(%)_tot', 'D(%)_tot', 'E(%)_tot',
                                      'A(%)_spec', 'D(%)_spec', 'E(%)_spec')),
              name = 'Variance_components'),

    mxAlgebra(cbind(diag2vec(V_tot),
                    diag2vec(V_A),
                    diag2vec(V_D),
                    diag2vec(V_E)),
              dimnames = list(vars, c('Total', 'A', 'D', 'E')),
              name = 'Raw_variance'),

    mxAlgebra(sqrt(I * V_tot),
              name = 'SD'),
    lapply(c('A', 'D', 'E'),
           function(x) {
             list(
               mxAlgebraFromString(glue('sqrt(I * V_{x})'),
                                   name = glue('SD_{x}')),
               mxAlgebraFromString(glue('solve(SD_{x}) %&% S_{x}'),
                                   name = glue('Sst_{x}')),
               mxAlgebraFromString(glue('solve(SD_{x}) %*% A_{x} %*% SD_{x}'),
                                   name = glue('Ast_{x}')),
               mxAlgebraFromString(glue('solve(SD) %&% S_{x}'),
                                   name = glue('Sst2_{x}')),
               mxAlgebraFromString(glue('solve(SD) %*% A_{x} %*% SD'),
                                   name = glue('Ast2_{x}')),
               mxAlgebraFromString(glue('abs(S_{x})'),
                                   name = glue('Sabs_{x}')),
               mxAlgebraFromString(glue('abs(A_{x} %*% (I * V_{x}))'),
                                   name = glue('Aabs_{x}'))
             )
           }),

    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}st_A[{row}, {col}],
                             {matrix}st_D[{row}, {col}],
                             {matrix}st_E[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('A', 'D', 'E')),
      name = 'Standardized_paths'),

    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}_A[{row}, {col}],
                             {matrix}_D[{row}, {col}],
                             {matrix}_E[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('A', 'D', 'E')),
      name = 'Raw_paths'),

    mxAlgebra(A_A %*% (I * V_A) + A_D %*% (I * V_D) + A_E %*% (I * V_E),
              name = 'A_ph'),
    mxAlgebra(Aabs_A + Aabs_D + Aabs_E,
              name = 'Aabs_ph'),
    mxAlgebra(S_A + S_D + S_E,
              name = 'S_ph'),
    mxAlgebra(Sabs_A + Sabs_D + Sabs_E,
              name = 'Sabs_ph'),
    mxAlgebra(solve(SD) %*% A_ph %*% SD,
              name = 'Ast_ph'),
    mxAlgebra(solve(SD) %&% S_ph,
              name = 'Sst_ph'),
    mxAlgebraFromString(glue(
      'rbind(',
      paste(glue_data(params,
                      'cbind({matrix}_ph[{row}, {col}],
                             {matrix}st_ph[{row}, {col}],
                            ({matrix}abs_A/{matrix}abs_ph)[{row}, {col}],
                            ({matrix}abs_D/{matrix}abs_ph)[{row}, {col}],
                            ({matrix}abs_E/{matrix}abs_ph)[{row}, {col}])'),
            collapse = ', '),
      ')'),
      dimnames = list(param_names, c('raw', 'st', 'A(%)', 'D(%)', 'E(%)')),
      name = 'Phenotypic_paths')
  )

  twmodel <- assign_auto_labels(twmodel)
  output_tables(twmodel) <- c('Variance_components',
                              'Raw_variance',
                              'Raw_paths',
                              'Standardized_paths',
                              'Phenotypic_paths')
  return(twmodel)
}
