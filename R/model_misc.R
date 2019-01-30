#' @import OpenMx
#' @import mlth.data.frame

#' @rdname output_tables
#' @title Define and extract \code{MxModel} output tables.
#'
#' @description
#' The functions to define and extract the output tables for \code{MxModels}.
#' The names of the output tables are stored as the \code{output_tables}
#' option of an \code{MxModel} object.
#'
#' @param model \code{MxModel} object.
#' @param output_tables,value character vector, the names of the output tables
#'        defined inside the model as \code{MxAlgebra}.
#'
#' @return \code{def_output_tables} returns an \code{MxModel}
#'         (\code{`output_tables<-`} is an alias).
#'
#'         \code{output_table} returns \code{character} vector of the names of
#'         the output tables.
#'
#'         \code{get_output_tables} returns the list of the evaluated output
#'         tables. If the model output includes confidence intervals, it append
#'         the intervals. Each table is either \code{data.frame} or
#'         \code{mlth.data.frame}.
#'
#' @examples
#' @export
def_output_tables <- function(model, output_tables = character(0)) {
  # The names of the output tables are stored in the options slot of
  # MxModel because this way they stay unchanged after the model
  # has been fit
  if (length(output_tables > 0))
  if (!any(output_tables %in% names(model)))
    stop('some output tables are not defined in the model. check names(model)')

  model@options$output_tables <- output_tables
  return(model)
}

#' @rdname output_tables
#' @export
`output_tables<-` <- function(model, value) {
  def_output_tables(model, value)
}

#' @rdname output_tables
#' @export
get_output_tables <- function(model) {
  output_tables <- output_tables(model)
  pars <- lapply(output_tables, mxEvalByName, model)
  names(pars) <- output_tables
  cis <- get_ci(model)

  outp <- list()
  for (i in output_tables) {
    if (length(cis[[i]]) == 0) {
      outp[[i]] <- pars[[i]]
    } else {
      outp[[i]] <- as.mlth.data.frame(Map(list,
                                          est = as.data.frame(pars[[i]]),
                                          lCI = as.data.frame(cis[[i]]$lbound),
                                          uCI = as.data.frame(cis[[i]]$ubound)))
      names(outp[[i]]) <- colnames(pars[[i]])
      row.names(outp[[i]]) <- rownames(pars[[i]])
    }
  }
  return(outp)
}


#' @rdname output_tables
#' @export
output_tables <- function(model){
  model@options$output_tables
}

#' @rdname ref_models
#' @title Define saturared and independence models.
#'
#' @description The functions to define, fit and store the reference models.
#' The reference models are stored in the \code{options} slot of the target
#' model.
#'
#' @param model \code{MxModel} object
#' @param ref_models the list of length two of \code{MxModel}s
#'        (saturated and independence) or \code{NULL}.
#' @param ... arguments passed to \code{mxRefModels} or \code{mxRun}
#'        (e.g., \code{run = TRUE}).
#'
#' @return \code{def_ref_models} returns the \code{MxModel}
#'         (\code{`ref_models<-`} is an alias)).
#'
#'         \code{ref_models} returns the list of \code{MxModel}s.
#'
#'         \code{run_ref_models} runs the parameter estimation for the
#'         reference models and return the \code{MxModel}.
#'
#' @note
#' The idea is to store the reference models together with the original model,
#' so the function \code{fit_stats} could address them and compute fit
#' statistics. The user can define the reference models and pass them to the
#' function (as a list of two \code{MxModels}). When the argument
#' \code{ref_models} is empty or there is no Saturated or Independence model
#' in it, the function makes use of \code{mxRefModels}.
#'
#' Don't forget to re-define the reference models when you modify the
#' model's data.
#'
#' @examples
#' @export
def_ref_models <- function(model,
                           ref_models = NULL,
                           ...) {
  if (length(ref_models) == 0)
    ref_models <- mxRefModels(model, ...)
  # TODO: Check that there are only two models: Saturated and Independence
  model@options$ref_models <- ref_models
  return(model)
}

#' @rdname ref_models
#' @export
`ref_models<-` <- function(model, value) {
  def_ref_models(model, value)
}

#' @rdname ref_models
#' @export
ref_models <- function(model) {
  model@options$ref_models
}

#' @rdname ref_models
#' @export
run_ref_models <- function(model, ...) {
  sat_model <- ref_models(model)$Saturated
  sat_model <- mxRun(sat_model, ...)
  ind_model <- ref_models(model)$Independence
  ind_model <- mxRun(ind_model, ...)
  ref_models(model) <- list(Saturated = sat_model, Independence = ind_model)
  return(model)
}

#' @rdname fit_stats
#' @title Compute fit statistics.
#'
#' @description
#' This function assemble the model parameters and fit statistics into
#' a table. It also compares nested models.
#'
#' @param model \code{MxModel} with the reference models included.
#' @param nested_models the list of nested models (\code{MxModel}s) to be
#'        compared against the target model.
#'
#' @return \code{data.frame}
#'
#' @note
#' The function uses \code{summary.MxModel} and \code{mxCompare} to assembe
#' the following parameters and statistics:
#' \itemize{
#'    \item{ep}~--- number of estimated parameters;
#'    \item{minus2LL}~--- \eqn{-2*log-likelihood};
#'    \item{df}~--- degrees of freedom;
#'    \item{AIC}~--- Akaike Information Criterion;
#'    \item{BIC}~--- Bayesian Information Criterion;
#'    \item{CFI}~--- Comparative Fit Index;
#'    \item{TLI}~--- Tucker-Lewis Index;
#'    \item{RMSEA}~--- Root Mean Square Error of Approximation;
#'    \item{diffLL}~--- log-likelihood difference between two models,
#'       chi-squared;
#'    \item{diffdf}~--- the difference in degrees of freedom between
#'       two models, chi-square df;
#'    \item{p}~--- p-value of the chi-square test.
#' }
#'
#' The references will be added.
#'
#' The nested models are compared against the target model by means of the
#' chi-square test (the same way the target model is compared against the
#' saturated model).
#'
#' This function does not perform any computation, it is a wrapper for the
#' \code{OpenMx} functions aiming to simplify the preparation of the output
#' from SEM.
#'
#' @examples
#' @export
fit_stats <- function(model, nested_models = NULL) {
  refmodels <- ref_models(model)
  modelsum <- summary(model, refModels = refmodels)
  outp <- data.frame(model    = model@name,
                     base     = 'Saturated',
                     ep       = modelsum$estimatedParameters,
                     minus2LL = modelsum$Minus2LogLikelihood,
                     df       = modelsum$degreesOfFreedom,
                     AIC      = modelsum$AIC.Mx,
                     BIC      = modelsum$BIC.Mx,
                     CFI      = modelsum$CFI,
                     TLI      = modelsum$TLI,
                     RMSEA    = modelsum$RMSEA,
                     diffLL   = modelsum$Chi,
                     diffdf   = modelsum$ChiDoF,
                     p        = modelsum$p)
  if (length(nested_models) > 0) {
    if (length(names(nested_models)) == 0)
      names(nested_models) <- paste('Model', 1:length(nested_models))
    nested_table <- mxCompare(model, nested_models)
    for (i in 1:length(nested_models)) {
      modelsum <- summary(nested_models[[i]], refModels = refmodels)
      outp <- rbind(outp,
                    data.frame(model    = names(nested_models)[i],
                               base     = model@name,
                               ep       = modelsum$estimatedParameters,
                               minus2LL = modelsum$Minus2LogLikelihood,
                               df       = modelsum$degreesOfFreedom,
                               AIC      = modelsum$AIC.Mx,
                               BIC      = modelsum$BIC.Mx,
                               CFI      = modelsum$CFI,
                               TLI      = modelsum$TLI,
                               RMSEA    = modelsum$RMSEA,
                               diffLL   = nested_table[i+1, 'diffLL'],
                               diffdf   = nested_table[i+1, 'diffdf'],
                               p        = nested_table[i+1, 'p']))
    }
  }
  return(outp)
}

#' @title Define nested models.
#'
#' @description
#' Define nested models by adding \code{MxConstraint}s or modifying the labels
#' or other elements of the target \code{MxModel}.
#'
#' @param model \code{MxModel} object.
#' @param ... named \code{OpenMx} objects to add or replace in the target
#'        model (in a list if there are more than one).
#' @param run should the models be run? Default is \code{FALSE}.
#'
#' @return The named list of \code{MxModel}s.
#'
#' @examples
#' @export
def_nested_models <- function(model, ..., run = FALSE) {
  # TODO: Fix the bug when the nested models names lacking
  dots <- list(...)
  if (run)
    lapply(dots, function(x) mxRun(mxModel(model, x)))
  else
    lapply(dots, function(x) mxModel(model, x))
}

#' @rdname ci
#' @title Define confidence intervals.
#'
#' @description
#' The function to define and extract confidence intervals. This information
#' is stored in the \code{options} slot of an \code{MxModel}.
#'
#' @param model \code{MxModel} object
#' @param ci the character vector pointing out the parameters to
#'           compute the intervals for, see \code{\link[OpenMx]{mxCI}}.
#' @param ... the arguments passed to mxCI.
#'
#' @return \code{MxModel}
#'
#' @note
#' These functions aim to facilitate the assembly of parameters and their CIs.
#' \code{get_ci} is used by \code{\link{get_output_tables}}.
#'
#' @examples
#'
#' @export
def_ci <- function(model, ci = character(0), ...) {
  if (length(ci) > 0)
    mxModel(model, mxCI(ci, ...))
}

#' @rdname ci
#' @export
get_ci <- function(model, ci = output_tables(model)) {
  cis <- model@output$confidenceIntervals
  if (length(ci) == 0 || length(cis) == 0)
    return(NULL)

  outp <- list()
  for (i in ci) {
    cinames <- grep(i, rownames(cis), value = TRUE)
    if (length(cinames) > 0) {
      pars <- mxEvalByName(i, model)
      lci <- uci <- matrix(NA, nrow = nrow(pars), ncol = ncol(pars),
                           dimnames = dimnames(pars))
      ind <- expand.grid(1:nrow(pars), 1:ncol(pars))
      for (j in 1:nrow(ind)) {
        ciname <- grep(sprintf('%s[%d,%d]', i, ind[j, 1], ind[j, 2]),
                       cinames, value = TRUE, fixed = TRUE)
        if (length(ciname) > 0){
          lci[ind[j, 1], ind[j, 2]] <- cis[ciname, 'lbound']
          uci[ind[j, 1], ind[j, 2]] <- cis[ciname, 'ubound']
        }
      }
      outp[[i]] <- list(lbound = lci, ubound = uci)
    }
  }
  return(outp)
}

#' @title Define standardized parameters for the RAM model.
#'
#' @description
#' The function adds the algebra defining the standardized parameters of a
#' RAM model and a fancy table with raw and standardized estimates.
#'
#' @param model \code{MxModel} object of type RAM.
#'
#' @return \code{MxModel} object.
#'
#' @note
#' The standardization is based on the variances from the expected covariation
#' matrix of the model. One-way paths are standardized as
#' \eqn{b_XY * \sqrt{V_Y}/\sqrt{V_X}} (CHECK!!!). Two-way paths are
#' standardized as \eqn{Cov_XY / (\sqrt{V_X}\sqrt{V_Y})}.
#'
#' The function defines two output tables: \code{Parameter_table} includes
#' one-way and two-way paths and \code{Means} include the means. The user can
#' identify the confidence intervals for the whole table or for its part,
#' see examples for \code{\link{def_ci}}.
#'
#' @examples
#'
#' @export
def_stand_params <- function(model) {
  stopifnot(class(model) == 'MxRAMModel')
  vars <- c(c(model@manifestVars, model@latentVars))
  nv <- length(vars)

  # Standardize parameters.
  model <- mxModel(
    model,
    mxMatrix(
      type = 'Iden',
      ncol = nv,
      nrow = nv,
      name = 'I'
    ),
    mxAlgebra(solve(I - A) %&% S, name =
                'V'),
    mxAlgebra(sqrt(I * V), name = 'SDs'),
    mxAlgebra(t(solve(SDs)) %*% S %*% solve(SDs), name = 'Sst'),
    mxAlgebra(solve(SDs) %*% A %*% SDs, name = 'Ast')
  )

  params <- omxLocateParameters(model)
  params <- params[params$matrix == 'A' |
                     params$matrix == 'S' & params$row >= params$col, ]
  params_loc <- sprintf('%s[%d, %d]', params$matrix, params$row, params$col)
  st_params_loc <- sprintf('%sst[%d, %d]', params$matrix, params$row, params$col)

  params_alg <- paste0('cbind(rbind(',
                       paste(params_loc,collapse = ', '
                       ), '), rbind(',
                       paste(st_params_loc,collapse = ', '
                       ), '))')

  param_dimnames <- list(paste(vars[params$col],
                               sapply(
                                 params$matrix,
                                 switch,
                                 A = '-->',
                                 S = '<->'
                               ),
                               vars[params$row]),
                         c('raw estimate', 'st estimate'))

  model <- mxModel(model,
                   eval(parse(text = paste0('mxAlgebra(',
                                            params_alg,
                                            ', dimnames = param_dimnames,
                                            name = "Parameter_table")'))),
                   mxAlgebra(t(M), dimnames = list(vars, 'Mean'), name = 'Means'))
  output_tables(model) <- c('Parameter_table', 'Means')

  return(model)
}
