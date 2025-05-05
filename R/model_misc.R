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
get_output_tables <- function(model, tables = character(0), CIs = TRUE) {
  if (length(tables) == 0)
    tables <- output_tables(model)

  pars <- lapply(tables, mxEvalByName, model)
  names(pars) <- tables
  cis <- if (CIs)
    get_ci(model)
  else NULL

  outp <- list()
  for (i in tables) {
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
#' @title Define saturated and independence models
#'
#' @description A function to define and run two reference models: Saturated and Independence.
#' The reference models are stored in the \code{options} slot of the target model.

#' @param model \code{MxModel} object
#' @param ref_models a named list of length two of \code{MxModels} (Saturated and Independence).
#'        If \code{NULL}, \code{mxRefModels} is used to define reference models.
#' @param run should we fit the reference models? Default is \code{TRUE}.
#' @param ... arguments passed to \code{mxRefModels}.
#'
#' @return the original model with two reference models appended in the \code{option} slot.
#'
#' @note
#' If you provide the list with two reference models, the models have to be stored as \code{'Saturated'}
#' and \code{'Independence'} object within this list. Otherwise they will be thrown away.
#'
#' The reference models summarize the data, so don't forget to update the reference models
#' if you modify the data.
#'
#' @export
ref_models <- function(model, ref_models = NULL, run = TRUE, ...) {
  if (length(ref_models) == 0)
    ref_models <- mxRefModels(model, ...)

  ref_models <- ref_models[c('Saturated', 'Independence')]

  if (run)
    ref_models <- lapply(ref_models, mxRun)

  model@options$ref_models <- ref_models
  return(model)
}

#' @rdname ref_models
#' @export
get_ref_models <- function(model) {
  model@options$ref_models
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
  refmodels <- get_ref_models(model)
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
  out_list <- list()
  for (i in names(dots)) {
    if (run)
      out_list[[i]] <-
        mxRun(
          mxRename(
            mxModel(model, dots[[i]]),
            i))
    else
      out_list[[i]] <-
        mxRename(
          mxModel(model, dots[[i]]),
          i)
  }

  out_list
  # if (run)
  #   lapply(dots, function(x) mxRun(mxRename(mxModel(model, x)), )
  # else
  #   lapply(dots, function(x) mxModel(model, x))
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
    mxAlgebra(solve(SDs) %*% S %*% solve(SDs), name = 'Sst'),
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
                                            name = "Parameter_table")'))))

  if (length(model$M) > 0) {
    model <- mxModel(model,
                     mxAlgebra(t(M), dimnames = list(vars, 'Mean'),
                               name = 'Means'))
    output_tables(model) <- c('Parameter_table', 'Means')
  } else {
    output_tables(model) <- c('Parameter_table')
  }



  return(model)
}

#' @export
process_data <- function(
    data,
    vars,
    data_type = 'raw'
) {
  if (data_type == 'raw') {
    return(
      mxData(
        observed = data,
        type = 'raw'
      )
    )
  } else {
    data$type <- data_type
    data$observed <- data$observed[vars, vars]

    if (length(data$means == 0) || any(is.na(data$means))) {
      data$means <- rep(0, length(vars))
    }
    return(do.call('mxData', data))
  }
}

#' @title Assign the labels to free parameters in MxModel
#'
#' @description
#'
#' @param M \code{MxModel} object
#'
#' @return \code{MxModel} object.
#'
#' @note
#' The labels are assigned in the form "matrix_row_col"
#' or "model_matrix_row_col". In symmetric matrices, the labels
#' are the same upper and lower triangle.
#'
#' @export
assign_auto_labels <- function(M) {
  P <- omxLocateParameters(M)
  one_model <- length(unique(P$model) == 1)
  for (p in split(P, 1:nrow(P))) {
    M[[p$model]][[p$matrix]]$labels[p$row, p$col] <-
      ifelse(
        one_model,
        paste(p$matrix, p$row, p$col, sep = '_'),
        paste(p$model, p$matrix, p$row, p$col, sep = '_')
      )
  }
  for (m in unique(P$matrix)) {
    if (class(M[[m]]) == 'SymmMatrix') {
      m_labs <- M[[m]]$labels
      M[[m]]$labels[upper.tri(m_labs)] <-
        t(m_labs)[upper.tri(m_labs)]
    }
  }
  return(M)
}

#' @title Extract the bootstrapped confidence interval
#'
#' @description
#' If the bootstrap has been applied to the model, this function helps to
#' extract the bootstrapped confidence interval and present it in a table
#' combined with the parameters' estimates. Since the extraction of CI for the
#' computed parameters can take long time, it also allows saving the tables
#' to cache. If the cache file already exists, it will load the tables from it
#'
#' @param model_boot bootstrapped \code{MxModel}
#' @param table_name a character vector of the names of \code{MxAlgebra} or
#'                   \code{MxMatrix} in the model for which the CI is estimated.
#'                   If this is \code{character(0)}, it will try to extract the
#'                   names of the output tables from the model using
#'                   \code{output_tables()} (a character vector stored in
#'                   \code{@output$output_tables})
#' @param interval the interval's quantile (default is 95% CI)
#' @param cache_file a string indicating the file where the parameter tables
#'                   will be stored (.rds)
#' @param force force the computation regardless of whether the cache file
#'              exists or not
#'
#' @return a table in the form of \code{mlth.data.frame}
#'
#'
#' @export
get_boot_ci <- function(
    model_boot,
    table_name = character(0),
    interval = c(0.025, 0.975),
    cache_file = character(0),
    force = FALSE
) {
  require(mlth.data.frame)

  if (length(table_name) == 0)
    table_name <- output_tables(model_boot)

  if (length(table_name) == 0)
    stop('Please provide table_name')

  if (length(cache_file) != 0 && file.exists(cache_file) && !force) {
    ## If cache file exists, read from cache
    outp <- readRDS(cache_file)
    message(
      'Bootstrapped CI loaded from ',
      cache_file
    )

    M1 <- attr(outp, 'original_model')()
    M1@output <- list()
    M2 <- model_boot
    M2@output <- list()

    if (!isTRUE(all.equal(M1, M2)))
      stop(
        'The provided model is different from the original model\n',
        'Restart with force = TRUE'
      )
    outp <- outp[table_name]
  } else {
    ## Check that the folder for the cache file exists before running
    ## the computation
    if (length(cache_file) != 0) {
      cache_folder <- gsub('\\/[^\\/]*$', '', cache_file)
      if (!dir.exists(cache_folder))
        stop(
          'The cache folder ',
          cache_folder,
          ' does not exist'
        )
    }

    ## Extract the CIs from the bootstrapped model
    outp <- list()

    for (i in table_name) {
      message('Processing ', i)
      est <- mxEvalByName(
        i, model_boot
      )
      ci <- mxBootstrapEvalByName(
        i, model_boot, bq = interval
      )

      tab <- vector('list', ncol(est))
      if (length(colnames(est)) !=0 ) {
        names(tab) <- colnames(est)
      } else {
        names(tab) <- paste0('[,', 1:length(tab), ']')
      }
      for (j in 1:ncol(est)) {
        tab[[j]] <- c(
          list(est = est[, j, drop = FALSE]),
          apply(
            ci[(nrow(est) * (j - 1) + 1):(nrow(est) * j), -1, drop = FALSE],
            2, identity,
            simplify = FALSE
          )
        )
      }
      outp[[i]] <- as.mlth.data.frame(tab)
      if (length(rownames(est) != 0)) {
        row.names(outp[[i]]) <- rownames(est)
      } else {
        row.names(outp[[i]]) <- paste0('[', nrow(est), ',]')
      }
    }

    attr(outp, 'original_model') <- function() model_boot

    ## Write the bootstrapped intervals to cache
    if (length(cache_file) != 0) {
      saveRDS(
        outp,
        file = cache_file
      )
      message(
        'The tables with bootstrapped CIs are saved in ',
        cache_file
      )
    }
  }

  return(outp)
}

#' @title Run bootstrap on an MxModel
#'
#' @description
#' It allows running boostrap and saving it to a cache file. If the cache file
#' already exists, it will load the bootstrapped model from cache. It can run
#' in chunks giving an update when each chunk is finished. It measures the time
#' needed to process one chunk and gives an estimation of the total execution
#' time.
#'
#' @param model MxModel
#' @param cache_file a character vector indicating the cache file (.rds)
#' @param N number of bootstrap samples
#' @param force force bootstrap regardless of whether the cache file exists
#' @param by_chunk run the bootstrap in chunks
#' @param chunk_size the size of one chunk
#'
#' @return MxModel
#'
#' @note
#' The bootstrap is only supported for the models with the real (raw) data.
#'
#' There may be different reasons why an the cached model might not match the
#' provided model. From my experience: 1) The parameter estimates may be very
#' unstable, meaning that each run results in very different estimates. There is
#' no good reason for caching or even using such a model. 2) The whole data
#' table provided to the model is stored in the model, including the variables
#' that may not be actually used. Therefore, when the variables are added or
#' changed, the model will not match, even if it does not use these
#' variables. So it is recommended to feed to the model only those variables
#' that are actually used.
#'
#' @export
#'
#'

run_bootstrap <- function(
    ## FIXME: when loading from the cache, check that this is the same model
    ## Same for the bootstrapped CIs
    ## FIXME: It depends on glue for no good reason
  model,
  cache_file = character(0),
  N = 500,
  force = FALSE,
  by_chunk = TRUE,
  chunk_size = 100
) {
  ## If the cached bootstrapped file exists, load from cache
  if (length(cache_file) != 0 && file.exists(cache_file) && !force) {
    model_boot <- readRDS(cache_file)
    message(
      'Bootstrapped model loaded with ',
      model_boot@compute$replications,
      ' replications loaded from ',
      cache_file
    )

    M1 <- attr(model_boot, 'original_model')()
    M1@output <- list()
    M2 <- model
    M2@output <- list()

    if (!all(all.equal(M1, M2) == TRUE))
      stop(
        'The provided model is different from the original model\n',
        'Restart with force = TRUE\n',
        'Details:\n',
        all.equal(M1, M2)
      )
  } else {
    ## Check that the folder for the cache file exists berofe running
    ## the bootstrap
    if (length(cache_file) != 0) {
      cache_folder <- gsub('\\/[^\\/]*$', '', cache_file)
      if (!dir.exists(cache_folder))
        stop(
          'The cache folder ',
          cache_folder,
          ' does not exist'
        )
    }

    model_boot <- model

    if (by_chunk) {
      ## Boostrap run in chunks
      ## The first chunk is timed to estimate the total execution time
      message(
        'Running bootstrap with ',
        N,
        ' replications in the chunks of ',
        chunk_size
      )

      ## Start the bootstrap
      for (i in 1:ceiling(N / chunk_size)) {
        if (i == 1)
          then <- proc.time()[3]

        model_boot <- mxBootstrap(
          model_boot,
          replications = min(i * chunk_size, N)
        )

        message(
          model_boot@compute$replications,
          ' replications completed'
        )
        if (i == 1 & N > chunk_size) {
          now <- proc.time()[3]
          t <- floor(now - then)
          t_full <- floor(t * N / chunk_size)
          message(
            'The first chunk was completed in ',
            paste0(t %/% 3600, 'h'),
            paste0(t %/% 60, 'm'),
            paste0(t %% 60, 's'),
            '\n',
            'The estimated total execution time is ',
            paste0(t_full %/% 3600, 'h'),
            paste0(t_full %/% 60, 'm'),
            paste0(t_full %% 60, 's')
          )
        }
      }
    } else {
      ## Bootstrap in a single run
      message(
        'Running bootstrap with ',
        N,
        ' replications'
      )

      ## Start the boostrap
      model_boot <- mxBootstrap(
        model_boot,
        replications = N
      )
    }

    attr(model_boot, 'original_model') <- function() model

    ## Write the bootstrapped model to cache
    if (length(cache_file) != 0) {
      saveRDS(
        model_boot,
        file = cache_file
      )
      message(
        'The bootstrapped model is saved in ',
        cache_file
      )
    }
  }

  return(model_boot)
}
