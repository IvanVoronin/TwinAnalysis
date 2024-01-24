#' @export
cross_lag <- function(data, definition, data_type = 'raw')
{
  vars <- unlist(definition)
  nv <- length(vars)

  if (data_type == 'raw' || length(data$means) > 0) {
    model_means <- mxPath(
      from = 'one',
      to = vars,
      connect = 'all.bivariate',
      values = 0)
  } else {
    model_means <- NULL
  }

#   if (data_type == 'raw') {
#     data <- data[, vars]
#   } else {
# #    data <- list(c(data, type = data_type))
#   }

  model_data <- process_data(
    data,
    vars = vars,
    data_type = data_type
  )

  model <- mxModel(
    name = 'crosslag',
    type = 'RAM',
    manifestVars = vars,
    # Variances for all variables.
    mxPath(
      from = vars,
      arrows = 2,
      connect = 'single',
      values = 1,
      lbound = 0
    ),

    # Means.
    model_means,

    # Data
    model_data
  )

  # Covariances at each time point.
  for (i in definition)
    if (length(i) > 1)
      model <- mxModel(model,
                        mxPath(
                          from = i,
                          arrows = 2,
                          connect = 'unique.bivariate',
                          values = .5
                        ))

  # One way paths.
  for (i in 1:(length(definition) - 1))
    model <- mxModel(model,
                      mxPath(
                        from = definition[[i]],
                        to = definition[[i + 1]],
                        connect = 'all.pairs',
                        values = 0.1
                      ))

  model <- def_stand_params(model)

  return(model)
}
