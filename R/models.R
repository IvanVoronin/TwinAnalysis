#' @export
cross_lag <- function(data, definition)
{
  vars <- unlist(definition)
  nv <- length(vars)
  model <- mxModel(
    name = 'crosslag',
    type = 'RAM',
    manifestVars = vars,
    # Variances for all variables.
    mxPath(
      from = vars,
      arrows = 2,
      connect = 'single',
      values = 1
    ),

    # Means.
    mxPath(
      from = 'one',
      to = vars,
      connect = 'all.bivariate',
      values = 0
    ),

    mxData(data[, vars], type = 'raw')
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
