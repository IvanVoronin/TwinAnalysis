#' @rdname univ_descriptives
#' @title Compute the univariate twin descriptives
#' @description
#' Compute number of valid observations, means, SDs and cross-twin correlations
#' for MZ and DZ twins. Can be used for any data on groups.
#'
#' @param dat the data.frame with the variables zyg and paste(vars, rep(1:2, each = nv))
#' @param zyg the name of the zygosity variable (character string)
#' @param vars the names of the variables to compute the descriptives for.
#'             If vars = c('X', 'Y', 'Z'), the function looks for
#'             c('X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2')
#' @details
#'
#' @examples
#'
#' @export
get_univ_descriptives <- function(dat, zyg, vars) {
  nv <- length(vars)
  selvars <- paste0(vars, rep(1:2, each = nv))
  zyg <- factor(dat[, zyg])

  descr <- lapply(vars, function(y) {
    outp <- do.call('rbind', by(dat, zyg, function(x) {
      tw1 <- x[, paste0(y, 1)]
      tw2 <- x[, paste0(y, 2)]
      ct <- cor.test(tw1, tw2)
      mlth.data.frame(Twin1 = data.frame(n = sum(!is.na(tw1)),
                                         M = mean(tw1, na.rm = TRUE),
                                         SD = sd(tw1, na.rm=TRUE)),
                      Twin2 = data.frame(n = sum(!is.na(tw2)),
                                         M = mean(tw2, na.rm = TRUE),
                                         SD = sd(tw2, na.rm=TRUE)),
                      Cor = data.frame(r = ct$estimate,
                                       lowerCI = ct$conf.int[1],
                                       upperCI = ct$conf.int[2]))
      }))

    row.names(outp) <- levels(zyg)
    return(outp)
    })
  names(descr) <- vars
  return(descr)
}


#' @rdname cross_trait_cors
#' @title Compute the cross-twin cross-trait correlations
#' @description
#' Compute the phenotypic correlations across twins and across traits,
#' for MZ and DZ twins. Can be used for any data on groups.
#'
#' @param dat the data.frame with the variables zyg and paste(vars, rep(1:2,
#'            each = nv))
#' @param zyg the name of the zygosity variable (character string)
#' @param vars the names of the variables to compute the descriptives for.
#'             If vars = c('X', 'Y', 'Z'), the function looks for
#'             c('X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2')
#' @param compact return the MZ and DZ correlations in one table, MZ
#'                correlations on the upper triangle and DZ correlations on
#'                the lower triangle. Works only when the zygosity variable
#'                has exactly two levels
#'
#' @examples
#'
#' @export
get_cross_trait_cors <- function(dat, zyg, vars, compact = FALSE) {
  nv <- length(vars)
  selvars <- paste0(vars, rep(1:2, each = nv))
  zyg <- factor(dat[, zyg])

  cors <- by(dat[, selvars], zyg, function(x) {
    cor(x, use = 'pairwise')
  })

  zyglevels <- levels(zyg)
  if (compact && length(zyglevels) == 2) {
    cors1 <- cors[[1]]
    cors2 <- cors[[2]]
    cors1[row(cors1) > col(cors1)] <- cors2[row(cors2) > col(cors2)]
    diag(cors1) <- NA

    cors<-as.data.frame(cors1)
    cors<-cbind(row.names(cors), cors)
    names(cors)[1] <- paste(rev(zyglevels), collapse = ' \\ ')
    row.names(cors) <- NULL
  }

  return(cors)
}

#' @rdname fisher_z
#' @title Perform the Fisher test for two correlations
#'
#' @description
#' The Fisher r-to-z transformation implied to compare the correlations
#' in two groups
#'
#' @param dat the data.frame with the variables zyg and paste(vars, rep(1:2,
#'            each = nv))
#' @param zyg the name of the zygosity variable (character string),
#'            must have exactly two levels
#' @param vars the names of the variables to perform the test.
#'             If vars = c('X', 'Y', 'Z'), the function looks for
#'             c('X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2')
#' @param compare the labels of groups to compare, must be either character(0)
#'                or length of 2
#'
#' @examples
#'
#' @export
get_fisher_z <- function(dat, zyg, vars, compare = character(0)) {
  zyg <- factor(dat[, zyg])
  zyglevels <- levels(zyg)

  if (length(compare) == 0)
    compare <- zyglevels

  if (length(compare) != 2)
    stop('compare must have the length of 0 or 2\n',
         'when compare has the length of 0, zyg myst have exactly two levels')

  if (any(!compare %in% zyglevels))
    stop('compare must be a subset of zyg levels')

  fisher_z <- function(r1, r2, n1, n2)
    (atanh(r1) - atanh(r2)) / ((1 / (n1 - 3)) + (1 / (n2 - 3)))^0.5

  outp <- do.call('rbind', lapply(vars, function(x) {
    cor1 <- cor.test(dat[zyg == compare[1], paste0(x, 1)],
                     dat[zyg == compare[1], paste0(x, 2)])
    cor2 <- cor.test(dat[zyg == compare[2], paste0(x, 1)],
                     dat[zyg == compare[2], paste0(x, 2)])
    r1 <- cor1$estimate
    r2 <- cor2$estimate
    n1 <- cor1$parameter + 2
    n2 <- cor2$parameter + 2

    z <- fisher_z(r1, r2, n1, n2)
    p <- 2 * (1 - pnorm(abs(z)))

    data.frame(r1, n1, r2, n2, z, p)
  }))

  dimnames(outp) <- list(vars,
                          c(paste(c('r', 'n'),
                                  rep(compare, each = 2),
                                  sep = '_'),
                            'z', 'p'))
  return(outp)
}
