% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_misc.R
\name{def_nested_models}
\alias{def_nested_models}
\title{Define nested models.}
\usage{
def_nested_models(model, ..., run = FALSE)
}
\arguments{
\item{model}{\code{MxModel} object.}

\item{...}{named \code{OpenMx} objects to add or replace in the target
model (in a list if there are more than one).}

\item{run}{should the models be run? Default is \code{FALSE}.}
}
\value{
The named list of \code{MxModel}s.
}
\description{
Define nested models by adding \code{MxConstraint}s or modifying the labels
or other elements of the target \code{MxModel}.
}
