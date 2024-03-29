% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_linearization_functions.R
\name{hessian_linear_pred}
\alias{hessian_linear_pred}
\title{Hessian matrix with respect to linear predictors.}
\usage{
hessian_linear_pred(allocation_prob, y, linear_pred, n_comp)
}
\arguments{
\item{allocation_prob}{A vector of component allocation probabilities.}

\item{y}{A vector of response variable values.}

\item{linear_pred}{A vector of the linear predictors in both components and mixture weights models.}

\item{n_comp}{Number of mixture components.}
}
\value{
The hessian matrix with respect to the linear predictors.
}
\description{
Computing the posterior gradients with respect to the linear predictors.
}
