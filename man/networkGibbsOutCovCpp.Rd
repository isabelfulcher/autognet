% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{networkGibbsOutCovCpp}
\alias{networkGibbsOutCovCpp}
\title{Run Gibbs sampler for covariates to use in network causal effect estimation}
\usage{
networkGibbsOutCovCpp(
  tau,
  rho,
  nu,
  ncov,
  R,
  N,
  burnin,
  rho_mat,
  adjacency,
  weights,
  cov_mat,
  group_lengths,
  group_functions,
  additional_nu
)
}
\arguments{
\item{tau}{A numeric vector for the intercept terms in the covariate model}

\item{rho}{A numeric vector for the correlation terms in the covariate model}

\item{nu}{A numberic matrix for the neighbor terms in the covariate model}

\item{ncov}{An integer for the number of covariates}

\item{R}{An integer indicating the number of iterations for the Gibbs}

\item{N}{An integer indicating the size of the interconnected network}

\item{burnin}{An integer indicating when to start saving values in the chain}

\item{rho_mat}{A numeric matrix for rho terms}

\item{adjacency}{A binary matrix indicating connected units}

\item{weights}{A numeric vector indicating the number of neighbors for each node}

\item{cov_mat}{A numeric matrix for starting values for each covariate}

\item{group_lengths}{An integer vector indicating the number of categories for each variable}

\item{group_functions}{An integer vector indicating the type of variable}

\item{additional_nu}{An integer (0/1) specifying whether neighbor cross terms will be evaluated (i.e. non-zero)}
}
\value{
A list of numeric matrices that contain the covariate values and neighbor covariate
values for each person at that specific point in the chain
}
\description{
Given the specific inputs, determine covariate
values using a Gibbs sampling procedure.
}
