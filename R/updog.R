#' \code{updog} Flexible Genotyping for Autopolyploids
#' 
#' Implements empirical Bayes approaches to genotype
#' autopolyploids from next generation sequencing data while
#' accounting for allelic bias, overdispersion, and sequencing
#' error. The main function is \code{\link{flexdog}}, which 
#' allows the specification
#' of many different genotype distributions. An experimental
#' function that takes into account varying leves of relatedness
#' is implemented in \code{\link{mupdog}}. 
#'
#' @section \code{updog} Functions:
#' \describe{
#'   \item{\code{\link{flexdog}}}{The main function that
#'       fits an empirical Bayes approach to genotype autopolyploids
#'       from next generation sequencing data.}
#'   \item{\code{\link{mupdog}}}{An experimental approach to genotype
#'       autopolyploids that accounts for varying levels of
#'       relatedness between the individuals in the sample.}
#'   \item{\code{\link{rgeno}}}{simulate the genotypes of a sample
#'       from one of the models allowed in \code{\link{flexdog}}.}
#'   \item{\code{\link{rflexdog}}}{Simulate from the 
#'       \code{\link{flexdog}} model.}
#'   \item{\code{\link{plot.flexdog}}}{Plotting the output of
#'       \code{\link{flexdog}}.}
#'   \item{\code{\link{plot.mupdog}}}{Plotting the output of 
#'       \code{\link{mupdog}}.}
#'   \item{\code{\link{summary.mupdog}}}{Providing some summaries 
#'       of the output of \code{\link{mupdog}}.}
#' }
#'
#' @section \code{updog} Datasets:
#' \describe{
#'   \item{\code{\link{snpdat}}}{A small example dataset for using
#'       \code{\link{flexdog}}.}
#'   \item{\code{\link{uitdewilligen}}}{A small example dataset 
#'       for using \code{\link{mupdog}}.}
#'   \item{\code{\link{mupout}}}{The output from fitting 
#'       \code{\link{mupdog}} to \code{\link{uitdewilligen}}.}
#' }
#'
#' @useDynLib updog
#' @importFrom Rcpp sourceCpp
#'
#' @docType package
#' @name updog
#'
#' @author David Gerard
NULL


#' The Beta-Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the beta-binomial distribution when parameterized
#' by the mean \code{mu} and the overdispersion parameter \code{rho}
#' rather than the typical shape parameters.
#'
#'
#' Let \eqn{\mu} and \eqn{\rho} be the mean and overdispersion parameters.
#' Let \eqn{\alpha} and \eqn{\beta} be the usual shape parameters of
#' a beta distribution. Then we have the relation
#' \deqn{\mu = \alpha/(\alpha + \beta),}
#' and
#' \deqn{\rho = 1/(1 + \alpha + \beta).}
#' This necessarily means that
#' \deqn{\alpha = \mu (1 - \rho)/\rho,}
#' and
#' \deqn{\beta = (1 - \mu) (1 - \rho)/\rho.}
#'
#' @param x,q A vector of quantiles.
#' @param p A vector of probabilities.
#' @param n The number of observations.
#' @param size A vector of sizes.
#' @param mu Either a scalar of the mean for each observation,
#'     or a vector of means of each observation, and thus
#'     the same length as \code{x} and \code{size}. This must
#'     be between 0 and 1.
#' @param rho Either a scalar of the overdispersion parameter
#'     for each observation, or a vector of overdispersion
#'     parameters of each observation, and thus the same length as
#'     \code{x} and \code{size}. This must be between 0 and 1.
#' @param log,log_p A logical vector either of length 1 or the same
#'     length as \code{x} and \code{size}. This determines whether
#'     to return the log probabilities for all observations
#'     (in the case that its length is 1) or for
#'     each observation (in the case that
#'     its length is that of \code{x} and \code{size}).
#'
#' @name betabinom
#'
#' @author David Gerard
#'
#'
NULL