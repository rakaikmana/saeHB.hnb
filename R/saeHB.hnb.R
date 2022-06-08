#' saeHB.hnb : Small Area Estimation under Hurdle Negative Binomial Model using Hierarchical Bayesian Method
#'
#' Provides function and datasets for area level of Small Area Estimation under Hurdle Negative Binomial Model using Hierarchical Bayesian (HB) Method. For the reference, see Rao and Molina (2015), Hilbe (2011), and Andika, et al. (2019)
#'
#' @section Author(s):
#' Raka Ikmana, Azka Ubaidillah
#'
#' \strong{Maintainer}: Raka Ikmana \email{221810548@@stis.ac.id}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{HurdleNB}}}{This function gives small area estimator under Hurdle Negative Binomial Model and is implemented to variable of interest (y) that assumed to be a HNB Distribution. The value of variable of interest must be a non-negative data count. This model can be used to handle overdispersion and excess zero in data.}
#' }
#'
#' @section Reference:
#' \itemize{
#'    \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New Jersey: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#'    \item{Hilbe, J. M. (2011). Negative Binomial Regression 2nd Edition. New York : Cambridge University Press. <doi:10.1017/CBO9780511973420>.}
#'    \item{Ntzoufras, I. (2009). Bayesian Modelling Using WinBUGS. New Jersey :  John Wiley & Sons, Inc.}
#'
#' }
#'
#' @docType package
#' @name saeHB.hnb
#'
#' @import stringr
#' @import coda
#' @import rjags
#' @import stats
#' @import grDevices
#' @import graphics
#'
NULL
