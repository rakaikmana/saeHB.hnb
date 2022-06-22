#' @title Small Area Estimation using Hierarchical Bayesian under Hurdle Negative Binomial Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Hurdle Negative Binomial Distribution. The value of variable of interest must be a non-negative data count. This model can be used to handle overdispersion and excess zero in data.
#' @param formula Formula that describe the fitted model
#' @param iter.update Number of updates with default \code{3}
#' @param coef.nonzero Optional vector for the mean of the prior distribution of the model coefficients \eqn{(\beta)} for variable of interest \eqn{(y)} which value is zero count
#' @param var.coef.nonzero Optional vector of variance of coefficient non-zero count
#' @param coef.zero Optional vector for the mean of the prior distribution of the model coefficients \eqn{(\gamma)} for variable of interest \eqn{(y)} which value is zero count
#' @param var.coef.zero Optional vector for variance of coefficient zero count
#' @param iter.mcmc Number of total iterations per chain with default \code{2000}
#' @param thin Thinning rate, must be a positive integer with default \code{1}
#' @param burn.in Number of iterations to discard at the beginning with default \code{1000}
#' @param tau.u Variance of random effect area for non-zero count of variable interest with default \code{1}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{alpha}{Dispersion parameter}
#'    \item{plot}{Trace, Density, Autocorrelation Function Plot of MCMC samples}
#'
#' @examples
#' ##For data without any non-sampled area
#' data(dataHNB)     # Load dataset
#'
#' \donttest{result <- HurdleNB(y ~ x1 + x2, data = dataHNB)}
#'
#' \donttest{result$Est}          # Small Area mean Estimates
#' \donttest{result$refVar}       # Estimated random effect variances
#' \donttest{result$coefficient}  # Estimated model coefficient
#' \donttest{result$alpha}        # Estimated dispersion parameter
#'
#' # Load library 'coda' to execute the plot
#' # autocorr.plot(result$plot[[3]])    # Generate ACF Plot
#' # plot(result$plot[[3]])             # Generate Density and Trace plot
#'
#' ## For data with non-sampled area use dataHNBNs
#'
#' @import stats
#' @import rjags
#' @import stringr
#' @import grDevices
#' @import graphics
#' @import coda
#'
#' @export

HurdleNB <- function(formula,iter.update=3, iter.mcmc=2000,
                     coef.nonzero, var.coef.nonzero, coef.zero, var.coef.zero,
                     thin = 1, burn.in =1000, tau.u = 1, data){

  result <- list(Est = NA, refVar = NA, coefficient = NA, alpha = NA, plot = NA)

  formuladata <- model.frame(formula,data,na.action=NULL)
  if (any(is.na(formuladata[,-1])))
    stop("Auxiliary Variables contains NA values.")
  auxVar <- as.matrix(formuladata[,-1])
  nvar <- ncol(auxVar) + 1

  if (!missing(var.coef.nonzero)){

    if( length(var.coef.nonzero) != nvar ){
      stop("length of vector var.coef nonzero does not match the number of regression coefficients, the length must be ",nvar)
    }

    tau.b.value = 1/var.coef.nonzero
  } else {
    tau.b.value = 1/rep(1,nvar)
  }

  if (!missing(coef.nonzero)){
    if( length(coef.nonzero) != nvar ){
      stop("length of vector coef nonzero does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.b.value = coef.nonzero
  } else {
    mu.b.value = rep(0,nvar)
  }

  if (!missing(var.coef.zero)){

    if( length(var.coef.zero) != nvar ){
      stop("length of vector var.coef zero does not match the number of regression coefficients, the length must be ",nvar)
    }

    tau.g.value = 1/var.coef.zero
  } else {
    tau.g.value = 1/rep(1,nvar)
  }

  if (!missing(coef.zero)){
    if( length(coef.zero) != nvar ){
      stop("length of vector coef zero does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.g.value = coef.zero
  } else {
    mu.g.value = rep(0,nvar)
  }

  if (iter.update < 3){
    stop("the number of iteration updates at least 3 times")
  }

  #Fungsi Tersampel
  if (!any(is.na(formuladata[,1]))){

    formuladata <- as.matrix(na.omit(formuladata))
    x <- model.matrix(formula,data = as.data.frame(formuladata))
    n <- nrow(formuladata)
    mu.b = mu.b.value
    mu.g = mu.g.value
    tau.b = tau.b.value
    tau.g = tau.g.value
    tau.ua = tau.ub = 1
    a.var.u = 1
    alpha = 1
    tau.aa = tau.ab = 0.001

    for (i in 1:iter.update){
      dat <- list("n"= n,  "nvar"= nvar, "zeros"= rep(0,n), "y" = formuladata[,1], "x"=as.matrix(x[,-1]),
                  "mu.b"=mu.b,"mu.g"=mu.g, "tau.b"=tau.b, "tau.g"=tau.g,
                  "tau.ua"=tau.ua, "tau.ub"=tau.ub, "tau.aa"=tau.aa, "tau.ab"=tau.ab)

      inits <- list(u = rep(0, n), b = mu.b, g = mu.g, tau.u = tau.u, alpha=alpha)

      cat("model{
       #Likelihood using zero trick
          for (i in 1:n){
              zeros[i] ~ dpois(zeros.mean[i])
              zeros.mean[i] <- -ll[i] + 1000000
              B[i] ~ dbern(phi[i])

              logztnb[i] <- loggam(y[i]+alpha)
                            - loggam(y[i]+1)
                            - loggam(alpha)
                            + y[i]*(log(mu[i])-log(mu[i]+alpha))
                            + alpha*(log(alpha)-log(mu[i]+alpha))
                            + log((1-(alpha/(mu[i]+alpha))^alpha)^(-1))


              z[i]  <- step(y[i] - 0.0001)
              l1[i] <- (1 - z[i])*log(phi[i])
              l2[i] <- z[i]*(log(1-phi[i]) + logztnb[i])
              ll[i] <- l1[i] + l2[i]

              #Model regresi
              log(mu[i])   <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
              logit(phi[i]) <- g[1] + sum(g[2:nvar]*x[i,])

              mu.exp[i] <- mu[i]*(1-phi[i])*(1-(alpha/(mu[i]+alpha))^alpha)^(-1)

              #Random effect area
              u[i] ~ dnorm(0,tau.u)

          }

          #Priors
          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					    g[k] ~ dnorm(mu.g[k],tau.g[k])
          }

					alpha ~ dgamma(tau.aa, tau.ab)

          tau.u ~ dgamma(tau.ua, tau.ub)

          a.var.u <- 1/tau.u


        }", file="hnb.txt")

      jags.m <- jags.model( file = "hnb.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("hnb.txt")
      params <- c("mu.exp","a.var.u","alpha","b","g","tau.u")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps=summary(samps1)

      a.var.u = result_samps$statistics[1]
      alpha   = result_samps$statistics[2]

      beta    = result_samps$statistics[3:(nvar+2),1:2]
      gamma   = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]

      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
        mu.g[i]  = gamma[i,1]
        tau.b[i] = 1/(gamma[i,2]^2)
      }

      tau.aa  = result_samps$statistics[2,1]^2/result_samps$statistics[2,2]^2
      tau.ab  = result_samps$statistics[2,1]/result_samps$statistics[2,2]^2


      tau.ua  = result_samps$statistics[2*nvar+n+3,1]^2/result_samps$statistics[2*nvar+n+3,2]^2
      tau.ub  = result_samps$statistics[2*nvar+n+3,1]/result_samps$statistics[2*nvar+n+3,2]^2


    }
    result_samps=summary(samps1)

    b.varnames <- list()
    g.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <- str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
      idx.g.varnames <- as.character(i-1)
      g.varnames[i] <- str_replace_all(paste("g[",idx.g.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc <- samps1[,c(3:(nvar+2))]
    colnames(result_mcmc[[1]]) <- b.varnames

    result_mcmc1 <- samps1[,c((nvar+3):(2*nvar+2))]
    colnames(result_mcmc1[[1]]) <- g.varnames

    a.var.u = result_samps$statistics[1]

    a.var   = a.var.u

    alpha = result_samps$statistics[2]

    beta = result_samps$statistics[3:(nvar+2),1:2]
    rownames(beta) <- b.varnames

    gamma = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]
    rownames(gamma) <- g.varnames

    mu = result_samps$statistics[(2*nvar+3):(2*nvar+n+2),1:2]

    Estimation = data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles)
    q_beta    <- Quantiles[3:(nvar+2),]
    q_gamma   <- Quantiles[(nvar+3):(2*nvar+2),]
    q_mu      <- Quantiles[(2*nvar+3):(2*nvar+n+2),]

    rownames(q_beta)  <- b.varnames
    rownames(q_gamma) <- g.varnames
    beta  <- cbind(beta,q_beta)
    gamma <- cbind(gamma,q_gamma)
    coef  <- rbind(beta,gamma)

    Estimation <- data.frame(Estimation,q_mu)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")

  } else {
    #NONSAMPLED
    formuladata <- as.data.frame(formuladata)

    x <- as.matrix(formuladata[,2:nvar])
    n <- nrow(formuladata)

    mu.b = mu.b.value
    mu.g = mu.g.value
    tau.b = tau.b.value
    tau.g = tau.g.value
    tau.ua = tau.ub = 1
    alpha = 1
    a.var.u = 1
    tau.aa = tau.ab = 0.001


    formuladata$idx <- rep(1:n)
    data_sampled    <- na.omit(formuladata)
    data_nonsampled <- formuladata[-data_sampled$idx,]

    r  = data_nonsampled$idx
    n1 = nrow(data_sampled)
    n2 = nrow(data_nonsampled)

    for (i in 1:iter.update){
      dat <- list("n1"= n1, "n2"=n2,"nvar"=nvar,"zeros_sampled"=rep(0,n1),
                  "y_sampled" = data_sampled[,1],
                  "x_sampled"=as.matrix(data_sampled[,2:nvar]),
                  "x_nonsampled"=as.matrix(data_nonsampled[,2:nvar]),
                  "mu.b"=mu.b,"mu.g"=mu.g, "tau.b"=tau.b, "tau.g"=tau.g,
                  "tau.ua"=tau.ua, "tau.ub"=tau.ub,
                  "tau.aa"=tau.aa,"tau.ab"=tau.ab)

      inits <- list(u = rep(0,n1), uT = rep(0,n2), b = mu.b, g = mu.g,
                    tau.u = tau.u, alpha=alpha)

      cat("model{
       #Likelihood using zero trick
          for (i in 1:n1){
              zeros_sampled[i] ~ dpois(zeros.mean[i])
              zeros.mean[i] <- -ll[i] + 100000000
              B[i] ~ dbern(phi[i])

              logztnb[i] <- loggam(y_sampled[i]+alpha)
                            - loggam(y_sampled[i]+1)
                            - loggam(alpha)
                            + y_sampled[i]*(log(mu[i])-log(mu[i]+alpha))
                            + alpha*(log(alpha)-log(mu[i]+alpha))
                            + log((1-(alpha/(mu[i]+alpha))^alpha)^(-1))


              z[i]  <- step(y_sampled[i] - 0.0001)
              l1[i] <- (1 - z[i])*log(phi[i])
              l2[i] <- z[i]*(log(1-phi[i]) + logztnb[i])
              ll[i] <- l1[i] + l2[i]

              #Model regresi
              log(mu[i])   <- b[1] + sum(b[2:nvar]*x_sampled[i,]) + u[i]
              logit(phi[i]) <- g[1] + sum(g[2:nvar]*x_sampled[i,])

              mu.exp[i] <- mu[i]*(1-phi[i])*(1-(alpha/(mu[i]+alpha))^alpha)^(-1)

              #Random effect area
              u[i] ~ dnorm(0,tau.u)


          }

          for (j in 1:n2){

              #Model regresi

              log(muT[j])    <- b[1] + sum(mu.b[2:nvar]*x_nonsampled[j,]) + uT[j]
              logit(phiT[j]) <- g[1] + sum(mu.g[2:nvar]*x_nonsampled[j,])

              mu.nonsampled[j] <- muT[j]*(1-phiT[j])*(1-(alpha/(muT[j]+alpha))^alpha)^(-1)

              #Random effect area
              uT[j] ~ dnorm(0,tau.u)


          }

          #Priors
          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					    g[k] ~ dnorm(mu.g[k],tau.g[k])
          }

					alpha ~ dgamma(tau.aa,tau.ab)

          tau.u ~ dgamma(tau.ua, tau.ub)

          a.var.u <- 1/tau.u



        }", file="hnb.txt")

      jags.m <- jags.model( file = "hnb.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("hnb.txt")
      params <- c("mu.exp","mu.nonsampled","a.var.u","alpha","b","g","tau.u")
      samps  <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps = summary(samps1)

      a.var.u = result_samps$statistics[1]
      alpha   = result_samps$statistics[2]

      beta    = result_samps$statistics[3:(nvar+2),1:2]
      gamma   = result_samps$statistics[nvar+3:(2*nvar+2),1:2]

      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
        mu.g[i]  = gamma[i,1]
        tau.b[i] = 1/(gamma[i,2]^2)
      }

      tau.aa  = result_samps$statistics[2,1]^2/result_samps$statistics[2,2]^2
      tau.ab  = result_samps$statistics[2,1]/result_samps$statistics[2,2]^2



      tau.ua  = result_samps$statistics[2*nvar+n+3,1]^2/result_samps$statistics[2*nvar+n+3,2]^2
      tau.ub  = result_samps$statistics[2*nvar+n+3,1]/result_samps$statistics[2*nvar+n+3,2]^2


    }
    result_samps = summary(samps1)

    b.varnames <- list()
    g.varnames <- list()

    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i]  <- str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
      idx.g.varnames <- as.character(i-1)
      g.varnames[i]  <- str_replace_all(paste("g[",idx.g.varnames,"]"),pattern=" ", replacement="")
    }

    #beta and gamma for plot
    result_mcmc <- samps1[,c(3 :(nvar+2))]
    colnames(result_mcmc[[1]]) <- b.varnames

    result_mcmc1 <- samps1[,c((3+nvar) :(2*nvar+2))]
    colnames(result_mcmc1[[1]]) <- g.varnames

    #random effect area
    a.var.u = result_samps$statistics[1]
    a.var   = a.var.u

    #dispersion params
    alpha = result_samps$statistics[2]

    #regression coef
    beta  = result_samps$statistics[3:(nvar+2),1:2]
    gamma = result_samps$statistics[(nvar+3):(2*nvar+2),1:2]
    rownames(beta)  <- b.varnames
    rownames(gamma) <- g.varnames

    #estimation
    mu            = result_samps$statistics[(2*nvar+3):(2+2*nvar+n1),1:2]
    mu.nonsampled = result_samps$statistics[(3+2*nvar+n1):(2+2*nvar+n),1:2]

    Estimation      = matrix(rep(0,n),n,2)
    Estimation[r,]  = mu.nonsampled
    Estimation[-r,] = mu
    Estimation      = as.data.frame(Estimation)

    #Quantiles
    Quantiles       <- as.data.frame(result_samps$quantiles)
    q_beta          <- (Quantiles[3:(nvar+2),])
    q_gamma         <- (Quantiles[(nvar+3):(2*nvar+2),])
    q_mu            <- (Quantiles[(2*nvar+3):(2*nvar+2+n1),])
    q_mu.nonsampled <- (Quantiles[(3+2*nvar+n1):(2+2*nvar+n),])
    q_Estimation    <- matrix(0,n,5)

    for (i in 1:5){
      q_Estimation[r,i]  <- q_mu.nonsampled[,i]
      q_Estimation[-r,i] <- q_mu[,i]
    }

    rownames(q_beta)  <- b.varnames
    rownames(q_gamma) <- g.varnames
    beta  <- cbind(beta,q_beta)
    gamma <- cbind(gamma,q_gamma)
    coef  <- rbind(beta,gamma)

    Estimation <- data.frame(Estimation,q_Estimation)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
  }


  result$Est                  = Estimation
  result$refVar               = a.var
  result$coefficient          = coef
  result$alpha                = alpha
  result$plot                 = list(graphics.off(), par(mar=c(2,2,2,2)),
                                     autocorr.plot(result_mcmc,col="brown2",lwd=2),
                                     plot(result_mcmc,col="brown2",lwd=2),
                                     autocorr.plot(result_mcmc1,col="brown2",lwd=2),
                                     plot(result_mcmc1,col="brown2",lwd=2))

  return(result)

}
