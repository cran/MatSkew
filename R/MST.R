#' @keywords internal
M_update <- function(X, estep, N) {
  const <- mean(estep[[1]]) * estep[[2]] - 1
  denom<-sum(estep[[2]])*mean(estep[[1]])-N
  Mnew <- Reduce("+", lapply(1:N, function(i) X[[i]] * const[i]))/denom
  return(Mnew)
}
#' @keywords internal
A_update <- function(X, estep, N) {
  const <- mean(estep[[2]]) - estep[[2]]
  denom<-sum(estep[[2]])*mean(estep[[1]])-N
  Anew <- Reduce("+", lapply(1:N, function(i) X[[i]] * const[i]))/denom
}
#' @keywords internal
nu_update <- function(estep) {
  dfnew <- try(uniroot(function(v) log(v/2) + 1 - digamma(v/2) - mean(estep[[2]]) - mean(estep[[3]]),
                       lower = 0.1, upper = 1000)$root,silent=F)
  nu <- dfnew
  if (is.double(dfnew)) {

    if (dfnew > 200) {
      nu <- 200
    }

    if (dfnew < 2) {
      nu <- 2
    }

  } else {
    w.message <- "Difficulty updating degrees of freedom"
    #print(w.message)
    #print(estep)
    prob <- 1
    return(list(w.message, prob=prob))
  }

  return(list(nu=nu,prob=0))
}
#' @keywords internal
Sigma_update <- function(X, M, A, Psi, estep, N, p, tol_sing = sqrt(.Machine$double.eps)) {
  Dev <- lapply(1:N, function(i) (X[[i]] - M))
  Dev_Sum <- Reduce("+", lapply(1:N, function(i) estep[[2]][i] * Dev[[i]] %*% solve(Psi) %*%
                                  t(Dev[[i]])))
  Crossterm1 <- Reduce("+", lapply(1:N, function(i) Dev[[i]] %*% solve(Psi) %*% t(A)))
  Crossterm2 <- Reduce("+", lapply(1:N, function(i) A %*% solve(Psi) %*% t(Dev[[i]])))
  Sigmanew <- (Dev_Sum - Crossterm1 - Crossterm2 + sum(estep[[1]]) * A %*% solve(Psi) %*%
                 t(A))/(N * p)

  detSig <- det(Sigmanew)

  if (detSig < tol_sing) {
    w.message <- "Singular Sigma Matrix"
    return(list(w.message, prob = 1, Sigmanew))
  }
  return(list(Sigmanew, prob = 0))
}
#' @keywords internal
Psi_update <- function(X, M, A, Sigma, estep, N, n, tol_sing = sqrt(.Machine$double.eps)) {
  Dev <- lapply(1:N, function(i) (X[[i]] - M))
  Dev_Sum <- Reduce("+", lapply(1:N, function(i) estep[[2]][i] * t(Dev[[i]]) %*% solve(Sigma) %*%
                                  Dev[[i]]))
  Crossterm1 <- Reduce("+", lapply(1:N, function(i) t(A) %*% solve(Sigma) %*% (Dev[[i]])))
  Crossterm2 <- Reduce("+", lapply(1:N, function(i) t(Dev[[i]]) %*% solve(Sigma) %*% (A)))
  Psinew <- (Dev_Sum - Crossterm1 - Crossterm2 + sum(estep[[1]]) * t(A) %*% solve(Sigma) %*%
               (A))/(N * n)
  #astar<-ncol(Psinewtest)/sum(diag(Psinewtest))
  detPsi <- det(Psinew)

  if (detPsi < tol_sing) {
    w.message <- "Singular Psi Matrix"
    return(list(w.message, prob = 1, Psinew))
  }

  return(list(Psinew, prob = 0))
}
#' @keywords internal
Deltaupdate <- function(X, M, Sigma, Psi, nu) {
  Delta <- lapply(X, function(x) sum(diag((solve(Sigma) %*% (x - M) %*% solve(Psi) %*% t(x -
                                                                                           M)))) + nu)
  Deltat <- unlist(Delta)
  return(Deltat)
}
#' @keywords internal
rhoupdate <- function(A, Sigma, Psi, N) {
  rho <- sum(diag(solve(Sigma) %*% A %*% solve(Psi) %*% t(A)))
  rhot <- rep(rho, N)
  return(rho)
}

#' @keywords internal
Estep_updates <- function(X, M, A, Sigma, Psi, N, nu, n, p) {
  lambda <- -((nu + n * p)/2)
  Delta <- Deltaupdate(X, M, Sigma, Psi, nu)
  rho <- rhoupdate(A, Sigma, Psi, N)
  #print(rho)
  #print(Delta)
  ai<-NULL
  bi<-NULL
  ci<-NULL
  for (i in 1:N){
    ai[i] <- myEgig(lambda, Delta[i], rho, func = "x")
    bi[i] <- myEgig(lambda, Delta[i], rho, func = "1/x")
    ci[i] <- myEgig(lambda, Delta[i], rho, func = "logx")
  }
  if (any(is.na(ci))){
    w.message<-"Problem Calculating log moment"
    #print(w.message)
    return(list(w.message,prob=1))
  }
  return(list(ai, bi, ci,prob=0))
}

#' @keywords internal
MVskewt <- function(X, M, A, Sigma, Psi, nu, n, p, N) {
  lambda = -(nu + n * p)/2
  Delta <- sum(diag((solve(Sigma) %*% (X - M) %*% solve(Psi) %*% t(X - M)))) + nu
  rho <- sum(diag(solve(Sigma) %*% A %*% solve(Psi) %*% t(A)))
  BessX<-besselK(sqrt(Delta*rho),lambda)
  if (BessX==0){
    BessX<-0.0001
  }
  dens <-log(2)+(nu/2)*log(nu/2)-log(2*pi)*(n*p/2)-log(det(Sigma))*(p/2)-log(det(Psi))*(n/2)-log(gamma(nu/2))+(lambda/2)*(log(Delta)-log(rho))+log(BessX)+sum(diag(solve(Sigma)%*%(X-M)%*%solve(Psi)%*%t(A)))

  return(list(dens=dens,lambda=lambda,rho=rho,Delta=Delta))
}

#' Matrix Skew t Parameter Estimation
#'
#' Performs paramter estimation for the matrix variate skew-t distribution using an ECM algorithm.
#'@param X A list of matrices of the same size
#'@param Tol The tolerance of the ECM algorithm. Defaults to 0.001
#'@param max_iter The maximum number of iterations. Defaults to 1000
#'@return Returns a list with elements M (the estimate of the location), A (the estimate of the skewness), nu (the estimate of the degrees of freedom), Sigma (the estimate of Sigma), Psi (the estimate of Psi), loglik (a vector of log likelihood values), flag (returns TRUE if a numerical issue occured, FALSE otherwise).
#'@export
#'@examples
#'data(SimX)
#'Fit_st<-Fit_Skewt(SimX)
Fit_Skewt <- function(X, Tol = 0.001, max_iter = 1000) {
  n1 <- dim(X[[1]])[1]
  p <- dim(X[[1]])[2]
  N <- length(X)
  M_in<-Reduce("+",X)/N
  A_in<-matrix(0.1,n1,p)
  Sigma_in<-diag(1,n1,n1)
  Psi_in<-diag(1,p,p)
  nu_in<-50

  #Initialize ai, bi and ci

  estep <- Estep_updates(X, M_in, A_in, Sigma_in, Psi_in, N, nu_in, n1, p)
  if (estep$prob==1){
    w.message<-"Problem with log moment on initialization"
    #print(w.message)
    return(w.message=w.message,flag=T)

  }

  Sigma <- Sigma_in
  Psi <- Psi_in

  conv <- 0
  prob <- 0
  iter <- 1
  likCM1 <- NULL
  likCM2 <- NULL
  likCM3 <- NULL

  while (conv == 0) {
    #Update Parameters
    M <- M_update(X, estep, N)
    A <- A_update(X, estep, N)
    dfnew <- nu_update(estep)

    if (dfnew$prob==0) {
      nu <- dfnew$nu
    } else {
      w.message<-"Problem with updating the degrees of freedom."
      return(list(w.message=w.message,flag=TRUE))
    }

    Sigmatest <- Sigma_update(X, M, A, Psi, estep, N, p)
    if (Sigmatest$prob == 0) {
      Sigma <- Sigmatest[[1]]
    } else {
      w.message<-"Problem with Sigma Update."
      return(list(w.message=w.message,flag=TRUE))
    }

    Psitest <- Psi_update(X, M, A, Sigma, estep, N, n1)
    if (Psitest$prob == 0) {
      Psi <- Psitest[[1]]
    } else {
      return(list(Psitest, iter,flag=TRUE))
    }

    #Update E step
    estep <- Estep_updates(X, M, A, Sigma, Psi, N, nu, n1, p)
    if (estep$prob==1){
      w.message<-"Problem with the E-Step."
      return(list(w.message=w.message,flag=TRUE))
    }


    #Calculate Likelihood
    Dens <- unlist(lapply(X, function(x) MVskewt(x, M, A, Sigma, Psi, nu, n1, p, N)$dens))

    likCM3[iter] <- sum(Dens)
    if (is.na(likCM3[iter])){
      w.message<-"NA in likelihood!"
      #print(w.message)
      return(list(w.message=w.message,flag=TRUE))
    }

    if (is.infinite(likCM3[iter])){
      w.message<-"Infinite likelihood!"
      #print(w.message)
      #print(Dens)

      return(list(w.message=w.message,flag=TRUE))
    }


    if (iter >1){
      if ((likCM3[iter] - likCM3[iter-1]) < 0) {
        w.message <- "Decreasing Likelihood!"
        #print(w.message)
        return(list(w.message,flag=TRUE))
      }
    }

    if (iter > 3) {
      if ((likCM3[iter - 1] - likCM3[iter - 2]) == 0) {
        conv <- 1
      } else {
        ak <- (likCM3[iter] - likCM3[iter - 1])/(likCM3[iter - 1] - likCM3[iter - 2])
        linf <- likCM3[iter - 1] + (likCM3[iter] - likCM3[iter - 1])/(1 - ak)
        if (abs(linf - likCM3[iter - 1]) < Tol) {
          conv <- 1
          #print(iter)
        }
      }
    }
    iter <- iter + 1
    if (iter > max_iter) {
      w.message <- paste("Did Not Converge after ", iter - 1, " iterations. Consider increasing max_iter or Tol",
                         sep = "")
      return(list(w.message,flag=TRUE))
    }
  }
  Sigma<-Sigma/Sigma[1,1]
  Psi<-Sigma[1,1]*Psi
  return(list(M=M, A=A, nu=nu, Sigma=Sigma, Psi=Psi,likelihood=likCM3,flag=FALSE))
}

#' @keywords internal
myEgig<-function (lambda, chi, psi, func = c("x", "logx", "1/x"))
{

  if (func == "x") {

    alpha.bar <- sqrt(chi * psi)
    term1 <- 0.5 * log(chi/psi)
    term2 <- mybessel(abs(lambda + 1), alpha.bar)
    term3 <- mybessel(abs(lambda), alpha.bar)
    return(exp(term1 + term2 - term3))
  }
  else if (func == "logx") {
    alpha.bar <- sqrt(chi * psi)

    Kderiv <- besderiv(lambda,alpha.bar)
    BesselX<-mybessel(abs(lambda),alpha.bar)
    res<-0.5 * log(chi/psi) + Kderiv

    return(res)
  }
  else if (func == "1/x") {
    alpha.bar <- sqrt(chi * psi)
    term1 <- -0.5 * log(chi/psi)
    term2 <- mybessel(abs(lambda - 1), alpha.bar)
    term3 <- mybessel(abs(lambda), alpha.bar)
    return(exp(term1 + term2 - term3))
  }
}

#' @keywords internal
mybessel<-function(nu,z){
  bes<-tryCatch(besselK(z,nu),error=function(e) return(NULL))

  if (is.null(bes)) {
    bes2 <- 0.5 * (log(pi) - log(2) - log(nu)) - nu * log(exp(1) * z) + nu *
      log(2 * nu)
    return(bes2)
  } else if (bes == 0) {
    bes3 <-
      tryCatch(
        besselK(z,nu,expon.scaled = T),
        error = function(e)
          NULL
      )
    if (is.null(bes3)) {
      bes4 <- 0.5 * (log(pi) - log(2) - log(nu)) - nu * log(exp(1) * z) + nu *
        log(2 * nu)
      return(bes4)
    } else if (bes3 < 1e-300) {
      bes5 <- 1e-300
      return(bes5)
    } else  {
      return(log(bes3) - z)
    }
  } else {
    return(log(bes))
  }
}
#' @keywords internal
besderiv<-function(nu,z){
  eps=0.001
  bessdev<-(mybessel(abs(nu+eps),z)-mybessel(abs(nu),z))/eps
  return(bessdev)
}

#'Simulated Data
#'
#' This is a simulated dataset with 100 observations from 4 by 3 matrix skew-t distribution.
#'
#' @docType data
#' @usage data(SimX)
#' @keywords data
"SimX"
