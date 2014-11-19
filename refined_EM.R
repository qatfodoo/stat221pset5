source("iid_EM.R")
library(mvtnorm)
library(numDeriv)

## Refined EM algorithm

EMref <- function(y, V=diag(1, 17), I.par=16, A=RoutingMatrix(I.par=16), m.step=1000) {
  
  ## Compute local datasets
  h <- 5
  w <- 2 * h + 1
  T.y <- dim(y)[2]
  N.t <- T.y - 2 * h
  
  win <- replicate(N.t, matrix(data=0, nrow=7, ncol=w), simplify=F)
  
  for (t in 1:5) {#N.t) {
    win[[t]] <- y[, t:(t + w - 1)]
  }
  
  # theta for each t
  theta.est <- replicate(N.t, list(theta.t=matrix(data=0, nrow=(I.par + 1), ncol=1),
                                   q.list=rep(0, m.step)), simplify=F)
  
  # Initialize MAP and curvature 
  eta <- log(matrix(data=c(rep(1000, I.par), 1), nrow=(I.par + 1), ncol=1))
  sigma.hat <- diag(1, (I.par + 1))
  
  for (t in 1:N.t) {
    
    # windowed obs
    y.t <- win[[t]]
    
    # Initial guess
    theta.t <- matrix(data=c(rep(1000, I.par), 1), nrow=(I.par + 1), ncol=1)
    q.list <-c()
    
    for (k in 1:m.step) {
      print(k)
      lambda <- theta.t[1:I.par, ]
      phi <- theta.t[I.par + 1, ]
      sigma <- phi * diag(lambda^2)
      
      m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
        (y.t - replicate(w, (A %*% lambda), simplify=T))
      R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
      
      # Current expectation and gradient
      min_Q.k <- function(theta.t) {
        return( - Q.ref(theta.t, m, R, eta, sigma.hat, V))
      }
      opt.k <- optim(theta.t, min_Q.k, method=c("L-BFGS-B"), 
                     lower=0.1, upper=Inf, control=list(maxit=1000))
      theta.t <- opt.k$par
      Q.par <- Q.ref(theta.t, m, R, eta, sigma.hat, V)
      q.list <- c(q.list, Q.par) # Sequence of conditional expectation
      
    }
    
    plot(1:length(q.list), q.list)
    
    # Update eta and sigma
    eta.min1 <- eta
    sigmahat.min1 <- sigma.hat
    eta <- log(theta.t)
    g.t <- function(e) {
      return(g(e, eta.min1, sigmahat.min1, V, A, y.t))
    }
    sigma.hat <- solve(hessian(g.t, eta))
    
    theta.est[[t]] <- list(theta.t=theta.t, q.list=q.list)
  }
  
}

Q.ref <- function(theta, m, R, eta, sigma.hat, V) {
  
  Q.iid <- Q(theta, m, R)
  log.prior <- dmvnorm(c(log(theta)), mean=c(eta), sigma=(sigma.hat + V), log=T)
  
  return(Q.iid + log.prior)
}

g <- function(eta, eta.min1, sigmahat.min1, V, A, y.t) {
  
  log.prior <- dmvnorm(c(log(theta)), mean=c(eta), sigma=(sigmahat.min1 + V), log=T)
  theta <- exp(eta)
  lambda <- theta[1:16, ]
  phi <- theta[17, ]
  sigma <- phi * diag(lambda^2)
  log.lik <- sum(dmvnorm(t(y.t), mean=c(A %*% lambda), sigma=(A %*% sigma %*% t(A)), log=T))
  
  return (log.prior + log.lik)

}