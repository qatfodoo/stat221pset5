source("iid_EM.R")
library(mvtnorm)

## Refined EM algorithm

EMref <- function(y, V, I.par=16, A=RoutingMatrix(I.par=16), m.steps=1000) {
  
  ## Compute local datasets
  h <- 5
  w <- 2 * h + 1
  T.y <- dim(y)[2]
  N.t <- T.y - 2 * h
  
  win <- replicate(N.t, matrix(data=0, nrow=7, ncol=w), simplify=F)
  for (t in 1:N.t) {
    win[[t]] <- y[, t:(t + w - 1)]
  }
  
  # Initialize MAP and curvature 
  eta <- log(theta)
  sigma.hat <- diag(1, I.par)
  
  for (t in 1:N.t) {
    
    # Initial guess
    theta <- matrix(data=1, nrow=(I.par + 1), ncol=1)
    q.list <-c()
    
    for (k in 1:m.step) {
      print(k)
      lambda <- theta[1:I.par, ]
      phi <- theta[I.par + 1, ]
      sigma <- phi * diag(lambda^2)
      
      m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
        (y - replicate(T.y, (A %*% lambda), simplify=T))
      R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
      
      # Current expectation and gradient
      min_Q.k <- function(theta) {
        return( - Q.ref(theta, m, R, eta, sigma.hat, V))
      }
      opt.k <- optim(theta, min_Q.k, method=c("L-BFGS-B"), 
                     lower=0.1, upper=Inf, control=list(maxit=1000))
      theta <- opt.k$par
      Q.par <- Q.ref(theta, m, R, eta, sigma.hat, V)
      q.list <- c(q.list, Q.par) # Sequence of conditional expectation
      
      
      
    }
    
    plot(1:length(q.list), q.list)
    
    eta <- log(theta)
    sigma.hat <- 
    
    return(list(theta=theta, q.list=q.list))
  }
  
}

Q.ref <- function(theta, m, R, eta, sigma.hat, V) {
  
  Q.iid <- Q(theta, m, R)
  log.prior <- log(dmvnorm(log(theta)), mean=eta, sigma=(sigma.hat + V))
  
  return(Q.iid + log.prior)
}

g <- function(eta, eta.min1, sigma.hat, V, A, y.t) {
  
  log.prior <- log(dmvnorm(log(theta)), mean=eta, sigma=(sigma.hat + V))
  theta <- exp(eta)
  lambda <- theta[1:16, ]
  phi <- theta[17, ]
  sigma <- phi * diag(lambda^2)
  log.lik <- log(dmvnorm(y.t, mean=(A %*% lambda), sigma=(A %*% sigma %*% t(A))))
  
  return (log.prior + log.lik)

}