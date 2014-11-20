
## EM functions

# iid and local EM


## Implement the EM algorithm for the iid model and c=2

# EM using optim
iid_EM <- function(y, c.par=2, I.par=16, A=RoutingMatrix(16), m.step=1000) {
  
  T.y <- dim(y)[2]
  # Initial guess
  theta <- matrix(data=c(rep(1000, I.par), 1), nrow=(I.par + 1), ncol=1)
  q.list <-c()
  
  for (k in 1:m.step) {
    print(k)
    lambda <- theta[1:I.par, ]
    phi <- theta[I.par + 1, ]
    sigma <- phi * diag(lambda^c.par)
    
    m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
      (y - replicate(T.y, (A %*% lambda), simplify=T))
    R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
    
    # Current expectation and gradient
    min_Q.k <- function(theta) {
      return( - Q(theta, m, R))
    }
    
    opt.k <- optim(theta, min_Q.k, method=c("L-BFGS-B"), 
                   lower=0.1, upper=Inf, control=list(maxit=10000))
    theta <- opt.k$par
    Q.par <- Q(theta, m, R)
    q.list <- c(q.list, Q.par) # Sequence of conditional expectation
    
  }
  
  plot(1:length(q.list), q.list)
  
  return(list(theta=theta, q.list=q.list))
  
}


# Try to implement the EM, using fractional Newton-Raphson
iid_EM_NR <- function(y, c.par=2, I.par=16, A=RoutingMatrix(16), m.step=1000) {
  
  T.y <- dim(y)[2]
  # Initial guess
  theta <- matrix(data=7, nrow=(I.par + 1), ncol=1)
  q.list <-c()
  
  for (k in 1:m.step) {
    lambda <- theta[1:I.par, ]
    phi <- theta[I.par + 1, ]
    sigma <- phi * diag(lambda^c.par)
    
    m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
      (y - replicate(T.y, (A %*% lambda), simplify=T))
    R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
    
    a <- diag(R) + (1 / T.y) * diag(m %*% t(m))
    b <- rowMeans(m)
    
    theta <- FracNR(a, b, phi)
    
    Q.par <- Q(theta, m, R)
    q.list <- c(q.list, Q.par) # Sequence of conditional expectation
    
  }
  
  plot(1:length(q.list), q.list)
  
  return(theta)
  
}

# Routing matrix of router1 A
RoutingMatrix <- function(I.par=16) {
  
  n <- log(I.par, 2)
  orig.matlist <- matrix(data=0, nrow=n, ncol=I.par)
  for (i in 1:n) {
    for (k in (1 + (i-1) * n):(i*n)) {
      orig.matlist[i, k] <- 1
    }
  }
  dest.matlist <- matrix(data=0, nrow=(n - 1), ncol=I.par)
  for (i in 1:(n - 1)) {
    for (k in 1:I.par) {
      if ((k - i) %% n == 0) {
        dest.matlist[i, k] <- 1
      }
    }
  }
  A <- matrix(data=0, nrow=(2 * n - 1), ncol=I.par)
  A[1:n, ] <- orig.matlist
  A[(n + 1):(2 * n - 1), ] <- dest.matlist
  
  return (A)
}

# Fractional N-R on lambda, phi
FracNR <- function(a, b, phi0, I.par=16, n.step=1, tol=1e-2) {
  
  phi <- phi0
  lambda <- matrix(data=0, nrow=1, ncol=I.par)
  err <- 1
  n.iter <- 0
  
  for (k in 1:n.step) {
    #while (abs(err) > tol) {
    n.iter <- n.iter + 1
    
    # Optimal lambda for phi
    lambda <- (sqrt(b^2 + 4 * phi * a) - b) / (2 * phi)
    
    # f component
    f.k <- f(a, b, lambda, phi)
    f_Iplus1 <- f.k[I.par + 1]
    err <- f_Iplus1
    
    # Jacobian
    F.diag <- c(4 * phi / lambda + 2 * b, 0)
    F <- diag(F.diag)
    F[I.par + 1, 1:I.par] <- b / lambda
    F[1:I.par, I.par + 1] <- 2 * lambda^2
    # Compute F inverse
    F.inv <- solve(F)
    
    # NR step
    step <- - F.inv[I.par + 1, I.par + 1] * f_Iplus1
    t <- 1
    phi.next <- -1
    while (phi.next < 0) {
      phi.next <- phi + (1 / t) * step
      t <- t + 1 # Decrease the step while negative
    }
    phi <- phi.next
    
  }
  #lambda <- (sqrt(b^2 + 4 * phi * a) - b) / (2 * phi) # Updated lambdas
  
  #print("number of iters: ")
  #print(n.iter)
  return(matrix(data=c(lambda, phi), nrow=(I.par + 1), ncol=1))
  
}

# Conditional expectation function Q
Q <- function(theta, m, R) {
  
  T.y <- dim(m)[2]
  lambda <- theta[1:(length(theta) - 1)]
  phi <- theta[length(theta)]
  sigma <- phi * diag(lambda^2)
  
  res <- - (T.y / 2) * log(det(sigma))
  sig.inv <- solve(sigma)
  P <- sig.inv %*% R
  for (k in 1:dim(P)[2]) {
    res <- res - (T.y / 2) * P[k, k]
  }
  for (t in 1:T.y) {
    res <- res - (1 / 2) * t(m[, t] - lambda) %*% sig.inv %*% (m[, t] - lambda)
  }
  
  return(as.numeric(res))
  
}

# f function on the parameter
f <- function(a, b, lambda, phi) {
  
  f.head <- phi * lambda^2 + b * lambda - a
  f.last <- sum((lambda - b) / lambda)
  
  return(c(f.head, f.last))
}


## Implement refined model
library(mvtnorm)
library(numDeriv)

smoothed_EM <- function(y, c.par=2, V=diag(1, 17), I.par=16, A=RoutingMatrix(I.par=16), m.step=1000) {
  
  ## Compute local datasets
  h <- 5
  w <- 2 * h + 1
  n <- dim(y) + 1
  T.y <- dim(y)[2]
  N.t <- T.y - 2 * h
  
  win <- replicate(N.t, matrix(data=0, nrow=(n - 1), ncol=w), simplify=F)
  
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
    theta.t <- exp(eta)
    q.list <-c()
    
    for (k in 1:m.step) {
      print(k)
      lambda <- theta.t[1:I.par, ]
      phi <- theta.t[I.par + 1, ]
      sigma <- phi * diag(lambda^c.par)
      
      m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
        (y.t - replicate(w, (A %*% lambda), simplify=T))
      R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
      
      # Current expectation and gradient
      min_Q.k <- function(theta.t) {
        return( min(- Q.ref(theta.t, m, R, eta, sigma.hat, V),
                    .Machine$double.xmax))
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
    sigma.hat <- - solve(hessian(g.t, eta))
    
    theta.est[[t]] <- list(theta.t=theta.t, q.list=q.list)
  }
  
  return(theta.est)
  
}

Q.ref <- function(theta, m, R, eta, sigma.hat, V) {
  
  Q.iid <- Q(theta, m, R)
  log.prior <- dmvnorm(c(log(theta)), mean=c(eta), sigma=(sigma.hat + V), log=T)
  
  return(Q.iid + log.prior)
}

g <- function(eta, eta.min1, sigmahat.min1, V, A, y.t) {
  
  log.prior <- dmvnorm(c(eta), mean=c(eta.min1), sigma=(sigmahat.min1 + V), log=T)
  theta <- exp(eta)
  lambda <- theta[1:(dim(theta)[1] - 1), ]
  phi <- theta[dim(theta)[1], ]
  sigma <- phi * diag(lambda^2)
  log.lik <- sum(dmvnorm(t(y.t), mean=c(A %*% lambda), sigma=(A %*% sigma %*% t(A)), log=T))
  
  return (log.prior + log.lik)
  
}