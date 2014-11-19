
## Implement the EM algorithm for the iid model and c=2

# Try to implement the EM, using fractional Newton-Raphson
EMiid <- function(y, I.par=16, A=RoutingMatrix(16), m.step=1000) {
  
  T.y <- dim(y)[2]
  # Initial guess
  theta <- matrix(data=c(rep(10000, I.par), 1), nrow=(I.par + 1), ncol=1)
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
EMiid_NR <- function(y, I.par=16, A=RoutingMatrix(16), m.step=1000) {
  
  T.y <- dim(y)[2]
  # Initial guess
  theta <- matrix(data=7, nrow=(I.par + 1), ncol=1)
  q.list <-c()
  
  for (k in 1:m.step) {
    lambda <- theta[1:I.par, ]
    phi <- theta[I.par + 1, ]
    sigma <- phi * diag(lambda^2)
    
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
  
  orig.matlist <- matrix(data=0, nrow=4, ncol=I.par)
  for (i in 1:4) {
    for (k in (1 + (i-1) * 4):(i*4)) {
      orig.matlist[i, k] <- 1
    }
  }
  dest.matlist <- matrix(data=0, nrow=3, ncol=I.par)
  for (i in 1:3) {
    for (k in 1:I.par) {
      if ((k - i) %% 4 == 0) {
        dest.matlist[i, k] <- 1
      }
    }
  }
  A <- matrix(data=0, nrow=7, ncol=I.par)
  A[1:4, ] <- orig.matlist
  A[5:7, ] <- dest.matlist
  
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
  lambda <- theta[1:I.par]
  phi <- theta[I.par + 1]
  sigma <- phi * diag(lambda^2)
  #print(phi)
  
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
