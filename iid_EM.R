
## Implement the EM algorithm for the iid model and c=2



EMiid <- function(y, I.par=16, A=RoutingMatrix(16), m.step=100) {
  
  T.y <- dim(y)[2]
  # Initial guess
  theta <- matrix(data=1, nrow=(I.par + 1), ncol=1)
  
  for (k in 1:m.step) {
    lambda <- theta[1:I.par, ]
    phi <- theta[I.par + 1, ]
    sigma <- phi * diag(lambda^2)
    
    m <- lambda + sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*%
      (y - matrix(data=(A %*% lambda), nrow=7, ncol=T.y, byrow=F))
    R <- sigma - sigma %*% t(A) %*% solve(A %*% sigma %*% t(A)) %*% A %*% sigma
    
    a <- diag(R) + (1 / T.y) * diag(m %*% t(m))
    b <- rowMeans(m)
    
    theta <- FracNR(a, b, phi)
    
  }
  
  return(theta)
  
}

# Routing matrix of router1 A
RoutingMatrix <- function(I=16) {
  
  orig.matlist <- replicate(4, rep(0, I), simplify=F)
  for (i in 1:4) {
    for (k in (1 + (i-1) * 4):(i*4)) {
      orig.matlist[[i]][k] <- 1
    }
  }
  dest.matlist <- replicate(3, rep(0, I), simplify=F)
  for (i in 1:3) {
    for (k in 1:I) {
      if ((k - i) %% 4 == 0) {
        dest.matlist[[i]][k] <- 1
      }
    }
  }
  A <- matrix(data=0, nrow=7, ncol=I)
  A[1:4, ] <- matrix(unlist(orig.matlist), ncol=16, byrow=T)
  A[5:7, ] <- matrix(unlist(dest.matlist), ncol=16, byrow=T)
  
  return (A)
}

# Fractional N-R on lambda, phi
FracNR <- function(a, b, phi0, I.par=16, tol=1e-2) {
  
  phi <- phi0
  lambda <- matrix(data=0, nrow=1, ncol=I)
  err <- 1
  n.iter <- 0
  
  while (abs(err) > tol) {
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
  lambda <- (sqrt(b^2 + 4 * phi * a) - b) / (2 * phi) # Updated lambdas
  
  print("number of iters: ")
  print(n.iter)
  return(matrix(data=c(lambda, phi), nrow=(I.par + 1), ncol=1))
  
}

# f function on the parameter
f <- function(a, b, lambda, phi) {
  
  f.head <- phi * lambda^2 + b * lambda - a
  f.last <- sum((lambda - b) / lambda)
  
  return(c(f.head, f.last))
}