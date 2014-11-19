library(ggplot2)
source("iid_EM.R")
library(reshape2)

## Plotting of the OD estimates

load("./dat/1rout_y.dat")

I.par <- 16
## Compute local datasets
h <- 5
w <- 2 * h + 1
T.y <- dim(y)[2]
N.t <- T.y - 2 * h

theta <- matrix(data=0, nrow=(I.par + 1), ncol=N.t)

for (task.id in 1:30) {
  
  win.id <- seq(from=task.id, to=N.t, by=30)
  load(paste("./out/theta_part_", task.id, ".dat", sep=""))
  for (theta.t in theta.part) {
    theta[, theta.t$t] <- theta.t$theta
  }
}

lambda <- theta[1:I.par, ]

# Retrieve hour data
data <- read.table("./dat/1router_allcount.dat", header=T, sep=",")
loads <- data[which(grepl("src", data$nme) | grepl("dst", data$nme)), ]

time <- gsub("\\)", "", gsub("\\(", "", loads$time))
time <- strptime(time, format="%m/%d/%Y %H:%M:%S")
hour <- format(time, "%H:%M")
hour.num <- sapply(strsplit(hour,":"),
                   function(x) {
                     x <- as.numeric(x)
                     x[1]+x[2]/60
                   }
)
window.t <- unique(hour.num)[(h + 1):(T.y - h)]
c("src corp", "src local", "src switch", "src fddi",
  "dst corp", "dst local", "dst switch", "dst fddi")
# Origin dest 
orig.order <- c(rep("corp", 4), rep("local", 4), rep("switch", 4), rep("fddi", 4))
dest.order <- rep(c("corp", "local", "switch", "fddi"), 4)

lambda.vec <- c(lambda)
windowt.vec <- rep(window.t, each=I.par)
src.vec <- rep(orig.order, times=N.t)
dest.vec <- rep(dest.order, times=N.t)

lambda.est <- data.frame(lambda.t=lambda.vec, hour.num=windowt.vec, src=src.vec,
                         dest=dest.vec)

# Reorder for right plotting order
lambda.est$src <- factor(lambda.est$src,
                       levels = c("corp", "local", "switch", "fddi"))
lambda.est$dest <- factor(lambda.est$dest,
                       levels = c("fddi", "switch", "local", "corp"))

p <- ggplot(lambda.est, aes(hour.num, lambda.t, colour=src)) + 
  geom_line() +
  labs(title="Link loads for router 1") + xlab("hour of day") + ylab("bytes/sec")
p + facet_grid(src ~ dest)


## Plot linked loads predicted versus observed
A <- RoutingMatrix()

# Predicted load links
yhat.t <- matrix(data=0, nrow=8, ncol=N.t)
for (t in 1:N.t) {
  yhat.t[1:7, t] <- A %*% lambda[, t]
  yhat.t[8, t] <- sum(yhat.t[1:4, t]) - sum(yhat.t[5:7, t])
}

# y moving averages
y.movav <- matrix(data=0, nrow=8, ncol=N.t)
for (t in 1:N.t) {
  y.movav[1:7, t] <- rowMeans(y[, t:(t + w - 1)])
  y.movav[8, t] <- sum(y.movav[1:4, t]) - sum(y.movav[5:7, t])
}

side.order <- c("src", "src", "src", "src",
           "dst", "dst", "dst", "dst")
node.order <- c("corp", "local", "switch", "fddi",
                "corp", "local", "switch", "fddi")

yhat.vec <- c(yhat.t)
ymovav.vec <- c(y.movav)
side.vec <- rep(side.order, times=N.t)
node.vec <- rep(node.order, times=N.t)

mean.link <- data.frame(yhat=yhat.vec, ymovav=ymovav.vec, hour.num=windowt.vec,
                        side=side.vec, node=node.vec)

# Reorder for right plotting order
mean.link$side <- factor(mean.link$side,
                         levels = c("dst", "src"))
mean.link$node <- factor(mean.link$node,
                          levels = c("fddi", "switch", "local", "corp"))

ml <- melt(mean.link, id.vars=c("hour.num", "side", "node"))

p <- ggplot(ml, aes(hour.num, value, colour=variable)) + 
  geom_line() +
  labs(title="Link loads for router 1") + xlab("hour of day") + ylab("bytes/sec")
p + facet_grid(side ~ node)
