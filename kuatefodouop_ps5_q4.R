require(plyr)

## Use iid implementation of EM to derive the local iid EM

data <- read.table("./dat/1router_allcount.dat", header=T, sep=",")
loads <- data[which(grepl("src", data$nme) | grepl("dst", data$nme)), ]

# Compute vectors y of observed data

order <- c("src corp", "src local", "src switch", "src fddi",
           "dst corp", "dst local", "dst switch", "dst fddi")

dl <- dlply(loads, .(time), function(x) {
  df <- x[match(order, x$nme), ]
  return(df$value[1:7]) })

y <- matrix(unlist(dl), nrow=7, byrow=F)

save(y, file="./dat/1rout_y.dat")


load("./out/theta_iid.dat")

