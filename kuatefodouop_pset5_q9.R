require(plyr)

## Format data for router 2

loads <- read.table("./dat/2router_linkcount.dat", header=T, sep=",")
loads <- data[which(grepl("src", data$nme) | grepl("dst", data$nme)), ]

# Compute vectors y of observed data

order <- c("ori router5", "ori r4-local", "ori switch", "ori r4-others",
           "ori gw1", "ori gw2", "ori gw3", "ori gw-others",
           "dst router5", "dst r4-local", "dst switch", "dst r4-others",
           "dst gw1", "dst gw2", "dst gw3", "dst gw-others")

dl <- dlply(loads, .(time), function(x) {
  df <- x[match(order, x$nme), ]
  return(df$value[1:15]) })

y <- matrix(unlist(dl), nrow=15, byrow=F)

save(y, file="./dat/2rout_y.dat")
