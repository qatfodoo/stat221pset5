library(ggplot2)

## Pset 5 quest 1 & 2

data <- read.table("./dat/1router_allcount.dat", header=T, sep=",")

## Q1

# Limit data to link loads
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

load.df <- data.frame(side=matrix(unlist(strsplit(as.character(loads$nme), split=" ")), ncol=2, byrow=T)[, 1], 
                      node=matrix(unlist(strsplit(as.character(loads$nme), split=" ")), ncol=2, byrow=T)[, 2],
                      value=loads$value, hour.num=hour.num)

# Plot link loads in grid
load.df$side <- factor(load.df$side,
                       levels = c("src", "dst"))
load.df$node <- factor(load.df$node,
                       levels = c("corp", "local", "switch", "fddi"))

p <- ggplot(load.df, aes(hour.num, value, colour=side)) + 
  geom_line(aes(group=side)) + scale_colour_manual(values=c("magenta", "cyan")) +
  labs(title="Link loads for router 1") + xlab("hour of day") + ylab("bytes/sec")
p + facet_grid(node ~ .)

## Q2

load.df2 <- data.frame(side=matrix(unlist(strsplit(as.character(loads$nme), split=" ")), ncol=2, byrow=T)[, 1], 
                      node=matrix(unlist(strsplit(as.character(loads$nme), split=" ")), ncol=2, byrow=T)[, 2],
                      value=loads$value, hour=hour)

## Get restricted data frames
t1.low <- strptime("11:02", "%H:%M")
t1.high <- strptime("11:52", "%H:%M")
t2.low <- strptime("15:05", "%H:%M")
t2.high <- strptime("15:55", "%H:%M")

load.df_t1 = load.df2[which((strptime(load.df2$hour, format="%H:%M") >= t1.low)
                            & (strptime(load.df2$hour, format="%H:%M") <= t1.high)), ]


load.df_t2 = load.df2[which((strptime(load.df2$hour, format="%H:%M") >= t2.low)
                            & (strptime(load.df2$hour, format="%H:%M") <= t2.high)), ]

## Get means and sd

agg.mean_t1 <-aggregate(load.df_t1$value, by=list(load.df_t1$side, load.df_t1$node), 
                    FUN=mean, na.rm=TRUE)
agg.mean_t2 <-aggregate(load.df_t2$value, by=list(load.df_t2$side, load.df_t2$node), 
                        FUN=mean, na.rm=TRUE)

agg.var_t1 <- aggregate(load.df_t1$value, by=list(load.df_t1$side, load.df_t1$node), 
                     FUN=var, na.rm=TRUE)
agg.var_t2 <- aggregate(load.df_t2$value, by=list(load.df_t2$side, load.df_t2$node), 
                        FUN=var, na.rm=TRUE)

## Join and put on log scale
log.stats_t1 <- data.frame(time="time 11:30", side=agg.mean_t1[, 1], node=agg.mean_t1[, 2],
                           log.mean=log(agg.mean_t1[, 3], 10),
                           log.var= log(agg.var_t1[, 3], 10))

log.stats_t2 <- data.frame(time="time 15:30", side=agg.mean_t2[, 1], node=agg.mean_t2[, 2],
                           log.mean=log(agg.mean_t2[, 3], 10),
                           log.var= log(agg.var_t2[, 3], 10))

log.stats <- rbind(log.stats_t1, log.stats_t2)

# Get linear models
times <- c("time 11:30", "time 15:30")
lm_t1 <- lm(log.stats_t1$log.var ~ log.stats_t1$log.mean)
lm_t2 <- lm(log.stats_t2$log.var ~ log.stats_t2$log.mean)
coeff_lm <- data.frame(time=times, int=c(lm_t1$coefficients[1], lm_t2$coefficients[1]),
                       slope=c(lm_t1$coefficients[2], lm_t2$coefficients[2]))

# Slope 1
times <- c("time 11:30", "time 15:30")
int.1_t1 <- mean(log.stats_t1$log.var - log.stats_t1$log.mean)
int.1_t2 <- mean(log.stats_t2$log.var - log.stats_t2$log.mean)
coeff_lm.1 <- data.frame(time=times, int=c(int.1_t1, int.1_t2),
                       slope=c(1, 1))

# Slope 2
times <- c("time 11:30", "time 15:30")
int.2_t1 <- mean(log.stats_t1$log.var - 2 * log.stats_t1$log.mean)
int.2_t2 <- mean(log.stats_t2$log.var - 2 * log.stats_t2$log.mean)
coeff_lm.2 <- data.frame(time=times, int=c(int.2_t1, int.2_t2),
                         slope=c(2, 2))


p <- ggplot(log.stats, aes(log.mean, log.var, shape=interaction(side, node),
                           colour=interaction(side, node))) + 
  geom_point(aes(group=interaction(side, node))) +
  geom_abline(data=coeff_lm, aes(intercept=coeff_lm$int, slope=coeff_lm$slope),
              colour="red", size=1, linetype=1) +
  geom_abline(data=coeff_lm.1, aes(intercept=coeff_lm.1$int, slope=coeff_lm.1$slope),
              colour="blue", size=1, linetype=3) +
  geom_abline(data=coeff_lm.2, aes(intercept=coeff_lm.2$int, slope=coeff_lm.2$slope),
              colour="green", size=1, linetype=2) +
  scale_shape_manual(values=1:8) +
  labs(title="Local variances vs. local means for two time windows") + xlab("log10(mean)") + ylab("log10(var)")
p + facet_grid(. ~ time)

## Same for second data set

data2 <- read.table("./dat/2router_linkcount.dat", header=T, sep=",")

time2 <- gsub("\\)", "", gsub("\\(", "", data2$time))
time2 <- strptime(time2, format="%m/%d/%Y %H:%M:%S")
hour2 <- format(time2, "%H:%M")

data2.df <- data.frame(side=matrix(unlist(strsplit(as.character(data2$nme), split=" ")), ncol=2, byrow=T)[, 1], 
                       node=matrix(unlist(strsplit(as.character(data2$nme), split=" ")), ncol=2, byrow=T)[, 2],
                       value=data2$value, hour=hour2)

## Get restricted data frames
t1.low2 <- strptime("11:02", "%H:%M")
t1.high2 <- strptime("11:52", "%H:%M")
t2.low2 <- strptime("15:05", "%H:%M")
t2.high2 <- strptime("15:55", "%H:%M")

data2.df_t1 = data2.df[which((strptime(data2.df$hour, format="%H:%M") >= t1.low2)
                            & (strptime(data2.df$hour, format="%H:%M") <= t1.high2)), ]


data2.df_t2 = data2.df[which((strptime(data2.df$hour, format="%H:%M") >= t2.low2)
                            & (strptime(data2.df$hour, format="%H:%M") <= t2.high2)), ]

## Get means and sd

agg.mean2_t1 <-aggregate(data2.df_t1$value, by=list(data2.df_t1$side, data2.df_t1$node), 
                        FUN=mean, na.rm=TRUE)
agg.mean2_t2 <-aggregate(data2.df_t2$value, by=list(data2.df_t2$side, data2.df_t2$node), 
                        FUN=mean, na.rm=TRUE)

agg.var2_t1 <- aggregate(data2.df_t1$value, by=list(data2.df_t1$side, data2.df_t1$node), 
                        FUN=var, na.rm=TRUE)
agg.var2_t2 <- aggregate(data2.df_t2$value, by=list(data2.df_t2$side, data2.df_t2$node), 
                        FUN=var, na.rm=TRUE)

## Join and put on log scale
log.stats2_t1 <- data.frame(time="time 11:30", side=agg.mean2_t1[, 1], node=agg.mean2_t1[, 2],
                           log.mean=log(agg.mean2_t1[, 3], 10),
                           log.var= log(agg.var2_t1[, 3], 10))

log.stats2_t2 <- data.frame(time="time 15:30", side=agg.mean2_t2[, 1], node=agg.mean2_t2[, 2],
                           log.mean=log(agg.mean2_t2[, 3], 10),
                           log.var= log(agg.var2_t2[, 3], 10))

log.stats2 <- rbind(log.stats2_t1, log.stats2_t2)

# Get linear models
times <- c("time 11:30", "time 15:30")
lm2_t1 <- lm(log.stats2_t1$log.var ~ log.stats2_t1$log.mean)
lm2_t2 <- lm(log.stats2_t2$log.var ~ log.stats2_t2$log.mean)
coeff2_lm <- data.frame(time=times, int=c(lm2_t1$coefficients[1], lm2_t2$coefficients[1]),
                       slope=c(lm2_t1$coefficients[2], lm2_t2$coefficients[2]))

# Slope 1
times <- c("time 11:30", "time 15:30")
int2.1_t1 <- mean(log.stats2_t1$log.var - log.stats2_t1$log.mean)
int2.1_t2 <- mean(log.stats2_t2$log.var - log.stats2_t2$log.mean)
coeff2_lm.1 <- data.frame(time=times, int=c(int2.1_t1, int2.1_t2),
                         slope=c(1, 1))

# Slope 2
times <- c("time 11:30", "time 15:30")
int2.2_t1 <- mean(log.stats2_t1$log.var - 2 * log.stats2_t1$log.mean)
int2.2_t2 <- mean(log.stats2_t2$log.var - 2 * log.stats2_t2$log.mean)
coeff2_lm.2 <- data.frame(time=times, int=c(int2.2_t1, int2.2_t2),
                         slope=c(2, 2))


p <- ggplot(log.stats2, aes(log.mean, log.var, shape=interaction(side, node),
                           colour=interaction(side, node))) + 
  geom_point(aes(group=interaction(side, node))) +
  geom_abline(data=coeff2_lm, aes(intercept=coeff2_lm$int, slope=coeff2_lm$slope),
              colour="red", size=1, linetype=1) +
  geom_abline(data=coeff2_lm.1, aes(intercept=coeff2_lm.1$int, slope=coeff2_lm.1$slope),
              colour="blue", size=1, linetype=3) +
  geom_abline(data=coeff2_lm.2, aes(intercept=coeff2_lm.2$int, slope=coeff2_lm.2$slope),
              colour="green", size=1, linetype=2) +
  scale_shape_manual(values=1:16) +
  labs(title="Local variances vs. local means for two time windows") + xlab("log10(mean)") + ylab("log10(var)")
p + facet_grid(. ~ time)