source("iid_EM.R")

## question 1.4 on Odyssey

# 10 iterations of each chain, save plot and store chains
if (Sys.getenv("SLURM_JOB_ID") != "") { # Divide computation per tasks
  
  job.id <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
  print(paste("Job id", job.id, sep=": "))
  task.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  print(paste("Task id", task.id, sep=": "))
  
  gc() # garbage collection
  load("./dat/1rout_y.dat")
  
  I.par <- 16
  ## Compute local datasets
  h <- 5
  w <- 2 * h + 1
  T.y <- dim(y)[2]
  N.t <- T.y - 2 * h
  
  win <- replicate(N.t, matrix(data=0, nrow=7, ncol=w), simplify=F)
  for (t in 1:N.t) {
    win[[t]] <- y[, t:(t + w - 1)]
  }
  
  t1.sim <- as.numeric(Sys.time())
  
  m.step < 2000 # num of steps for EM
  
  win.id <- seq(from=task.id, to=N.t, by=30)
  theta.part <- replicate(length(win.id), list(t=0,
                            theta=matrix(data=0, nrow=I.par + 1, ncol=1),
                            q.list=rep(0, m.step)), simplify=F)
  
  for (i in 1:length(win.id)) {
    t <- win.id[i]
    y.t <- win[[t]]
    theta.part[[i]]$t <- t
    EM <- EMiid(y.t, m.step=m.step)
    theta.part[[i]]$theta <- EM$theta
    theta.part[[i]]$q.list <- EM$q.list
  }
  
  save(theta.part, file=paste("./out/theta_part_", task.id, ".dat", sep=""))
  
  t2.sim <- as.numeric(Sys.time())
  dt.sim <- (t2.sim - t1.sim) / 60 # dt in min
  print(paste(paste("EM fitting elapsed time (min), task", task.id, sep=" "), dt.sim, sep=": "))

}
