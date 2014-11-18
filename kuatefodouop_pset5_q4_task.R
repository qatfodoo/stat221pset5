source("iid_EM.R")

## question 1.4 on Odyssey

# 10 iterations of each chain, save plot and store chains
if (Sys.getenv("SLURM_JOB_ID") != "") { # Divide computation per tasks
  
  job.id <- as.numeric(Sys.getenv("SLURM_JOB_ID"))
  print(paste("Job id", job.id, sep=": "))
  task.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  print(paste("Task id", task.id, sep=": "))
  #data.id <- task.id %% 2 + 1 # Set of data considered
  #data.name <- data.names[data.id]
  #chain.id <- ceiling(task.id / 2) # Chain id for data set
  
  gc() # garbage collection
  load("./dat/1rout_y.dat")
  
  t1.sim <- as.numeric(Sys.time())
  
  theta <- EMiid(y, m.step=5)
  
  save(theta, "./out/theta_iid.dat")
  
  t2.sim <- as.numeric(Sys.time())
  dt.sim <- (t2.sim - t1.sim) / 60 # dt in min
  print(paste(paste("EM fitting elapsed time (min), task", task.id, sep=" "), dt.sim, sep=": "))

}
