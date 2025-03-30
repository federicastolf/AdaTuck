rm(list = ls())

library(tensorBF)
library(parallel)
source("Rfuncts/funcAdaTuck.R")
source("Rfuncts/AdaTuckgibbs.R")

#------------------------------# scenario 1 #--------------------------------#
Y_dim = c(30, 30, 10) #n1, n2, n3
core_dim = c(5, 5 ,5) # R1, R2, R3
N_sampl = 12000

param = list(a_eta = 1, b_eta = 0.3, a_theta = 2, b_theta=2,
             theta_inf = rep(0.05, length(Y_dim)), alpha = rep(3, length(Y_dim)),
             start_adapt = 200, c0_adapt = -1, c1_adapt = -5*10^(-4),
             a_tau = 2, b_tau = 2, a_rho = 10, b_rho = 10)
Nsim = 20

set.seed(4373)
seeds_g = sample.int(9000, Nsim)
all_data50 = all_data30 = vector("list", Nsim)
for(i in 1:Nsim){
  all_data50[[i]] = simSynthData(Y_dim = Y_dim, core_dim = core_dim,
                                 my_seed = seeds_g[i], missY = 0.5, sparseG = 0.4)
  all_data30[[i]] = simSynthData(Y_dim = Y_dim, core_dim = core_dim,
                                 my_seed = seeds_g[i], missY = 0.3, sparseG = 0.4)
}

## compute AdaTuck 
gibbs_out30 = parallel::mclapply(1:Nsim, gibbs_paralleltc, mc.cores = 20,
                                all_data = all_data30, N_sampl = N_sampl, param = param,
                                seeds_g = seeds_g)

gibbs_out50 = parallel::mclapply(1:Nsim, gibbs_paralleltc, mc.cores = 20,
                                 all_data = all_data50, N_sampl = N_sampl, param = param,
                                 seeds_g = seeds_g)

## Monte Carlo median of multi-rank 
burn_in = seq(8000+1, N_sampl)
rankstar5mmis3 = rankstar5mmis5 = matrix(NA, Nsim, 3)
for(s in 1:Nsim){
  rankstar5mmis3[s,] = apply(gibbs_out30[[s]]$b$Rstar[,burn_in], 1, median)
  rankstar5mmis5[s,] = apply(gibbs_out50[[s]]$b$Rstar[,burn_in], 1, median)
}

# compute MSE
MSE5miss30 = computeMSEadaTuck(all_data30, gibbs_out30, burn_in)
MSE5miss50 = computeMSEadaTuck(all_data50, gibbs_out50, burn_in)

#---# run TensorBf #---#
opts <- getDefaultOpts();
opts$iter.burnin = 8000
opts$iter.sampling = 4000
opts$verbose = 0

MSEBF5miss3 = parallel::mclapply(1:Nsim, TensorBF_paralleltc, mc.cores = 20,
                                 all_data = all_data30, opts = opts)

MSEBF5miss5 = parallel::mclapply(1:Nsim, TensorBF_paralleltc, mc.cores = 20,
                                 all_data = all_data50, opts = opts)

#------------------------------# scenario 2 #--------------------------------#
Y_dim = c(50, 40, 6) #n1, n2, n3
core_dim = c(10, 7 ,3) # R1, R2, R3

set.seed(8273)
seeds_g = sample.int(9000, Nsim)
all_data50 = all_data30 = vector("list", Nsim)
for(i in 1:Nsim){
  all_data50[[i]] = simSynthData(Y_dim = Y_dim, core_dim = core_dim,
                                 my_seed = seeds_g[i], missY = 0.5, sparseG = 0.4)
  all_data30[[i]] = simSynthData(Y_dim = Y_dim, core_dim = core_dim,
                                 my_seed = seeds_g[i], missY = 0.3, sparseG = 0.4)
}

## compute AdaTuck 
gibbs_out30 = parallel::mclapply(1:Nsim, gibbs_paralleltc, mc.cores = 20,
                                 all_data = all_data30, N_sampl = N_sampl, param = param,
                                 seeds_g = seeds_g)
gibbs_out50 = parallel::mclapply(1:Nsim, gibbs_paralleltc, mc.cores = 20,
                                 all_data = all_data50, N_sampl = N_sampl, param = param,
                                 seeds_g = seeds_g)

## Monte Carlo median of multi-rank 
rankstar10mmis3 = rankstar10mmis5 = matrix(NA, Nsim, 3)
for(s in 1:Nsim){
  rankstar10mmis3[s,] = apply(gibbs_out30[[s]]$b$Rstar[,burn_in], 1, median)
  rankstar10mmis5[s,] = apply(gibbs_out50[[s]]$b$Rstar[,burn_in], 1, median)
}

# compute MSE
MSE10miss30 = computeMSEadaTuck(all_data30, gibbs_out30, burn_in)
MSE10miss50 = computeMSEadaTuck(all_data50, gibbs_out50, burn_in)

#---# run TensorBf #---#
MSEBF10miss3 = parallel::mclapply(1:Nsim, TensorBF_paralleltc, mc.cores = 20,
                                  all_data = all_data30, opts = opts)

MSEBF10miss5 = parallel::mclapply(1:Nsim, TensorBF_paralleltc, mc.cores = 20,
                                  all_data = all_data50,opts = opts)
