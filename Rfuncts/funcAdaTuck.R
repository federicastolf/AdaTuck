library(rTensor)

###############################################################################
#-------------------------# functions for simulations #-----------------------#

simSynthData = function(Y_dim, core_dim,  my_seed, eta2 = 0.1, missY = 0, 
                        sparseG = 0, binary = FALSE){
  "---------------------------------------------------------------------
   It generate synthetic tensor data constructed from a Tucker decomposition
   arguments:
   Y_dim: dimensione of the array
   core_dim: dimensione of the core array
   eta2: variance of noisy tensor
   missY: percentage of missing values in the tensor
   sparseG: percentage of sparsity in the core tensor
   binary: TRUE if binary tensor data
  ----------------------------------------------------------------------"
  n_mode = length(Y_dim)
  U = vector("list", n_mode)
  set.seed(my_seed)
  for(k in 1:n_mode){
    nu = rgamma(core_dim[k],2,2)
    U[[k]] = matrix(rnorm(Y_dim[k]*core_dim[k], mean = 0, sd = 1), Y_dim[k], core_dim[k])
    U[[k]] = t(apply(U[[k]],1, function(x) x/nu))
  }
  G = as.tensor(array(rnorm(prod(core_dim), mean = 0, sd = 1), dim = core_dim))
  
  # core-tensor sparsity
  if(sparseG > 0){
    idx_0 = sample(x = 1:prod(core_dim), size = sparseG*prod(core_dim), replace=F)
    G@data[idx_0] = 0
  }
  
  if(binary==TRUE){
    Z = rTensor::ttl(G, U, ms = c(1,2,3))@data
    Y = array(as.numeric(Z>0), dim = Y_dim)
    # missing values
    yNA = Y
    if(missY>0){
      idx_na = sample(x = 1:prod(Y_dim), size = missY*prod(Y_dim), replace=F)
      yNA[idx_na] = NA
    }
    noisyY = NULL
  }
  else{
    Y = rTensor::ttl(G, U, ms = c(1,2,3))@data
    # add noise
    noisyY = Y + sqrt(eta2)*array(rnorm(prod(Y_dim)), dim = Y_dim)
    # missing values
    yNA = noisyY
    if(missY>0){
      idx_na = sample(x = 1:prod(Y_dim), size = missY*prod(Y_dim), replace=F)
      yNA[idx_na] = NA
    }
  }
  
  return(list(Y = Y, U = U, G = G@data, noisyY = noisyY, yNA = yNA))
}


gibbs_parallel = function(s, all_data, N_sampl, param, seeds_g){
  y1 = all_data[[s]]$yNA
  post_tucker = AdaTuck(y1, seeds_g[s], N_sampl, param)
  return(post_tucker)
}

gibbs_paralleltc = function(s, all_data, N_sampl, param, seeds_g){
  b = tryCatch(gibbs_parallel(s, all_data, N_sampl, param, seeds_g))
  c = "check"
  list(b = b, other = c)
}

computeMSEadaTuck = function(all_data, gibbs_out, burn_in){
  Nsim = length(all_data)
  MSE_adatuck = rep(NA, Nsim)
  for(s in 1:Nsim){
    y = all_data[[s]]$noisyY
    ind.na = which(is.na(all_data[[s]]$yNA)) 
    tucker_post = gibbs_out[[s]]$b
    MSE = array(0, length(burn_in))
    for(i in 1:length(burn_in)){
      err_sample = array(rnorm(prod(dim(y)), 0, 1/sqrt(tucker_post$eta_sq_inv[burn_in[i]])), 
                         dim = dim(y))
      Y_rc = ttl(tucker_post$g[[burn_in[i]]], 
                 list_matrix(tucker_post$u, burn_in[i], reverse=FALSE), ms = c(1,2,3))@data
      Y_rc2 = Y_rc + err_sample
      MSE[i] = mean((Y_rc2[ind.na] - y[ind.na])^2)
    }
    MSE_adatuck[s] = mean(MSE)
  }
  return(MSE_adatuck)
}

###############################################################################
#--# functions for predicting missing values using the competitor TensorBF #--#

TensorBF_parallel = function(s, all_data, opts){
  y1 = all_data[[s]]$yNA
  postBF = tensorBF(y1, opts = opts)
  yALL = all_data[[s]]$noisyY
  ind.na = which(is.na(all_data[[s]]$yNA))
  nsamples = length(postBF$posterior$X)
  MSE = rep(0, nsamples)
  for(i in 1:nsamples){
    YRec = rec_tensorBF(postBF$posterior, i)
    MSE[i] = mean((YRec[ind.na] - yALL[ind.na])^2)
  }
  cat(s)
  return(mean(MSE))
}

TensorBF_paralleltc = function(s, all_data, opts){
  b = tryCatch(TensorBF_parallel(s, all_data, opts))
  c = "check"
  list(b = b, other = c)
}

remake.tensor.mat <- function(model){
  return(make.tensor.mat(model))
}

make.tensor.mat <- function(model){
  K = ncol(model$X)
  kk = 1:K
  Yestim <- list()
  Yestim <- 0
  for(k in kk)
    Yestim <- Yestim + outer(outer(model$X[,k],model$W[,k]),model$U[,k])
  return(Yestim)
}

getPosteriorSample = function(post,s){
  samp <- list()
  samp$X <- post$X[[s]]
  samp$U <- post$U[[s]]
  samp$W <- post$W[[s]]
  samp$tau <- post$tau[[s]]
  samp$alpha <- post$alpha[[s]]
  samp$alphaU <- post$alphaU[[s]]
  samp$Z <- post$Z[[s]]
  return(samp)
}

rec_tensorBF = function(post, s){
  md = getPosteriorSample(post,s)
  YRec = remake.tensor.mat(md)
  return(YRec)
}
