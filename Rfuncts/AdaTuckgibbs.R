
library(rTensor)
library(mvnfast)
library(GIGrvg)
library(Hmisc)
library(truncnorm)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("Rfuncts/mat_mult.cpp")
source("Rfuncts/auxiliaryGibbs.R")

AdaTuck = function(y, my_seed, N_sampl, param, binary = FALSE, Rinit = NULL){
  "------------------------- --------------------------------------------------
  Gibbs sampler to obtain N_sampl posterior samples for AdaTuck model.
  ----------------------------------------------------------------------------"
  # checking for NA
  NAflag = sum(is.na(y))
  if(NAflag>0){ # there are NA values
    if(binary==FALSE){
      res = gibbsAdaTuckNA(y, my_seed, N_sampl, param, Rinit)
    }
    else{
      res = gibbsProbitAdaTuckNA(y, my_seed, N_sampl, param, Rinit)
    }
  }
  # no NA values
  else{
    if(binary==FALSE){
      res = gibbsAdaTuck(y, my_seed, N_sampl, param, Rinit)
    }
    else{
      res = gibbsProbitAdaTuck(y, my_seed, N_sampl, param, Rinit)
    }
  }
  return(res)
}


#----------------------------------# AdaTuck #---------------------------------#

gibbsAdaTuck = function(y, my_seed, N_sampl, param, Rinit){
  set.seed(my_seed)
  y = as.tensor(y)
  n_mode = y@num_modes # K
  n = y@modes # (n_1,..,n_K)
  unif_rvs = runif(N_sampl) # needed for adaptation
  
  ############### VARIABLES TO UPDATE ###################
  R = matrix(NA, nrow = n_mode, ncol = N_sampl) # numbers of factors
  Rstar = matrix(NA, nrow = n_mode, ncol = N_sampl) # number of active factors
  g = vector("list", N_sampl) # core tensor
  u = vector("list", n_mode) # factor matrices
  theta_inv = vector("list", n_mode)   # mode-specific factor matrices precisions
  omega = vector("list", n_mode) # stick-breaking weights
  s_augmented = vector("list", n_mode) # augmented data
  v = vector("list", n_mode) # stick-breaking construction beta variables
  for(k in 1:n_mode){
    u[[k]] = theta_inv[[k]] = omega[[k]] = s_augmented[[k]] = vector("list", N_sampl)
  }
  eta_sq_inv = rep(NA, N_sampl) # error variances
  tau = rep(NA, N_sampl) # global-scale variance core tensor
  nu = vector("list", n_mode) # local-scale variance core tensor
  rho = vector("list", n_mode) # core tensor shrinkage parameters
  
  ################# INITIALIZATION   ######################
  for(k in 1:n_mode){
    if(is.null(Rinit[k])){
      R[k,1] = floor((sort(n)[1] + n[k])/3)
    }
    else{
      R[k,1] = Rinit[k]
    }
    Rstar[k,1] = R[k,1] -1
    u[[k]][[1]] = matrix(rnorm(n[k]*R[k,1]), n[k], R[k,1])
    theta_inv[[k]][[1]] = rep(1, R[k,1])
    s_augmented[[k]][[1]] = rep(R[k,1], R[k,1])
    omega[[k]][[1]] = rep(1/R[k,1], R[k,1])
  }
  g[[1]] = as.tensor(array(rnorm(prod(n), sd=0.5), dim = R[,1]))
  eta_sq_inv[1] = 1
  tau[1] = 1
  nu[[1]] = rho[[1]] = array(1, dim = R[,1])
  
  ############### GIBBS SAMPLING ########################
  t0 = proc.time()
  for (l in 2:N_sampl){
    
    for(k in 1:n_mode){
      u[[k]][[l]] = matrix(NA, n[k], R[k,l-1])
    }
    
    ##  sample u
    for(k in 1:n_mode){
      yk = k_unfold(y, m = k)@data
      G_k = k_unfold(g[[l-1]], m=k)@data
      fctu = eigenMapMatMult(G_k, t(kronecker_list(listmat_u(u, k, l))))
      V_u = chol2inv(chol(diag(theta_inv[[k]][[l-1]]) + (fctu %*% t(fctu))*eta_sq_inv[[l-1]]))
      for(i in 1:n[k]){
        mu_u =  (V_u %*% (fctu %*%  yk[i,]))*eta_sq_inv[[l-1]]
        u[[k]][[l]][i,] = mvnfast::rmvn(1, mu_u, V_u)
      }
    }
    
    ## sample G 
    W = kronecker_list(list_matrix(u, l, reverse=TRUE))
    g[[l]] = sampleG(W, tau[l-1], nu[[l-1]], y, R[,l-1], eta_sq_inv[l-1])
    
    # sample tau
    tau[l] = sampleTau(param, R[,l-1], g[[l]], nu[[l-1]])
    
    # sample nu
    nu[[l]] = array(GIGrvg::rgig(prod(R[,l-1]), 1/2, vec(g[[l]])^2/tau[l],
                                 c(rho[[l-1]])^2), dim = R[,l-1])
    # sample rho
    rho[[l]] = array(rgamma(prod(R[,l-1]), param$a_rho + 1, param$b_rho + 
                              abs(vec(g[[l]]))/sqrt(tau[l])), dim = R[,l-1])
    
    ## sample eta
    Y_rc = ttl(g[[l]], list_matrix(u, l, reverse=FALSE), ms = c(1,2,3))@data
    eta_sq_inv[l] = rgamma(1, shape = param$a_eta + 0.5*prod(n), rate = param$b_eta + 
                             0.5*sum((y@data - Y_rc)^2))
    
    # sample s_augmented and active variables
    active = vector("list", n_mode)
    for(k in 1:n_mode){
      s_augmented[[k]][[l]] = sampleSag(u[[k]][[l]], R[k,l-1], param, k, n, 
                                        omega[[k]][[l-1]])
      
      # update latent multi-rank
      active[[k]] = which(s_augmented[[k]][[l]]>c(1:R[k,l-1]))
      if(length(active[[k]])==0){
        active[[k]] = c(1,2)
      }
      Rstar[k,l] = length(active[[k]])
      
      # stick breaking elements
      sampleOm = sampleOmega(s_augmented[[k]][[l]], R[k,l-1], param, k)
      omega[[k]][[l]] = sampleOm$omega
      v[[k]] = sampleOm$v
      
      ### loading precisions
      theta_inv[[k]][[l]] = rep(param$theta_inf[k]^(-1), R[k,l-1])
      theta_inv[[k]][[l]][active[[k]]] = rgamma(length(active[[k]]), 
                                                shape = param$a_theta + 0.5*n[k], rate =1)*
        1./(param$b_theta + 0.5 * colSums(u[[k]][[l]][,active[[k]],drop=F]^2))
    }
    
    ## update latent multi-rank
    R[,l] = R[,l-1]
    
    if(l %% 500==0) {
      cat(paste0("iteration: ", l, "\n"))
    }
    
    ###--------### adaptation ###--------###
    if(l >= param$start_adapt){
      for(k in 1:n_mode){
        if (unif_rvs[l]<=exp(param$c0_adapt + param$c1_adapt*l)){
          if (Rstar[k,l] < R[k,l-1]-1 & Rstar[k,l]>0){
            # set truncation to Rstar[k,l] and  keep only active columns
            R[k,l] = Rstar[k,l] + 1
            nrUpdate = rhoNuUpdateRstar(R[,l], param, rho[[l]], k, active[[k]], nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateRstar(param, k, R[,l], g[[l]], active[[k]])
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]][active[[k]]], 
                                    param$theta_inf[k]^(-1))
            omega[[k]][[l]] = c(omega[[k]][[l]][active[[k]]],
                                1 - sum(omega[[k]][[l]][active[[k]]]))
            u[[k]][[l]] = cbind(u[[k]][[l]][,active[[k]]],
                                rnorm(n[k], mean = 0, sd = sqrt(param$theta_inf[k])))
          }
          else if (R[k,l-1] < R[k,1]) {
            # increase truncation by 1 and extend all variables
            R[k,l] = R[k,l] + 1
            nrUpdate = rhoNuUpdateR(R[,l], param, rho[[l]], k, nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateR(param, k, R[,l], g[[l]])
            v[[k]][R[k,l-1]] = rbeta(1, shape1 = 1, shape2 = param$alpha[k])
            v[[k]] = c(v[[k]], 1)
            omega[[k]][[l]] = rep(NA, R[k,l])
            omega[[k]][[l]][1] = v[[k]][1]
            for (r in 2:R[k,l]){
              omega[[k]][[l]][r] = v[[k]][r]*prod(1-v[[k]][1:(r-1)])
            }
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]], param$theta_inf[k]^(-1))
            u[[k]][[l]] = cbind(u[[k]][[l]], rnorm(n[k], mean = 0,
                                                   sd = sqrt(param$theta_inf[k])))
          }
        }
      }
    }
  }
  
  runtime = proc.time()-t0
  output = list("y" = y, "my_seed" = my_seed, "N_sampl" = N_sampl, "R" = R,
                "Rstar" = Rstar, "g" = g, "u" = u, "theta_inv" = theta_inv,
                "eta_sq_inv" = eta_sq_inv, "s_augmented" = s_augmented, "tau" = tau,
                "nu" = nu, "rho" = rho, "runtime" = runtime)
  return(output)
}



#-----------------------------# AdaTuck with NA #------------------------------#


gibbsAdaTuckNA = function(y, my_seed, N_sampl, param, Rinit){
  set.seed(my_seed)
  y = as.tensor(y)
  n_mode = y@num_modes # K
  n = y@modes # (n_1,..,n_K)
  unif_rvs = runif(N_sampl) # needed for adaptation
  
  ############### VARIABLES TO UPDATE ###################
  R = matrix(NA, nrow = n_mode, ncol = N_sampl) # numbers of factors
  Rstar = matrix(NA, nrow = n_mode, ncol = N_sampl) # number of active factors
  g = vector("list", N_sampl) # core tensor
  u = vector("list", n_mode) # factor matrices
  theta_inv = vector("list", n_mode)   # mode-specific factor matrices precisions
  omega = vector("list", n_mode) # stick-breaking weights
  s_augmented = vector("list", n_mode) # augmented data
  v = vector("list", n_mode) # stick-breaking construction beta variables
  for(k in 1:n_mode){
    u[[k]] = theta_inv[[k]] = omega[[k]] = s_augmented[[k]] = vector("list", N_sampl)
  }
  eta_sq_inv = rep(NA, N_sampl) # error variances
  tau = rep(NA, N_sampl) # global-scale variance core tensor
  nu = vector("list", n_mode) # local-scale variance core tensor
  rho = vector("list", n_mode) # core tensor shrinkage parameters
  
  ################# INITIALIZATION   ######################
  for(k in 1:n_mode){
    if(is.null(Rinit[k])){
      R[k,1] = floor((sort(n)[1] + n[k])/3)
    }
    else{
      R[k,1] = Rinit[k]
    }
    Rstar[k,1] = R[k,1] -1
    u[[k]][[1]] = matrix(rnorm(n[k]*R[k,1]), n[k], R[k,1])
    theta_inv[[k]][[1]] = rep(1, R[k,1])
    s_augmented[[k]][[1]] = rep(R[k,1], R[k,1])
    omega[[k]][[1]] = rep(1/R[k,1], R[k,1])
  }
  g[[1]] = as.tensor(array(rnorm(prod(n), sd=0.5), dim = R[,1]))
  eta_sq_inv[1] = 1
  tau[1] = 1
  nu[[1]] = rho[[1]] = array(1, dim = R[,1])
  
  ##---## NA values
  get_na = is.na(y@data) 
  nna = prod(n) - sum(get_na)
  vec_na = which(get_na)
  
  ############### GIBBS SAMPLING ########################
  t0 = proc.time()
  for (l in 2:N_sampl){
    
    for(k in 1:n_mode){
      u[[k]][[l]] = matrix(NA, n[k], R[k,l-1])
    }
    
    # sample u
    for(k in 1:n_mode){
      yk = k_unfold(y, m = k)@data
      G_k = k_unfold(g[[l-1]], m=k)@data
      fctu = eigenMapMatMult(G_k, t(kronecker_list(listmat_u(u, k, l))))
      for(i in 1:n[k]){
        na_rm = which(is.na(yk[i,]))
        fctuNA = fctu[,-na_rm]
        V_u = chol2inv(chol(diag(theta_inv[[k]][[l-1]]) + 
                              (fctuNA %*% t(fctuNA))*eta_sq_inv[[l-1]]))
        mu_u =  (V_u %*% (fctuNA %*%  yk[i,-na_rm]))*eta_sq_inv[[l-1]]
        u[[k]][[l]][i,] = mvnfast::rmvn(1, mu_u, V_u)
      }
    }
    
    # sample G 
    W = kronecker_list(list_matrix(u, l, reverse=TRUE))[-vec_na,]
    g[[l]] = sampleGNA(W, tau[l-1], nu[[l-1]], y, R[,l-1], vec_na, eta_sq_inv[l-1])
    
    # sample tau
    tau[l] = sampleTau(param, R[,l-1], g[[l]], nu[[l-1]])
    
    # sample nu
    nu[[l]] = array(GIGrvg::rgig(prod(R[,l-1]), 1/2, vec(g[[l]])^2/tau[l],
                                 c(rho[[l-1]])^2), dim = R[,l-1])
    # sample rho
    rho[[l]] = array(rgamma(prod(R[,l-1]), param$a_rho + 1, param$b_rho + 
                              abs(vec(g[[l]]))/sqrt(tau[l])), dim = R[,l-1])
    
    ## sample eta
    Y_rc = ttl(g[[l]], list_matrix(u, l, reverse=FALSE), ms = c(1,2,3))@data
    eta_sq_inv[l] = rgamma(1, shape = param$a_eta + 0.5*prod(nna), rate = param$b_eta + 
                             0.5*sum((y@data[!get_na] - Y_rc[!get_na])^2))
    
    # sample s_augmented and active variables
    active = vector("list", n_mode)
    for(k in 1:n_mode){
      s_augmented[[k]][[l]] = sampleSag(u[[k]][[l]], R[k,l-1], param, k, n, 
                                        omega[[k]][[l-1]])
      
      ## update latent multi-rank
      active[[k]] = which(s_augmented[[k]][[l]]>c(1:R[k,l-1]))
      if(length(active[[k]])<2){
        active[[k]] = c(1,2)
      }
      Rstar[k,l] = length(active[[k]])
      
      # stick breaking elements
      sampleOm = sampleOmega(s_augmented[[k]][[l]], R[k,l-1], param, k)
      omega[[k]][[l]] = sampleOm$omega
      v[[k]] = sampleOm$v
      
      ### loading precisions
      theta_inv[[k]][[l]] = rep(param$theta_inf[k]^(-1), R[k,l-1])
      theta_inv[[k]][[l]][active[[k]]] = rgamma(length(active[[k]]), 
                                                shape = param$a_theta + 0.5*n[k], 
                                                rate =1)*
        1./(param$b_theta + 0.5 * colSums(u[[k]][[l]][,active[[k]],drop=F]^2))
    }
    
    ## update latent multi-rank
    R[,l] = R[,l-1]
    
    if(l %% 500==0) {
      cat(paste0("iteration: ", l, "\n"))
    }
    
    ###--------### adaptation ###--------###
    if(l >= param$start_adapt){
      for(k in 1:n_mode){
        if (unif_rvs[l]<=exp(param$c0_adapt + param$c1_adapt*l)){
          if (Rstar[k,l] < R[k,l-1]-1 & Rstar[k,l]>0){
            # set truncation to Rstar[k,l] and  keep only active columns
            R[k,l] = Rstar[k,l] + 1
            nrUpdate = rhoNuUpdateRstar(R[,l], param, rho[[l]], k, active[[k]], nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateRstar(param, k, R[,l], g[[l]], active[[k]])
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]][active[[k]]], 
                                    param$theta_inf[k]^(-1))
            omega[[k]][[l]] = c(omega[[k]][[l]][active[[k]]],
                                1 - sum(omega[[k]][[l]][active[[k]]]))
            u[[k]][[l]] = cbind(u[[k]][[l]][,active[[k]]],
                                rnorm(n[k], mean = 0, sd = sqrt(param$theta_inf[k])))
          }
          else if (R[k,l-1] < R[k,1]) {
            # increase truncation by 1 and extend all variables
            R[k,l] = R[k,l] + 1
            nrUpdate = rhoNuUpdateR(R[,l], param, rho[[l]], k, nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateR(param, k, R[,l], g[[l]])
            v[[k]][R[k,l-1]] = rbeta(1, shape1 = 1, shape2 = param$alpha[k])
            v[[k]] = c(v[[k]], 1)
            omega[[k]][[l]] = rep(NA, R[k,l])
            omega[[k]][[l]][1] = v[[k]][1]
            for (r in 2:R[k,l]){
              omega[[k]][[l]][r] = v[[k]][r]*prod(1-v[[k]][1:(r-1)])
            }
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]], param$theta_inf[k]^(-1))
            u[[k]][[l]] = cbind(u[[k]][[l]], rnorm(n[k], mean = 0,
                                                   sd = sqrt(param$theta_inf[k])))
          }
        }
      }
    }
  }
  
  runtime = proc.time()-t0
  output = list("y" = y, "my_seed" = my_seed, "N_sampl" = N_sampl, "R" = R,
                "Rstar" = Rstar, "g" = g, "u" = u, "theta_inv" = theta_inv,
                "eta_sq_inv" = eta_sq_inv, "s_augmented" = s_augmented, "tau" = tau,
                "nu" = nu, "rho" = rho, "runtime" = runtime)
  return(output)
}


#--------------------------------# Probit AdaTuck #----------------------------#

gibbsProbitAdaTuck = function(y, my_seed, N_sampl, param, Rinit){
  set.seed(my_seed)
  y = as.tensor(y)
  n_mode = y@num_modes # K
  n = y@modes # (n_1,..,n_K)
  unif_rvs = runif(N_sampl) # needed for adaptation
  
  ############### VARIABLES TO UPDATE ###################
  R = matrix(NA, nrow = n_mode, ncol = N_sampl) # numbers of factors
  Rstar = matrix(NA, nrow = n_mode, ncol = N_sampl) # number of active factors
  g = vector("list", N_sampl) # core tensor
  u = vector("list", n_mode) # factor matrices
  theta_inv = vector("list", n_mode)   # mode-specific factor matrices precisions
  omega = vector("list", n_mode) # stick-breaking weights
  s_augmented = vector("list", n_mode) # augmented data
  v = vector("list", n_mode) # stick-breaking construction beta variables
  for(k in 1:n_mode){
    u[[k]] = theta_inv[[k]] = omega[[k]] = s_augmented[[k]] = vector("list", N_sampl)
  }
  eta_sq_inv = rep(NA, N_sampl) # error variances
  tau = rep(NA, N_sampl) # global-scale variance core tensor
  nu = vector("list", n_mode) # local-scale variance core tensor
  rho = vector("list", n_mode) # core tensor shrinkage parameters
  q = vector("list", N_sampl) # augmented data for probit model
  
  ################# INITIALIZATION   ######################
  for(k in 1:n_mode){
    if(is.null(Rinit[k])){
      R[k,1] = floor((sort(n)[1] + n[k])/3)
      if(R[k,1]<5) R[k,1] = 5
    }
    else{
      R[k,1] = Rinit[k]
    }
    Rstar[k,1] = R[k,1] -1
    u[[k]][[1]] = matrix(rnorm(n[k]*R[k,1]), n[k], R[k,1])
    theta_inv[[k]][[1]] = rep(1, R[k,1])
    s_augmented[[k]][[1]] = rep(R[k,1], R[k,1])
    omega[[k]][[1]] = rep(1/R[k,1], R[k,1])
  }
  g[[1]] = as.tensor(array(rnorm(prod(n), sd=0.5), dim = R[,1]))
  eta_sq_inv[1] = 1
  tau[1] = 1
  nu[[1]] = rho[[1]] = array(1, dim = R[,1])
  q[[1]] = as.tensor(array(0, dim = n))
  
  # define bounds for truncated normal
  Tnbound = TNcomputebound(n, y)
  
  ############### GIBBS SAMPLING ########################
  t0 = proc.time()
  for (l in 2:N_sampl){
    for(k in 1:n_mode){
      u[[k]][[l]] = matrix(NA, n[k], R[k,l-1])
    }
    # sample q
    Zm_rc = ttl(g[[l-1]], list_matrix(u, l-1, reverse=FALSE), ms = c(1,2,3))@data
    q[[l]] = sampleq(Zm_rc, Tnbound$lb, Tnbound$ub, n, y)

    ##  sample u
    for(k in 1:n_mode){
      qk = k_unfold(q[[l]], m = k)@data
      G_k = k_unfold(g[[l-1]], m=k)@data
      fctu = eigenMapMatMult(G_k, t(kronecker_list(listmat_u(u, k, l))))
      V_u = chol2inv(chol(diag(theta_inv[[k]][[l-1]]) + (fctu %*% t(fctu))))
      for(i in 1:n[k]){
        mu_u =  (V_u %*% (fctu %*%  qk[i,]))
        u[[k]][[l]][i,] = mvnfast::rmvn(1, mu_u, V_u)
      }
    }
    
    ## sample G 
    W = kronecker_list(list_matrix(u, l, reverse=TRUE))
    g[[l]] = sampleG(W, tau[l-1], nu[[l-1]], q[[l]], R[,l-1])
    
    # sample tau
    tau[l] = sampleTau(param, R[,l-1], g[[l]], nu[[l-1]])
    
    # sample nu
    nu[[l]] = array(GIGrvg::rgig(prod(R[,l-1]), 1/2, vec(g[[l]])^2/tau[l],
                                 c(rho[[l-1]])^2), dim = R[,l-1])
    
    # sample rho
    rho[[l]] = array(rgamma(prod(R[,l-1]), param$a_rho + 1, param$b_rho +
                              abs(vec(g[[l]]))/sqrt(tau[l])), dim = R[,l-1])
    
    # sample s_augmented and active variables
    active = vector("list", n_mode)
    for(k in 1:n_mode){
      s_augmented[[k]][[l]] = sampleSag(u[[k]][[l]], R[k,l-1], param, k, n, 
                                        omega[[k]][[l-1]])
      
      # update latent multi-rank
      active[[k]] = which(s_augmented[[k]][[l]]>c(1:R[k,l-1]))
      if(length(active[[k]])==0){
        active[[k]] = c(1,2)
      }
      Rstar[k,l] = length(active[[k]])
      
      # stick breaking elements
      sampleOm = sampleOmega(s_augmented[[k]][[l]], R[k,l-1], param, k)
      omega[[k]][[l]] = sampleOm$omega
      v[[k]] = sampleOm$v
      
      ### loading precisions
      theta_inv[[k]][[l]] = rep(param$theta_inf[k]^(-1), R[k,l-1])
      theta_inv[[k]][[l]][active[[k]]] = rgamma(length(active[[k]]), 
                                                shape=param$a_theta+0.5*n[k],rate = 1)*
        1./(param$b_theta+0.5*colSums(u[[k]][[l]][, active[[k]], drop=F]^2))
    }
    
    ## update latent multi-rank
    R[,l] = R[,l-1]
    if(l %% 500==0) {
      cat(paste0("iteration: ", l, "\n"))
    }
    
    ###--------### adaptation ###--------###
    if(l >= param$start_adapt){
      for(k in 1:n_mode){
        if (unif_rvs[l]<=exp(param$c0_adapt + param$c1_adapt*l)){
          if (Rstar[k,l] < R[k,l-1]-1 & Rstar[k,l]>0){
            # set truncation to Rstar[k,l] and keep only active columns
            R[k,l] = Rstar[k,l] + 1
            nrUpdate = rhoNuUpdateRstar(R[,l], param, rho[[l]], k, active[[k]], nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateRstar(param, k, R[,l], g[[l]], active[[k]])
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]][active[[k]]], 
                                    param$theta_inf[k]^(-1))
            omega[[k]][[l]] = c(omega[[k]][[l]][active[[k]]],
                                1 - sum(omega[[k]][[l]][active[[k]]]))
            u[[k]][[l]] = cbind(u[[k]][[l]][,active[[k]]], 
                                rnorm(n[k], mean=0, sd=sqrt(param$theta_inf[k])))
          }
          else if (R[k,l-1] < R[k,1]) {
            # increase truncation by 1 and extend all variables
            R[k,l] = R[k,l] + 1
            nrUpdate = rhoNuUpdateR(R[,l], param, rho[[l]], k, nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateR(param, k, R[,l], g[[l]])
            v[[k]][R[k,l-1]] = rbeta(1, shape1 = 1, shape2 = param$alpha[k])
            v[[k]] = c(v[[k]], 1)
            omega[[k]][[l]] = rep(NA, R[k,l])
            omega[[k]][[l]][1] = v[[k]][1]
            for (r in 2:R[k,l]){
              omega[[k]][[l]][r] = v[[k]][r]*prod(1-v[[k]][1:(r-1)])
            }
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]], param$theta_inf[k]^(-1))
            u[[k]][[l]] = cbind(u[[k]][[l]], rnorm(n[k], mean = 0,
                                                   sd = sqrt(param$theta_inf[k])))
          }
        }
      }
    }
  }
  
  runtime = proc.time()-t0
  output = list("y" = y, "my_seed" = my_seed, "N_sampl" = N_sampl, "R" = R,
                "Rstar" = Rstar, "g" = g, "u" = u, "theta_inv" = theta_inv,
                "s_augmented" = s_augmented, "tau" = tau,
                "nu" = nu, "rho" = rho, "runtime" = runtime)
  return(output)
}


#--------------------------# Probit AdaTuck with NA #--------------------------#

gibbsProbitAdaTuckNA = function(y, my_seed, N_sampl, param, Rinit){
  set.seed(my_seed)
  y = as.tensor(y)
  n_mode = y@num_modes # K
  n = y@modes # (n_1,..,n_K)
  unif_rvs = runif(N_sampl) # needed for adaptation
  
  ############### VARIABLES TO UPDATE ###################
  R = matrix(NA, nrow = n_mode, ncol = N_sampl) # numbers of factors
  Rstar = matrix(NA, nrow = n_mode, ncol = N_sampl) # number of active factors
  g = vector("list", N_sampl) # core tensor
  u = vector("list", n_mode) # factor matrices
  theta_inv = vector("list", n_mode)   # mode-specific factor matrices precisions
  omega = vector("list", n_mode) # stick-breaking weights
  s_augmented = vector("list", n_mode) # augmented data
  v = vector("list", n_mode) # stick-breaking construction beta variables
  for(k in 1:n_mode){
    u[[k]] = theta_inv[[k]] = omega[[k]] = s_augmented[[k]] = vector("list", N_sampl)
  }
  eta_sq_inv = rep(NA, N_sampl) # error variances
  tau = rep(NA, N_sampl) # global-scale variance core tensor
  nu = vector("list", n_mode) # local-scale variance core tensor
  rho = vector("list", n_mode) # core tensor shrinkage parameters
  q = vector("list", N_sampl) # augmented data for probit model
  
  ################# INITIALIZATION   ######################
  for(k in 1:n_mode){
    if(is.null(Rinit[k])){
      R[k,1] = floor((sort(n)[1] + n[k])/3)
      if(R[k,1]<5) R[k,1] = 5
    }
    else{
      R[k,1] = Rinit[k]
    }
    Rstar[k,1] = R[k,1] -1
    u[[k]][[1]] = matrix(rnorm(n[k]*R[k,1]), n[k], R[k,1])
    theta_inv[[k]][[1]] = rep(1, R[k,1])
    s_augmented[[k]][[1]] = rep(R[k,1], R[k,1])
    omega[[k]][[1]] = rep(1/R[k,1], R[k,1])
  }
  g[[1]] = as.tensor(array(rnorm(prod(n), sd=0.5), dim = R[,1]))
  eta_sq_inv[1] = 1
  tau[1] = 1
  nu[[1]] = rho[[1]] = array(1, dim = R[,1])
  q[[1]] = as.tensor(array(0, dim = n))
  
  ##---## NA values
  get_na = is.na(y@data) # binary tensor with NA
  nna = prod(n) - sum(get_na) # number of observed values
  vec_na = which(get_na)
  
  # define bounds for truncated rnormal
  Tnbound = TNcomputeboundNA(nna, y, vec_na)
  
  ############### GIBBS SAMPLING ########################
  t0 = proc.time()
  for (l in 2:N_sampl){
    
    for(k in 1:n_mode){
      u[[k]][[l]] = matrix(NA, n[k], R[k,l-1])
    }
    
    # sample q
    Zm_rc = ttl(g[[l-1]], list_matrix(u, l-1, reverse=FALSE), ms = c(1,2,3))@data
    q[[l]] = sampleqNA(Zm_rc, Tnbound$lb, Tnbound$ub, nna, y, vec_na, n)
    
    ##  sample u
    for(k in 1:n_mode){
      qk = k_unfold(q[[l]], m = k)@data
      G_k = k_unfold(g[[l-1]], m=k)@data
      fctu = eigenMapMatMult(G_k, t(kronecker_list(listmat_u(u, k, l))))
      for(i in 1:n[k]){
        na_rm = which(is.na(qk[i,]))
        fctuNA = fctu[,-na_rm]
        V_u = chol2inv(chol(diag(theta_inv[[k]][[l-1]]) + (fctuNA %*% t(fctuNA))))
        mu_u =  (V_u %*% (fctuNA %*%  qk[i,-na_rm]))
        u[[k]][[l]][i,] = mvnfast::rmvn(1, mu_u, V_u)
      }
    }
    
    # sample G 
    W = kronecker_list(list_matrix(u, l, reverse=TRUE))[-vec_na,]
    g[[l]] = sampleGNA(W, tau[l-1], nu[[l-1]], q[[l]], R[,l-1], vec_na)
    
    # sample tau
    tau[l] = sampleTau(param, R[,l-1], g[[l]], nu[[l-1]])
    
    # sample nu
    nu[[l]] = array(GIGrvg::rgig(prod(R[,l-1]), 1/2, vec(g[[l]])^2/tau[l],
                                 c(rho[[l-1]])^2), dim = R[,l-1])
    
    # sample rho
    rho[[l]] = array(rgamma(prod(R[,l-1]), param$a_rho + 1, param$b_rho +
                              abs(vec(g[[l]]))/sqrt(tau[l])), dim = R[,l-1])
    
    # sample s_augmented and active variables
    active = vector("list", n_mode)
    for(k in 1:n_mode){
      s_augmented[[k]][[l]] = sampleSag(u[[k]][[l]], R[k,l-1], param, k, n, 
                                        omega[[k]][[l-1]])
      
      # update latent multi-rank
      active[[k]] = which(s_augmented[[k]][[l]]>c(1:R[k,l-1]))
      if(length(active[[k]])==0){
        active[[k]] = c(1,2)
      }
      Rstar[k,l] = length(active[[k]])
      
      # stick breaking elements
      sampleOm = sampleOmega(s_augmented[[k]][[l]], R[k,l-1], param, k)
      omega[[k]][[l]] = sampleOm$omega
      v[[k]] = sampleOm$v
      
      # loading precisions
      theta_inv[[k]][[l]] = rep(param$theta_inf[k]^(-1), R[k,l-1])
      theta_inv[[k]][[l]][active[[k]]] = rgamma(length(active[[k]]), shape = param$a_theta 
                                                + 0.5*n[k], rate = 1)*
        1./(param$b_theta + 0.5*colSums(u[[k]][[l]][, active[[k]], drop=F]^2))
    }
    
    ## update latent multi-rank
    R[,l] = R[,l-1]
    if(l %% 500==0) {
      cat(paste0("iteration: ", l, "\n"))
    }
    
    ###--------### adaptation ###--------###
    if(l >= param$start_adapt){
      for(k in 1:n_mode){
        if (unif_rvs[l]<=exp(param$c0_adapt + param$c1_adapt*l)){
          if (Rstar[k,l] < R[k,l-1]-1 & Rstar[k,l]>0){
            # set truncation to Rstar[k,l] and keep only active columns
            R[k,l] = Rstar[k,l] + 1
            nrUpdate = rhoNuUpdateRstar(R[,l], param, rho[[l]], k, active[[k]], nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateRstar(param, k, R[,l], g[[l]], active[[k]])
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]][active[[k]]], param$theta_inf[k]^(-1))
            omega[[k]][[l]] = c(omega[[k]][[l]][active[[k]]],
                                1 - sum(omega[[k]][[l]][active[[k]]]))
            u[[k]][[l]] = cbind(u[[k]][[l]][,active[[k]]],
                                rnorm(n[k], mean = 0, sd = sqrt(param$theta_inf[k])))
          }
          else if (R[k,l-1] < R[k,1]) {
            # increase truncation by 1 and extend all variables
            R[k,l] = R[k,l] + 1
            nrUpdate = rhoNuUpdateR(R[,l], param, rho[[l]], k, nu[[l]])
            rho[[l]] = nrUpdate$rho
            nu[[l]] = nrUpdate$nu
            g[[l]] = gUpdateR(param, k, R[,l], g[[l]])
            v[[k]][R[k,l-1]] = rbeta(1, shape1 = 1, shape2 = param$alpha[k])
            v[[k]] = c(v[[k]], 1)
            omega[[k]][[l]] = rep(NA, R[k,l])
            omega[[k]][[l]][1] = v[[k]][1]
            for (r in 2:R[k,l]){
              omega[[k]][[l]][r] = v[[k]][r]*prod(1-v[[k]][1:(r-1)])
            }
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]], param$theta_inf[k]^(-1))
            u[[k]][[l]] = cbind(u[[k]][[l]], rnorm(n[k], mean = 0,
                                                   sd = sqrt(param$theta_inf[k])))
          }
        }
      }
    }
  }
  
  runtime = proc.time()-t0
  output = list("y" = y, "my_seed" = my_seed, "N_sampl" = N_sampl, "R" = R,
                "Rstar" = Rstar, "g" = g, "u" = u, "theta_inv" = theta_inv,
                "s_augmented" = s_augmented, "tau" = tau,
                "nu" = nu, "rho" = rho, "runtime" = runtime)
  return(output)
}

