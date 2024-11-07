library(rTensor)
library(mvnfast)
library(GIGrvg)
library(Hmisc)
library(truncnorm)

################################################################################
##-----## Gibbs sampler with data augmentation for binary tensor data ##------##

ProbitAdaTuck = function(y, my_seed, N_sampl, param, Rinit = NULL){
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
  ubtn = rep(0, prod(n)) # upper bounds
  lbtn = rep(-Inf, prod(n)) # lower bounds
  yvecobs = vec(y)
  for(i in 1:length(yvecobs)){
    if(yvecobs[i]==1){
      ubtn[i] = Inf
      lbtn[i] = 0
    }
  }
  
  ############### GIBBS SAMPLING ########################
  t0 = proc.time()
  for (l in 2:N_sampl){
    
    for(k in 1:n_mode){
      u[[k]][[l]] = matrix(NA, n[k], R[k,l-1])
      theta_inv[[k]][[l]] = s_augmented[[k]][[l]] = omega[[k]][[l]] = rep(NA, R[k,l-1])
    }
    g[[l]] = array(NA, dim = R[,l-1])
    
    # sample q
    Zm_rc = ttl(g[[l-1]], list_matrix(u, l-1, reverse=FALSE), ms = c(1,2,3))@data
    zmean = vec(as.tensor(Zm_rc))
    qvec = rtruncnorm(prod(n), a=lbtn, b=ubtn, mean = zmean, sd = rep(1, prod(n)))
    qallv = vec(y)
    qallv[!is.na(qallv)] = qvec
    q[[l]] = as.tensor(array(qallv, dim=n))
    
    ##  sample u
    for(k in 1:n_mode){
      qk = k_unfold(q[[l]], m = k)@data
      G_k = k_unfold(g[[l-1]], m=k)@data
      fctu = G_k %*% t(kronecker_list(listmat_u(u, k, l)))
      V_u = chol2inv(chol(diag(theta_inv[[k]][[l-1]]) + (fctu %*% t(fctu))))
      for(i in 1:n[k]){
        mu_u =  (V_u %*% (fctu %*%  qk[i,]))
        u[[k]][[l]][i,] = mvnfast::rmvn(1, mu_u, V_u)
      }
    }
    
    ## sample G 
    W = kronecker_list(list_matrix(u, l, reverse=TRUE))
    Vg = diag(1/(tau[l-1]*c(nu[[l-1]]))) + t(W) %*% W
    var_g = chol2inv(chol(Vg))
    mu_g = (var_g %*% (t(W) %*% vec(q[[l]])))
    vecg = mvnfast::rmvn(1, mu_g, var_g)
    g[[l]] = as.tensor(array(vecg, dim=R[,l-1]))
    
    # sample tau
    ptau_gig = param$a_tau - prod(R[,l-1])/2
    chi_gig = sum(g[[l]]@data^2 /nu[[l-1]])
    if(ptau_gig < 0 & chi_gig == 0) chi_gig = 0.0001
    tau[l] = GIGrvg::rgig(1, ptau_gig, chi_gig, 2*param$b_tau)
    
    # sample nu
    nu[[l]] = array(GIGrvg::rgig(prod(R[,l-1]), 1/2, vec(g[[l]])^2/tau[l],
                                 c(rho[[l-1]])^2), dim = R[,l-1])
    
    # sample rho
    rho[[l]] = array(rgamma(prod(R[,l-1]), param$a_rho + 1,
                            param$b_rho + abs(vec(g[[l]]))/sqrt(tau[l])),
                     dim = R[,l-1])
    
    # sample s_augmented and active variables
    active = vector("list", n_mode)
    for(k in 1:n_mode){
      vec_L = colSums(u[[k]][[l]][,1:R[k,l-1], drop=F]^2)
      lonP_N = matrix(-0.5*vec_L/param$theta_inf[k] - 
                        rep(0.5*n[k]*log(2*pi*param$theta_inf[k]), R[k,l-1]),
                      R[k,l-1], R[k,l-1])
      logP_T = matrix(-(0.5*n[k]+param$a_theta)*log(1+0.5*vec_L/param$b_theta) +
                        rep(lgamma(0.5*n[k]+param$a_theta) - lgamma(param$a_theta) - 
                              0.5*n[k]*log(2*pi*param$b_theta), R[k,l-1]), R[k,l-1], R[k,l-1])
      lonP_Z = lonP_N
      lonP_Z[upper.tri(lonP_Z,diag = F)] = logP_T[upper.tri(lonP_Z, diag = F)]
      lonP_Z = lonP_Z + t(matrix(log(omega[[k]][[l-1]]), R[k,l-1], R[k,l-1]))
      max_pZ = matrix(apply(lonP_Z, 1, max), R[k,l-1], R[k,l-1])
      pr_z =  exp(lonP_Z - max_pZ)
      pr_Tot = apply(pr_z, 1, sum)
      pr_z =  pr_z/pr_Tot
      s_augmented[[k]][[l]] = as.vector(Hmisc::rMultinom(pr_z,1)) 
      
      # update latent multi-rank
      active[[k]] = which(s_augmented[[k]][[l]]>c(1:R[k,l-1]))
      if(length(active[[k]])==0){
        active[[k]] = c(1,2)
      }
      Rstar[k,l] = length(active[[k]])
      
      # stick breaking elements
      count_eq = colSums(as.matrix(s_augmented[[k]][[l]] == 
                                     t(c(1:R[k,l-1])*matrix(1, R[k,l-1], R[k,l-1]))))
      count_gr = rev(cumsum(rev(c(count_eq[-1], 0))))
      v[[k]] = c(rbeta(R[k,l-1]-1, shape1 = 1+count_eq[-R[k,l-1]], 
                       shape2 = param$alpha[k] + count_gr[-R[k,l-1]]), 1.)
      omega[[k]][[l]] = v[[k]] *c(1, cumprod(1-v[[k]][-R[k,l-1]]))
      
      ### loading precisions
      theta_inv[[k]][[l]] = rep(param$theta_inf[k]^(-1), R[k,l-1])
      theta_inv[[k]][[l]][active[[k]]] = rgamma(length(active[[k]]), shape = param$a_theta 
                                                + 0.5*n[k], rate = 1)*1./(param$b_theta +
                                                0.5*colSums(u[[k]][[l]][, active[[k]], drop=F]^2))
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
            # set truncation to Rstar[k,l] and subset all variables, keeping only active columns
            R[k,l] = Rstar[k,l] + 1
            rhosample = rgamma(prod(R[-k,l]), param$a_rho, param$b_rho)
            rho_tmp = rbind(k_unfold(as.tensor(rho[[l]]), m=k)@data[active[[k]],], rhosample)
            rho[[l]] = k_fold(rho_tmp, m = k, modes = R[,l])@data
            nusampl = rexp(prod(R[-k,l]), rhosample^2/2)
            nu_tmp = rbind(k_unfold(as.tensor(nu[[l]]), m=k)@data[active[[k]],], nusampl)
            nu[[l]] = k_fold(nu_tmp, m = k, modes = R[,l])@data
            gsampl = rep(param$theta_inf[k], prod(R[-k,l]))
            gtmp = rbind(k_unfold(g[[l]], m=k)@data[active[[k]],], gsampl)
            g[[l]] = k_fold(gtmp, m = k, modes = R[,l])
            theta_inv[[k]][[l]] = c(theta_inv[[k]][[l]][active[[k]]], param$theta_inf[k]^(-1))
            omega[[k]][[l]] = c(omega[[k]][[l]][active[[k]]],
                                1 - sum(omega[[k]][[l]][active[[k]]]))
            u[[k]][[l]] = cbind(u[[k]][[l]][,active[[k]]],
                                rnorm(n[k], mean = 0, sd = sqrt(param$theta_inf[k])))
          }
          else if (R[k,l-1] < R[k,1]) {
            # increase truncation by 1 and extend all variables, sampling from the prior/model
            R[k,l] = R[k,l] + 1
            rhosample = rgamma(prod(R[-k,l]), param$a_rho, param$b_rho)
            rho_tmp = rbind(k_unfold(as.tensor(rho[[l]]), m=k)@data, rhosample)
            rho[[l]] = k_fold(rho_tmp, m = k, modes = R[,l])@data
            nusampl = rexp(prod(R[-k,l]), rhosample^2/2)
            nu_tmp = rbind(k_unfold(as.tensor(nu[[l]]), m=k)@data, nusampl)
            nu[[l]] = k_fold(nu_tmp, m = k, modes = R[,l])@data
            gsampl = rep(param$theta_inf[k], prod(R[-k,l]))
            gtmp = rbind(k_unfold(g[[l]], m=k)@data, gsampl)
            g[[l]] = k_fold(gtmp, m = k, modes = R[,l])
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


###############################################################################
#####--------##### functions for tensor algebraic operations ####----------####


list_matrix = function(listm,  sample_idx, reverse=FALSE){
  "------------------------- --------------------------------------------------
  It returns the list of matrix for sample_idx-iteration 
  if reverse=T it reverse the ordere of the element in the list
  listm: list of matrix
  sample_idx: index in the gibbs of the elements for listm (can be vector or integer)
  ----------------------------------------------------------------------------"
  n_mode = length(listm)
  if(length(sample_idx)==1) sample_idx = rep(sample_idx, n_mode)
  m_it = vector("list", n_mode)
  for(i in 1:n_mode){
    m_it[[i]] = listm[[i]][[sample_idx[i]]]
  }
  m_it[sapply(m_it, is.null)] <- NULL
  if(reverse) m_it = rev(m_it)
  return(m_it)
}

listmat_u = function(listm, idx_u, idx_l, all_samel=F){
  n_mode = length(listm)
  m_it = vector("list", (n_mode-1))
  idx_dim = c(n_mode:1)
  idx_dim = idx_dim[ !idx_dim  %in% idx_u]
  idx_gibbs = rep((idx_l-1), (n_mode-1))
  if(!all_samel){
    if(idx_u>1){
      idx_gibbs[length(idx_gibbs) - ((idx_u-2):0)] = idx_l 
    }
  }
  for(i in 1:length(m_it)){
    m_it[[i]] = listm[[idx_dim[i]]][[idx_gibbs[i]]]
  }
  m_it[sapply(m_it, is.null)] <- NULL
  return(m_it)
}
