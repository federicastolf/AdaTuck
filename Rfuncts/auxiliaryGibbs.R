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


###############################################################################
#####-----------##### functions for Gibbs sampler updates ####-------------####

TNcomputebound = function(n, y){
  ubtn = rep(0, prod(n)) # upper bounds
  lbtn = rep(-Inf, prod(n)) # lower bounds
  yvecobs = vec(y)
  for(i in 1:length(yvecobs)){
    if(yvecobs[i]==1){
      ubtn[i] = Inf
      lbtn[i] = 0
    }
  }
  return(list("lb" = lbtn, "ub" = ubtn))
}

TNcomputeboundNA = function(nmiss, ydata, NAidx){
  ubtn = rep(0, nmiss) # upper bounds
  lbtn = rep(-Inf, nmiss) # lower bounds
  yvecobs = vec(ydata)[-NAidx]
  for(i in 1:length(yvecobs)){
    if(yvecobs[i]==1){
      ubtn[i] = Inf
      lbtn[i] = 0
    }
  }
  return(list("lb" = lbtn, "ub" = ubtn))
}

sampleq = function(Zm_rc, lbtn, ubtn, n, ydata){
  zmean = vec(as.tensor(Zm_rc))
  qvec = rtruncnorm(prod(n), a = lbtn, b = ubtn, mean = zmean,
                    sd = rep(1, prod(n)))
  qallv = vec(ydata)
  qallv[!is.na(qallv)] = qvec
  return(as.tensor(array(qallv, dim=n)))
}

sampleqNA = function(Zm_rc, lbtn, ubtn, nmiss, ydata, idxNA, ndim){
  zmean = vec(as.tensor(Zm_rc))[-idxNA]
  qvec = rtruncnorm(nmiss, a = lbtn, b = ubtn, mean = zmean, sd = rep(1, nmiss))
  qallv = vec(ydata)
  qallv[!is.na(qallv)] = qvec
  return(as.tensor(array(qallv, dim = ndim)))
}

sampleG = function(W, tau, nu, q, Rdim, eta_sq_inv = 1){
  Vg = diag(1/(tau*c(nu))) +  (eigenMapMatMult(t(W), W))*eta_sq_inv
  var_g = chol2inv(chol(Vg))
  mu_g = (var_g %*% (t(W) %*% vec(q)))*eta_sq_inv
  vecg = mvnfast::rmvn(1, mu_g, var_g)
  return(as.tensor(array(vecg, dim = Rdim)))
}

sampleTau = function(param, Rdim, g, nu){
  ptau_gig = param$a_tau - prod(Rdim)/2
  chi_gig = sum(g@data^2 /nu)
  if(ptau_gig < 0 & chi_gig == 0) chi_gig = 0.0001
  taus = GIGrvg::rgig(1, ptau_gig, chi_gig, 2*param$b_tau)
  return(taus)
}

sampleSag = function(uk, Rdim, param, k, n, omega){
  vec_L = colSums(uk[,1:Rdim, drop=F]^2)
  lonP_N = matrix(-0.5*vec_L/param$theta_inf[k] - 
                    rep(0.5*n[k]*log(2*pi*param$theta_inf[k]), Rdim), Rdim, Rdim)
  logP_T = matrix(-(0.5*n[k]+param$a_theta)*log(1+0.5*vec_L/param$b_theta) +
                    rep(lgamma(0.5*n[k]+param$a_theta) - lgamma(param$a_theta) - 
                          0.5*n[k]*log(2*pi*param$b_theta), Rdim), Rdim, Rdim)
  lonP_Z = lonP_N
  lonP_Z[upper.tri(lonP_Z,diag = F)] = logP_T[upper.tri(lonP_Z, diag = F)]
  lonP_Z = lonP_Z + t(matrix(log(omega), Rdim, Rdim))
  max_pZ = matrix(apply(lonP_Z, 1, max), Rdim, Rdim)
  pr_z =  exp(lonP_Z - max_pZ)
  pr_Tot = apply(pr_z, 1, sum)
  pr_z =  pr_z/pr_Tot
  sAg = as.vector(Hmisc::rMultinom(pr_z,1)) 
  return(sAg)
} 

sampleOmega = function(s_augmented, Rdim, param, k){
  count_eq = colSums(as.matrix(s_augmented == t(c(1:Rdim)*matrix(1, Rdim, Rdim))))
  count_gr = rev(cumsum(rev(c(count_eq[-1], 0))))
  v = c(rbeta(Rdim-1, shape1 = 1+count_eq[-Rdim], shape2 = param$alpha[k] 
              + count_gr[-Rdim]), 1.)
  omega = v*c(1, cumprod(1-v[-Rdim]))
  return(list("v" = v, "omega" = omega))
}

rhoNuUpdateRstar = function(Rl, param, rhoT, k, activek, nuT){
  rhosample = rgamma(prod(Rl[-k]), param$a_rho, param$b_rho)
  rho_tmp = rbind(k_unfold(as.tensor(rhoT), m=k)@data[activek,], rhosample)
  rs = k_fold(rho_tmp, m = k, modes = Rl)@data
  nusampl = rexp(prod(Rl[-k]), rhosample^2/2)
  nu_tmp = rbind(k_unfold(as.tensor(nuT), m=k)@data[activek,], nusampl)
  ns = k_fold(nu_tmp, m = k, modes = Rl)@data
  return(list("rho" = rs, "nu" = ns))
} 

gUpdateRstar = function(param, k, Rl, gT, activek){
  gsampl = rep(param$theta_inf[k], prod(Rl[-k]))
  gtmp = rbind(k_unfold(gT, m=k)@data[activek,], gsampl)
  gs = k_fold(gtmp, m = k, modes = Rl)
  return(gs)
}

gUpdateR = function(param, k, Rl, gT){
  gsampl = rep(param$theta_inf[k], prod(Rl[-k]))
  gtmp = rbind(k_unfold(gT, m=k)@data, gsampl)
  gs = k_fold(gtmp, m = k, modes = Rl)
  return(gs)
}

rhoNuUpdateR = function(Rl, param, rhoT, k,  nuT){
  rhosample = rgamma(prod(Rl[-k]), param$a_rho, param$b_rho)
  rho_tmp = rbind(k_unfold(as.tensor(rhoT), m=k)@data, rhosample)
  rs = k_fold(rho_tmp, m = k, modes = Rl)@data
  nusampl = rexp(prod(Rl[-k]), rhosample^2/2)
  nu_tmp = rbind(k_unfold(as.tensor(nuT), m=k)@data, nusampl)
  ns = k_fold(nu_tmp, m = k, modes = Rl)@data
  return(list("rho" = rs, "nu" = ns))
}

sampleGNA = function(W, tau, nu, q, Rdim, vec_na, eta_sq_inv = 1){
  Vg = diag(1/(tau*c(nu))) +  (eigenMapMatMult(t(W), W))*eta_sq_inv
  var_g = chol2inv(chol(Vg))
  mu_g = (var_g %*% (t(W) %*% vec(q)[-vec_na]))*eta_sq_inv
  vecg = mvnfast::rmvn(1, mu_g, var_g)
  return(as.tensor(array(vecg, dim = Rdim)))
}

