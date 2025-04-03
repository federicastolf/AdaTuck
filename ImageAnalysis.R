rm(list = ls())

library(tidyverse)
library(tensorBF)
# BiocManager::install(c('EBImage'))
library(EBImage)
source("Rfuncts/AdaTuckgibbs.R")
source("Rfuncts/funcAdaTuck.R")

load("data/imageData.Rdata")

IMCat = as.Image(imData)
colorMode(IMCat)=2
display(IMCat, method="raster")

##-------------------------## missing data 50% ##-----------------------------##
mseed = 345
missY = 0.5
set.seed(mseed)
idx_na50 = sample(x = 1:prod(dim(imData)), size = missY*prod(dim(imData)), replace=F)
yNA50 = imData
yNA50[idx_na50] = NA

#-----# AdaTuck #-----#
param = list(a_eta = 1, b_eta = 0.3, a_theta = 2, b_theta=2,
             theta_inf = rep(0.05, length(imData)), alpha = rep(3, length(imData)),
             start_adapt = 100, c0_adapt = -1, c1_adapt = -5*10^(-4),
             a_tau = 2, b_tau = 2, a_rho = 10, b_rho = 10)
N_sampl = 12000
fitIm50 = AdaTuck(yNA50, mseed, N_sampl, param, Rinit = c(30, 35, 3))

burn_in = seq(8000+1, N_sampl)
MSEIm50 = array(0, length(burn_in))
Yrec50 = matrix(NA, length(burn_in), length(idx_na50))
for(i in 1:length(burn_in)){
  Y_rc = ttl(fitIm50$g[[burn_in[i]]], 
             list_matrix(fitIm50$u, burn_in[i], reverse=FALSE), ms = c(1,2,3))@data
  MSEIm50[i] = mean((Y_rc[idx_na50] - imData[idx_na50])^2)
  Yrec50[i,] = Y_rc[idx_na50]
}

#---------# image reconstruction #---------#
Yrec50mean = apply(Yrec50,2,mean)

YrecIm50 = imData
YrecIm50[idx_na50] = Yrec50mean
YrecIm50[YrecIm50<0] = 0
YrecIm50[YrecIm50>1] = 1
Imrec50 = as.Image(YrecIm50)
colorMode(Imrec50)=2
display(Imrec50, method="raster")


#-----# TensorBF #------#
opts <- getDefaultOpts();
opts$iter.burnin = 8000
opts$iter.sampling = 4000
opts$verbose = 0

fitImBF50 = tensorBF(yNA50, opts = opts)

nsamples50 = length(fitImBF50$posterior$X)
MSE50BF = rep(0, nsamples50)
for(i in 1:nsamples50){
  YRec = rec_tensorBF(fitImBF50$posterior, i)
  MSE50BF[i] = mean((YRec[idx_na50] - imData[idx_na50])^2)
}
summary(MSE50BF)


##-------------------------## missing data 70% ##-----------------------------##
mseed = 345
missY = 0.5
set.seed(mseed)
idx_na70 = sample(x = 1:prod(dim(imData)), size = missY*prod(dim(imData)), replace=F)
yNA70 = imData
yNA70[idx_na70] = NA

#-----# AdaTuck #-----#

fitIm70 = AdaTuck(yNA70, mseed, N_sampl, param, Rinit = c(30, 35, 3))

MSEIm70 = array(0, length(burn_in))
Yrec70 = matrix(NA, length(burn_in), length(idx_na70))
for(i in 1:length(burn_in)){
  err_sample = array(rnorm(prod(dim(imData)), 0, 1/sqrt(fitIm70$eta_sq_inv[burn_in[i]])), 
                     dim = dim(imData))
  Y_rc = ttl(fitIm70$g[[burn_in[i]]], 
             list_matrix(fitIm70$u, burn_in[i], reverse=FALSE), ms = c(1,2,3))@data
  MSEIm70[i] = mean((Y_rc[idx_na70] - imData[idx_na70])^2)
  Yrec70[i,] = Y_rc[idx_na70]
}
summary(MSEIm70)


#---------# image reconstruction #---------#
Yrec70mean = apply(Yrec70,2,mean)

YrecIm70 = imData
YrecIm70[idx_na70] = Yrec70mean
YrecIm70[YrecIm70<0] = 0
YrecIm70[YrecIm70>1] = 1
Imrec70 = as.Image(YrecIm70)
colorMode(Imrec70)=2
display(Imrec70, method = "raster")



#-----# TensorBF #------#
fitImBF70 = tensorBF(yNA70, opts = opts)

nsamples70 = length(fitImBF70$posterior$X)
MSE70BF = rep(0, nsamples70)
for(i in 1:nsamples70){
  YRec = rec_tensorBF(fitImBF70$posterior, i)
  MSE70BF[i] = mean((YRec[idx_na70] - imData[idx_na70])^2)
}
summary(MSE70BF)


