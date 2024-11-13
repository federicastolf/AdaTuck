library(pheatmap)
library(tidyverse)
library(gridExtra)

rm(list=ls())
source("Rfuncts/AdaTuckgibbs.R")
load("data/LicoriceData.Rdata")

# standardize data
Ystd = array(scale(c(Licorice)), dim(Licorice))

mseed = 345
param = list(a_eta = 1, b_eta = 0.3, a_theta = 2, b_theta=2,
             theta_inf = rep(0.05, length(Licorice)), alpha = rep(3, length(Licorice)),
             start_adapt = 100, c0_adapt = -1, c1_adapt = -5*10^(-4),
             a_tau = 2, b_tau = 2, a_rho = 10, b_rho = 10)
N_sampl = 12000
fitLicorice = AdaTuck(Ystd, mseed, N_sampl, param, Rinit = c(10, 40, 8))

## remove burnin
N_sampl = fitLicorice$N_sampl
burn_in = seq(8000+1, N_sampl, by=4)

#---------------------------# median of multi-rank #--------------------------#
apply(fitLicorice$Rstar[,burn_in],1,median)

#-------------------# plot posterior samples #-----------------------------#
Ypost_rc = matrix(NA, length(burn_in), length(Licorice))
for(i in 1:length(burn_in)){
  err_sample = array(rnorm(prod(dim(Licorice)), 0, 
                           1/sqrt(fitLicorice$eta_sq_inv[burn_in[i]])), dim = dim(Licorice))
  Y_rc = ttl(fitLicorice$g[[burn_in[i]]], list_matrix(fitLicorice$u, burn_in[i], 
                                                      reverse=FALSE), ms = c(1,2,3))@data
  Y_rc2 = Y_rc + err_sample
  Ypost_rc[i,] = c(Y_rc2)
}
Ypostmean = apply(Ypost_rc, 2, mean) # median
Ypostlb = apply(Ypost_rc, 2, quantile, probs = 0.05) # lower CI
Ypostub = apply(Ypost_rc, 2, quantile, probs = 0.95) # upper CI

# transform to tensor
Ypostmeant = array(Ypostmean, dim = dim(Licorice))
Ypostlbt = array(Ypostlb, dim = dim(Licorice))
Ypostubt = array(Ypostub, dim = dim(Licorice))

##---------------------## plot sensor response curves ##----------------------##
idgood = 15 
idbad = 2
s1 = 3
s2 = 11
s3 = 5
idgood2 = 12 
idbad2 = 7
dm = cbind.data.frame(c(Ypostmeant[idgood,,s1], Ypostmeant[idbad,,s1], Ypostmeant[idgood,,s2], 
                        Ypostmeant[idbad,,s2], Ypostmeant[idgood,,s3], Ypostmeant[idbad,,s3], 
                        Ypostmeant[idgood2,,s1], Ypostmeant[idbad2,,s1], Ypostmeant[idgood2,,s2],
                        Ypostmeant[idbad2,,s2], Ypostmeant[idgood2,,s3], Ypostmeant[idbad2,,s3]),
                      c(Ypostlbt[idgood,,s1], Ypostlbt[idbad,,s1], Ypostlbt[idgood,,s2],
                        Ypostlbt[idbad,,s2], Ypostlbt[idgood,,s3], Ypostlbt[idbad,,s3],
                        Ypostlbt[idgood2,,s1], Ypostlbt[idbad2,,s1], Ypostlbt[idgood2,,s2],
                        Ypostlbt[idbad2,,s2], Ypostlbt[idgood2,,s3], Ypostlbt[idbad2,,s3]),
                      c(Ypostubt[idgood,,s1], Ypostubt[idbad,,s1], Ypostubt[idgood,,s2],
                        Ypostubt[idbad,,s2], Ypostubt[idgood,,s3], Ypostubt[idbad,,s3],
                        Ypostubt[idgood2,,s1], Ypostubt[idbad2,,s1], Ypostubt[idgood2,,s2],
                        Ypostubt[idbad2,,s2], Ypostubt[idgood2,,s3], Ypostubt[idbad2,,s3]),
                      rep(c(rep("Good",241),rep("Bad", 241)), 6), 
                      c(rep(c(rep("good",241),rep("bad", 241)), 3),
                        rep(c(rep("good1",241),rep("bad1", 241)), 3)),
                      rep(c(rep("Sensor 3", 482), rep("Sensor 11", 482), 
                            rep("Sensor 5", 482)),2), rep(c(1:241),12))
colnames(dm) = c("mean", "lb","ub","Quality", "quality2", "sensor", "time")

pSig = ggplot(data = dm) +
  geom_line(aes(x = time, y = lb, color = Quality, group = quality2), linetype = "dotted") +
  geom_line(aes(x = time, y = ub, color = Quality, group = quality2), linetype = "dotted") +
  geom_ribbon(aes(x = time, ymax = ub, ymin = lb, fill = Quality, group = quality2), 
              alpha = 0.3) +
  geom_line(aes(x = time, y = mean, color = Quality, group = quality2), size=1) +
  facet_wrap(~sensor, scales = "free") +
  xlab("Time") + ylab("Sensor signal") +
  theme_light() + 
  theme(legend.position = "top", strip.text = element_text(size=13),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
pSig
# ggsave(plot = pSig, filename = "postSignal.png",  width = 9, height = 3.5)

#-------------------# plot U^(3) #-----------------------------#
m1 = fitLicorice$u[[3]][[12000]][,c(1:7)]

paletteLength = 50
myColor = colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks = c(seq(min(m1), 0, length.out=ceiling(paletteLength/2) + 1), 
             seq(max(m1)/paletteLength, max(m1), 
                 length.out=floor(paletteLength/2)))

phmap1 = pheatmap(m1, cluster_rows = F, cluster_cols = F, color=myColor,
                  breaks=myBreaks)


#--------------------# plot all sensors for one sample #-----------------------#
ms = ls = us = sn = NULL
idAlls = 15
for(i in 1:dim(Licorice)[3]){
  ms = c(ms, Ypostmeant[idAlls,,i])
  ls = c(ls, Ypostlbt[idAlls,,i])
  us = c(us, Ypostubt[idAlls,,i])
  sn = c(sn, rep(i,dim(Licorice)[2]))
}
dallSensor = cbind.data.frame(ms,ls, us, sn)
dallSensor$sn = as.factor(dallSensor$sn)
dallSensor$time = rep(c(1:dim(Licorice)[2]), dim(Licorice)[3])

dallSensor1 = dallSensor %>% filter(sn!="3")
pSensor = ggplot(data = dallSensor) +
  geom_line(aes(x = time, y = ms, group = sn), size=0.5) +
  xlab("Time") + ylab("Sensor signal") +
  theme_light() + 
  theme(legend.position = "top", strip.text = element_text(size=13), 
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
pSensor

map_all = grid.arrange(grobs = list(phmap1[[4]], pSensor), nrow=1)

# ggsave(plot = map_all, filename = "postSensro.png",  width = 10, height = 5)
