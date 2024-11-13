library(pheatmap)
library(tidyverse)
library(reshape2)
library(gridExtra)

rm(list = ls())

source("Rfuncts/AdaTuckgibbs.R")
load("data/FishNorthSea.Rdata")

N_sampl = 20000
burnin = 10000
param = list(a_theta = 2, b_theta=2, theta_inf = rep(0.05, 3), 
             alpha = rep(3, 3), start_adapt = 400, c0_adapt = -1,
             c1_adapt = -5*10^(-4), a_tau = 2, b_tau = 2, a_rho = 5, b_rho = 5)
mseed = 1782

FishTuck = AdaTuck(yFish, mseed, N_sampl, param, binary=T)

## remove burnin
idx = seq(burnin, N_sampl, by=5)

## label tensor data
species_lab = row.names(yFish)
years_lab = colnames(yFish)
site_lab = dimnames(yFish)[[3]]

#---------# median of multi-rank #---------------#
apply(FishTuck$Rstar[,idx],1,median)

###------------### posterior marginal probabilities ###----------------------###
pi_rc = matrix(NA, length(idx), length(yFish))
for(i in 1:length(idx)){
  Y_rc = ttl(FishTuck$g[[idx[i]]], list_matrix(FishTuck$u, idx[i], reverse=FALSE),
             ms = c(1,2,3))@data
  pi_rc[i,] = c(pnorm(Y_rc))
}
pimean = apply(pi_rc,2,mean) # median 
pilb = apply(pi_rc,2,quantile, probs=0.05) # lower CI
piub = apply(pi_rc,2,quantile, probs=0.95) # upper CI
# transform to tensor
pimeant = array(pimean, dim=dim(yFish))
pilbt = array(pilb, dim=dim(yFish))
piubt = array(piub, dim=dim(yFish))
row.names(pimeant) = row.names(pilbt) = row.names(piubt) = species_lab
colnames(pimeant) = colnames(pilbt) =  colnames(piubt) = years_lab


#### --------------------- #### TREND plot ### ------------------------- ###
# 57 fish ---> Squalus acanthias
# c("7"=E,"6"=SE,"5"=S)
d57 = as.data.frame(pimeant[57,,])
dl57 = as.data.frame(pilbt[57,,])
du57 = as.data.frame(piubt[57,,])
colnames(d57) = colnames(dl57) = colnames(du57) = c(1:4, "S","SE","E")
d57$year = dl57$year = du57$year = rownames(d57)
d57long = melt(d57)
dl57long = melt(dl57)
du57long = melt(du57)
d57long = cbind.data.frame(d57long, dl57long$value, du57long$value)
colnames(d57long)[c(4:5)] = c("lb", "ub")
d57long$species = rep("Squalus acanthias", NROW(d57long))
d57long = d57long %>% filter(variable %in% c("E","SE","S"))
colnames(d57long)[2] = "Site"

# 3  fish --> Anarhichas lupus
#(2  = "C", 3 = "NW", 4="SW")
d3 = as.data.frame(pimeant[3,,])
dl3 = as.data.frame(pilbt[3,,])
du3 = as.data.frame(piubt[3,,])
colnames(d3) = colnames(dl3) = colnames(du3) = c("1", "C", "NW", "SW",c(5:7))
d3$year = dl3$year =du3$year =  rownames(d57)
d3long = melt(d3)
dl3long = melt(dl3)
du3long = melt(du3)
d3long = cbind.data.frame(d3long, dl3long$value, du3long$value)
colnames(d3long)[c(4:5)] = c("lb", "ub")
d3long$species = rep("Anarhichas lupus", NROW(d3long))
d3long = d3long %>% filter(variable %in% c("C","NW","SW"))
colnames(d3long)[2] = "Site"

p57 = ggplot(data = d57long) +
  geom_line(aes(x = year, y = lb, color = Site, group = Site), linetype = "dotted") +
  geom_line(aes(x = year, y = ub, color = Site, group = Site), linetype = "dotted") +
  geom_ribbon(aes(x = year, ymax = ub, ymin = lb, fill = Site, group = Site), alpha = 0.4) +
  geom_point(aes(x = year, y = value, color = Site, group = Site),size=2.5) + 
  geom_line(aes(x = year, y = value, color = Site, group = Site), size=1) +
  theme_light() + 
  ylim(c(0,1)) +
  xlab("Year") + ylab("") +
  scale_x_discrete(breaks = seq(1985,2015, by=4)) +
  ggtitle("Squalus acanthias") + ylab(expression(pi))  +
  theme(text = element_text(size = 13), legend.key.size = unit(0.9,"line"),
        plot.title = element_text(hjust = 0.5, size=11),
        strip.text = element_text(size=16), legend.position = "bottom",
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
p57

p3 = ggplot(data = d3long) +
  geom_line(aes(x = year, y = lb, color = Site, group = Site), linetype = "dotted") +
  geom_line(aes(x = year, y = ub, color = Site, group = Site), linetype = "dotted") +
  geom_ribbon(aes(x = year, ymax = ub, ymin = lb, fill = Site, group = Site), alpha = 0.4) +
  geom_point(aes(x = year, y = value, color = Site, group = Site),size=2.5) + 
  geom_line(aes(x = year, y = value, color = Site, group = Site), size=1) +
  theme_light() + 
  ylim(c(0,1)) +
  xlab("Year") + ylab(expression(pi)) +
  scale_x_discrete(breaks = seq(1985,2015, by=4)) +
  ggtitle("Anarhichas lupus") +
  theme(text = element_text(size = 13), legend.key.size = unit(0.9,"line"),
        plot.title = element_text(hjust = 0.5, size=11), 
        strip.text = element_text(size=16), legend.position = "bottom",
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
p3

plotTrend = grid.arrange(p57,p3, nrow=1)
# ggsave(plot = plotTrend, filename = "Trend.png",  width = 9, height = 3.8)


###--------------------### plot frontal slices ###-----------------------###

# load Trait data
load("data/TraitFish.Rdata")

dpim6 = as.data.frame(pimeant[,,6])
dpim6$Species = species_lab
dpim6 = inner_join(dpim6, Trait, by = "Species")
dpim6 = dpim6 %>% arrange(Biogeography)
dpim1 = as.data.frame(pimeant[,,1])
dpim1$Species = species_lab
dpim1 = inner_join(dpim1, Trait, by = "Species")
dpim1 = dpim1 %>% arrange(Biogeography)

lab_col = c("1985", rep("",4), "1990",rep("",4),"1995", rep("",4),"2000", rep("",4),
            "2005", rep("",4),"2010" ,rep("",4),"2015")

rownames(dpim1) = dpim6$Species
annotation_row = data.frame(dpim1$Biogeography)
colnames(annotation_row)[1] = "Biogeography"
rownames(annotation_row) =  dpim6$Species
lbpi = c("0","0.1","","","","0.5","","","","0.9","1")
rownames(dpim6) = dpim6$Species

pheat6 = pheatmap(dpim6[,1:31], treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
                  cluster_cols = F, border_color ="NA", legend=T, show_colnames = T, 
                  color = colorRampPalette(c("lightblue1", "navyblue"))(100),
                  labels_col = lab_col, fontsize_row = 9.6, fontsize_col = 15,
                  annotation_row = annotation_row, annotation_names_row = F, 
                  annotation_legend = F, legend_breaks = seq(0,1,0.1), legend_labels = lbpi)

pheat1 = pheatmap(dpim1[,1:31], treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
                   color=colorRampPalette(c("lightblue1", "navyblue"))(100),
                  border_color ="NA", legend=T, show_colnames = T, cluster_cols = F,
                  labels_col = lab_col, fontsize_row = 9.6, fontsize_col = 15,
                  annotation_row = annotation_row, annotation_names_row = F, 
                  annotation_legend = F, legend_breaks = seq(0,1,0.1), legend_labels = lbpi)

map_all = grid.arrange(grobs = list(pheat1[[4]], pheat6[[4]]))
# ggsave(plot = map_all, filename = "pheat.png",  width = 17, height = 15)


