library(dplyr)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

#bring in mcmc chains
setwd("~/Desktop/Rscripts/Data/Fish")

chainsA <- readRDS("Joint_GW_full_A2_2.RDA")
chainsB <- readRDS("Joint_GW_full_B2_2.RDA")
chainsC <- readRDS("Joint_GW_full_C2_2.RDA")

#remove chains that did not converge
chainsAconverged = chainsA[-c(1:3)]
chainsBconverged = chainsB[-c(2)]
chainsCconverged = chainsC[-c(1,6,7)]

fishchains = c(chainsAconverged,chainsBconverged,chainsCconverged)

# #get best fit
# get_best_fit = function(chain.list){
#   L = length(chain.list)
#   chain.scores = numeric()
#   for(i in 1:L){
#     chain.scores[i] = max(chain.list[[i]]$log.p)
#   }
#   list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
#        chain.list[[which.max(chain.scores)]]$cov.jump,
#        max(chain.list[[which.max(chain.scores)]]$log.p))
#   
# }
# 
# samps = get_best_fit(fishchains)
# parameters = samps[[1]]
# variances = samps[[2]]

f = numeric()
f_J = numeric()
f_N = numeric()
h = numeric()
h_J = numeric()
h_N = numeric()

for(i in 1:length(fishchains)){
  
  f = c(f,fishchains[[i]]$samples[50001:250000,"f"]) 
  f_J = c(f_J,fishchains[[i]]$samples[50001:250000,"f"]*fishchains[[i]]$samples[50001:250000,"f_J"]) 
  f_N = c(f_N,fishchains[[i]]$samples[50001:250000,"f"]*fishchains[[i]]$samples[50001:250000,"f_N"]) 
  h = c(h,fishchains[[i]]$samples[50001:250000,"h"]) 
  h_J = c(h_J,fishchains[[i]]$samples[50001:250000,"h"]*fishchains[[i]]$samples[50001:250000,"h_J"]) 
  h_N = c(h_N,fishchains[[i]]$samples[50001:250000,"h"]*fishchains[[i]]$samples[50001:250000,"h_N"]) 
  }

#set up correct SE calculation for mcmc chains 
library(coda)
SEM_mcmc <- function(x){
  ess <- effectiveSize(mcmc(x))
  sd(x) / sqrt(ess)
}

#####attack####
attackratesdf = data.frame(f,f_J,f_N)
colnames(attackratesdf) = c("A","J","N")

#convert from relative attack rates to total
meanattack = data.frame(as.list(colMeans(attackratesdf))) 

meanattack

meanattacklong = meanattack %>% pivot_longer(cols = everything(),names_to = "stage", values_to = "mean")

seattack = data.frame(as.list(apply(attackratesdf, 2, SEM_mcmc)))

seattacklong = seattack %>% pivot_longer(cols = everything(),names_to = "stage", values_to = "SE")

attack_Mean_SE = merge(seattacklong,meanattacklong,by="stage")

attack = ggplot(data=attack_Mean_SE,aes(x=stage,y=mean)) + geom_point(size = 4) + 
  geom_errorbar(data = attack_Mean_SE,
                aes(x = stage, ymin = mean - SE, ymax = mean + SE, width = 0.05)) + theme_classic() +
  labs(x = "Copepod Stage", y = expression('Average Attack Rate, Day, L' ^ -1),  ) + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

attack

####handle####
handlingtimedf = data.frame(h,h_J,h_N)
colnames(handlingtimedf) = c("A","J","N")

#convert from relative attack rates to total
meanhandle = data.frame(as.list(colMeans(handlingtimedf))) 

meanhandle

meanhandlelong = meanhandle %>% pivot_longer(cols = everything(),names_to = "stage", values_to = "mean")

sehandle = data.frame(as.list(apply(handlingtimedf, 2, SEM_mcmc)))

sehandlelong = sehandle %>% pivot_longer(cols = everything(),names_to = "stage", values_to = "SE")

handle_Mean_SE = merge(sehandlelong,meanhandlelong,by="stage")

handling = ggplot(data=handle_Mean_SE,aes(x=stage,y=mean)) + geom_point(size = 4) + 
  geom_errorbar(data = handle_Mean_SE,
                aes(x = stage, ymin = mean - SE, ymax = mean + SE, width = 0.05)) + theme_classic() +
  labs(x = "Copepod Stage", y = expression('Average Handling Time/Day x'*10^-5*'')) + 
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )  + scale_y_continuous(
    limits = c(0, 3e-5),
    breaks = seq(0, 3e-5, by = 1e-5),
    labels = function(x) x / 1e-5
  )
  
handling

ggarrange(attack,handling,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          labels = c("A", "B"),
          legend = "none",
          label.x = 0.2,
          label.y = 1,
          font.label = list(size=20) ,
          align = "v") 

library(cowplot)
plot_grid(attack+theme(plot.margin = unit(c(1, 1, 2, 2), "cm")), handling+theme(plot.margin = unit(c(1, 1, 1, 2), "cm")), labels = c('A', 'B'), label_size = 20,ncol=1)
