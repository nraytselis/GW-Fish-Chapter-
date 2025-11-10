library(dplyr)
library(tidyverse)
library(ggpubr)
library(BayesianTools)
library(mcmcplots)
library(bayesplot)
library(ggplot2)
library(dplyr)

#####Summary####
#This code takes all 24 chains from the fish experiment, removes the burn-in period (100000), and thins the chains before plotting the 
#relative mean attack rates and relative mean handling times +/- the 95% CI. 


#bring in mcmc chains
setwd("~/Desktop/Rscripts/Data/Fish")
chainsA <- readRDS("Joint_GW_full_A2.RDA")
chainsB <- readRDS("Joint_GW_full_B2.RDA")
chainsC <- readRDS("Joint_GW_full_C2.RDA")

fishchains = c(chainsA,chainsB,chainsC)

get_samples <- function(chain.list) {
  L <- length(chain.list)
  samples_list <- vector("list", L)  # empty list to store samples
  
  for (i in 1:L) {
    samples <- chain.list[[i]]$samples  # extract 'samples' from each chain
    samples_list[[i]] <- coda::mcmc(samples)  # convert to mcmc object
  }
  
  return(samples_list)
}

samples_allchains <- get_samples(fishchains)

# remove burn in, extract log.p from each chain and convert to mcmc objects
burn_in <- 10000
thin_interval <- 10

logp_chains <- lapply(fishchains, function(chain) {
  # Extract log.p vector
  logp <- chain$log.p
  
  # Remove burn-in
  logp_burned <- logp[(burn_in + 1):length(logp)]
  
  # Thin samples
  logp_thinned <- logp_burned[seq(1, length(logp_burned), by = thin_interval)]
  
  # Convert to mcmc object
  mcmc(logp_thinned)
})

logp_mcmc_list <- mcmc.list(logp_chains)


#denplot(samples_allchains) 
#denplot(samples_allchains, parms = c("i_P")) 
#denplot(samples_allchains, parms = c("h","h_N")) 
#mcmc_dens(samples_allchains)

#plot attack rates
samples_allchains_attack <- as.data.frame(do.call(rbind, samples_allchains)) %>%
  select(f, f_N, f_J) 


samples_allchains_attacksummary = samples_allchains_attack %>% summarise(meanfJ = mean(f_J),mean_fn = mean(f_N),
                                                                         lower_ci_fj = quantile(f_J, prob=0.025),
                                                                         upper_ci_fj = quantile(f_J, prob=0.975),
                                                                         lower_ci_fn = quantile(f_N, prob=0.025),
                                                                         upper_ci_fn = quantile(f_N, prob=0.975))
mean_fs = samples_allchains_attacksummary[1:2]
low_fs = samples_allchains_attacksummary[c(3,5)]
high_fs = samples_allchains_attacksummary[c(4,6)]
Stages = c("J", "N")
plotting_fs = data.frame(Stages, "f" = as.numeric(mean_fs),
                         "f_low" = as.numeric(low_fs),
                         "f_high" = as.numeric(high_fs))

palette = c(J = "deeppink2", N = "goldenrod")
attack = ggplot(data=plotting_fs, aes(x=Stages, y=f, ymin=f_low, ymax=f_high,group=Stages,color=Stages)) + geom_linerange(linewidth=1.5) +
  geom_point(size=5) + theme_classic(base_size = 20) + labs(y = "Relative Mean Attack Rate +/- 95% CI" ) + 
  ylim(0,1) + geom_hline(yintercept = 1,linewidth = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = palette)

attack


#do the same thing for handling time 
samples_allchains_handle <- as.data.frame(do.call(rbind, samples_allchains)) %>%
  select(h, h_N, h_J) 

samples_allchains_handlesummary = samples_allchains_handle %>% summarise(mean_hj=mean(h_J),mean_hn = mean(h_N),
                                                                         lower_ci_hj = quantile(h_J, prob=0.025),
                                                                         upper_ci_hj = quantile(h_J, prob=0.975),
                                                                         lower_ci_hn = quantile(h_N, prob=0.025),
                                                                         upper_ci_hn = quantile(h_N, prob=0.975))

mean_hs = samples_allchains_handlesummary[1:2]
low_hs = samples_allchains_handlesummary[c(3,5)]
high_hs = samples_allchains_handlesummary[c(4,6)]
Stages = c("J", "N")
plotting_hs = data.frame(Stages, "h" = as.numeric(mean_hs),
                         "h_low" = as.numeric(low_hs),
                         "h_high" = as.numeric(high_hs))


palette = c(J = "deeppink2", N = "goldenrod")

handle = ggplot(data=plotting_hs, aes(x=Stages, y=h, ymin=h_low, ymax=h_high,group=Stages,color=Stages)) + geom_linerange(linewidth=1.5) +
  geom_point(size=5) + theme_classic(base_size = 20) + labs(y = "Relative Mean Handling Time +/- 95% CI" ) +
  ylim(0,1) +  geom_hline(yintercept = 1, linewidth = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = palette)

handle


#plot both together 
ggarrange(attack,handle,
          nrow = 1, ncol = 2,
          labels = c("A","B"),
          legend = "none")






