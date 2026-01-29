# Loading packages
library(Matrix)
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(deSolve)
library(coda)
library(ggpubr)

#####ODE####

#The ODE has copepods that transition through N,J, A stages and adults copepods (A) that can become exposed and then infected.
#All copepod stages can also be eaten by fish, and if the copepod is infectious, it can lead to be infectious fish. 

Pond_ODE =function(t, y, parameters) {
  
  with(as.list(parameters),{
    N=y[1]; J=y[2]; A=y[3]; Es = y[4:(4+latent_stages - 1)]; I = y[4+latent_stages]; L3F = y[5+latent_stages] #Preds = y[5+latent_stages]; 
    VOL = 1
    
    Preds = Preds #ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
    Pred_A = f*(Preds/VOL)/(1 + f*h*(A + sum(Es) + I +f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + sum(Es) + I + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + sum(Es) + I + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A + sum(Es) + I)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es) + I)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A + sum(Es) + I)) #density dependence in maturation
    
    dNdt = b_M*A/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A)) - (m_N_c+d_N_c)*N - cann*(A + sum(Es) + I)*N - Pred_N*N 
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J - Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A - Pred_A*A - lambda*A
    
    #presence of adults in initial conditions automatically seeds the exposed class (Es) through lambda*A
    # development of all stages
    latent_progression = latent_rate*Es
    # lost to next stage   #death      #gained from last stage
    dEsdt = -latent_progression - d_A_c*Es + c(lambda*A, latent_progression[1:(latent_stages - 1)]) - Pred_A*Es
    
    dIdt = as.numeric(latent_progression[latent_stages]) - d_A_c*I - Pred_A*I
    
    #dPredsdt = 0 #convEff*(Pred_N*N + Pred_J*J + Pred_A*A + Pred_A*sum(Es)) - d_F*Preds 
    dL3Fdt <- Pred_A*I - d_W*L3F - d_F*L3F
    #dL3Fdt <- ifelse(Preds>0,Pred_A*I - d_W*L3F - d_F*L3F,0)  
    
    result = c(dNdt,dJdt,dAdt, dEsdt, dIdt,dL3Fdt) #dPredsdt
    
    
    return(list(result))
  }
  )
}

#bring in 24 mcmc chains
setwd("~/Desktop/Rscripts/Data/Fish")
chainsA <- readRDS("Joint_GW_full_A2_2.RDA")
chainsB <- readRDS("Joint_GW_full_B2_2.RDA")
chainsC <- readRDS("Joint_GW_full_C2_2.RDA")

#remove chains that did not converge
chainsAconverged = chainsA[-c(1:3)]
chainsBconverged = chainsB[-c(2)]
chainsCconverged = chainsC[-c(1,6,7)]

#combine converged chains
fishchains = c(chainsAconverged,chainsBconverged,chainsCconverged)

#get best fit
get_best_fit = function(chain.list){
  L = length(chain.list)
  chain.scores = numeric()
  for(i in 1:L){
    chain.scores[i] = max(chain.list[[i]]$log.p)
  }
  list(chain.list[[which.max(chain.scores)]]$samples[which.max(chain.list[[which.max(chain.scores)]]$log.p),],
       chain.list[[which.max(chain.scores)]]$cov.jump,
       max(chain.list[[which.max(chain.scores)]]$log.p))
  
}

pomp_parameter_trans = function(x){
  log_par_names = c("comp_d", "comp_b", "comp_m", "cann")
  x[log_par_names] = exp(x[log_par_names])
  x
}
pomp_parameter_backtrans = function(x){
  log_par_names = c("comp_d", "comp_b", "comp_m", "cann")
  x[log_par_names] = log(x[log_par_names])
  x
}

samps = get_best_fit(fishchains)
parameters = pomp_parameter_trans(samps[[1]])
#parameters = samps[[1]]
#variances = samps[[2]]

parameters["latent_stages"] = 60
parameters["latent_rate"] = 4.3
parameters["lambda"] = 0.032
parameters["d_W"] = 0.05
parameters["d_F"] = 0#0.05 # This should be zero, since you're holding fish constant
#parameters["convEff"] = 0.001 #how many fish can you build by eating one adult, temper for nauplii and juveniles (mass of n/mass over a)
parameters["Preds"] = 0.03
#parameters["Period"] = 100000

parameters = unlist(parameters)
Exposed_names = paste0("E", 1:parameters["latent_stages"])
Exposed_values = rep(0, times=parameters["latent_stages"])
names(Exposed_values) = Exposed_names
Exposed_values


Initial_conditions = c(N = 500, J=1000 , A=50, Exposed_values, I = 0, L3F = 0)
timespan = 365

#Run pond simulation
PondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters, hmax=1,
                         method="lsoda", func=Pond_ODE)) 

# Fitted model is black
plot(L3F ~ time, data=PondSim, typ="l")

parameters2 = parameters
parameters2[c("f_N", "f_J")] = 1
PondSim = dataparametersPondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters2, hmax=1,
                         method="lsoda", func=Pond_ODE)) 

# No attack preference is red
lines(L3F ~ time, data=PondSim, col="red")


parameters3 = parameters
parameters3[c("h_N", "h_J")] = 1
PondSim = dataparametersPondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters3, hmax=1,
                                                 method="lsoda", func=Pond_ODE)) 

lines(L3F ~ time, data=PondSim, type="l", lty=2, lwd=3)

parameters4 = parameters
parameters4[c("f_N", "f_J", "h_N", "h_J")] = 1
PondSim = dataparametersPondSim = data.frame(ode(y = Initial_conditions, times=1:timespan, parms=parameters4, hmax=1,
                                                 method="lsoda", func=Pond_ODE)) 

lines(L3F ~ time, data=PondSim, type="l", lty=2, col="red", lwd=3)

parameters



#####Pond sim over a fish density gradient#####
timespan = 365*2
preds_gradient <- seq(0, 0.13, by = 0.002)
Initial_conditions = c(N = 500, J=1000 , A=50, Exposed_values, I = 0, L3F = 0)

L3F_results = numeric()
I_results = numeric()
params = parameters
for(i in 1:length(preds_gradient)){
  params["Preds"] = preds_gradient[i]
  sim <- data.frame(
    ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = params,
      #hmax = 1,
      method = "lsoda",
      func = Pond_ODE
    ))
    L3F_results[i] <- sim[timespan, "L3F"]
    I_results[i] <- sim[timespan, "I"]
    print(i)
}

#####Pond sim no variation in f #####
timespan = 365*2
preds_gradient <- seq(0, 0.13, by = 0.002)
Initial_conditions = c(N = 500, J=1000 , A=50, Exposed_values, I = 0, L3F = 0)

L3F_results_noF = numeric()
I_results_noF = numeric()
params = parameters2
for(i in 1:length(preds_gradient)){
  params["Preds"] = preds_gradient[i]
  sim <- data.frame(
    ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = params,
      #hmax = 1,
      method = "lsoda",
      func = Pond_ODE
    ))
  L3F_results_noF[i] <- sim[timespan, "L3F"]
  I_results_noF[i] <- sim[timespan, "I"]
  print(i)
}

#####Pond sim no variation in h #####
timespan = 365*2
preds_gradient <- seq(0, 0.13, by = 0.002)
Initial_conditions = c(N = 500, J=1000 , A=50, Exposed_values, I = 0, L3F = 0)

L3F_results_noH = numeric()
I_results_noH = numeric()
params = parameters3
for(i in 1:length(preds_gradient)){
  params["Preds"] = preds_gradient[i]
  sim <- data.frame(
    ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = params,
      #hmax = 1,
      method = "lsoda",
      func = Pond_ODE
    ))
  L3F_results_noH[i] <- sim[timespan, "L3F"]
  I_results_noH[i] <- sim[timespan, "I"]
  print(i)
}

#####Pond sim no variation in h or f #####
timespan = 365*2
preds_gradient <- seq(0, 0.13, by = 0.002)
Initial_conditions = c(N = 500, J=1000 , A=50, Exposed_values, I = 0, L3F = 0)

L3F_results_noHF = numeric()
I_results_noHF = numeric()
params = parameters4
for(i in 1:length(preds_gradient)){
  params["Preds"] = preds_gradient[i]
  sim <- data.frame(
    ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = params,
      #hmax = 1,
      method = "lsoda",
      func = Pond_ODE
    ))
  L3F_results_noHF[i] <- sim[timespan, "L3F"]
  I_results_noHF[i] <- sim[timespan, "I"]
  print(i)
}

par(mfrow = c(1, 2))
plot(preds_gradient, L3F_results, type="l")
lines(preds_gradient, L3F_results_noF, col="red")
lines(preds_gradient, L3F_results_noH, lty="dashed", lwd=2)
lines(preds_gradient, L3F_results_noHF, lty="dashed", col="red", lwd=3)
plot(preds_gradient, I_results, type="l")
lines(preds_gradient, I_results_noF, col="red")
lines(preds_gradient, I_results_noH, lty="dashed", lwd=2)
lines(preds_gradient, I_results_noHF, lty="dashed", col="red", lwd=3)


#make df with preds_gradient 
predsdf = data.frame(as.list(preds_gradient))
predsdf = predsdf %>% pivot_longer(cols = everything(), names_to = "Preds")
predsdf = predsdf %>% select(2)
colnames(predsdf) = "Preds"


results_list <- list(
  L3F        = L3F_results,
  L3F_noF   = L3F_results_noF,
  L3F_noH   = L3F_results_noH,
  L3F_noHF  = L3F_results_noHF,
  I         = I_results,
  I_noF     = I_results_noF,
  I_noH     = I_results_noH,
  I_noHF    = I_results_noHF
)

for (i in names(results_list)) {
  predsdf[[i]] <- results_list[[i]]
}

predsdf <- predsdf %>%
  pivot_longer(
    cols = -Preds,
    names_to = "Model",
    values_to = "Value"
  )

L3F = predsdf %>% filter(grepl("L3F", Model)) 
InfectedCopes = predsdf %>% filter(grepl("I", Model)) 

palette <- c(
  "L3F" = "black",
  "L3F_noF"  ="red",
  "L3F_noHF" = "red" ,
  "L3F_noH" = "black"
)



L3F = ggplot(data = L3F, aes(x = Preds,y = Value, group = Model, color = Model, linetype = Model)) + geom_line(alpha = 0.5,aes(linewidth = Model)) + scale_color_manual(values = palette) +
  scale_linetype_manual(values = c(
    "L3F_noHF" = "dashed",
    "L3F_noH" = "dashed", 
    "L3F" = "solid",
    "L3F_noF" = "solid"
  )) + theme_classic() + labs(x = expression('Fish Density, L' ^ -1), y = expression('L3 in Fish (L3F), L' ^ -1)) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) + scale_linewidth_manual(values = c(
    "L3F_noHF" = 2.5,
    "L3F_noH"  = 2.5,
    "L3F"      = 1.5,
    "L3F_noF"  = 1.5
  )) 



palette2 <- c(
  "I" = "black",
  "I_noF"  ="red",
  "I_noHF" = "red" ,
  "I_noH" = "black"
)

I = ggplot(data = InfectedCopes, aes(x = Preds,y = Value, group = Model, color = Model, linetype = Model)) + geom_line(alpha = 0.5,aes(linewidth=Model)) + scale_color_manual(values = palette2) +
  scale_linetype_manual(values = c(
    "I_noHF" = "dashed",
    "I_noH" = "dashed", 
    "I" = "solid",
    "I_noF" = "solid"
  )) + theme_classic() + labs (x = expression('Fish Density, L' ^ -1), y = expression('L3 in Copepods (I), L' ^ -1))  +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) + scale_linewidth_manual(values = c(
    "I_noHF" = 2.5,
    "I_noH"  = 2.5,
    "I"      = 1.5,
    "I_noF"  = 1.5
  )) 


library(ggpubr)
ggarrange(L3F,I,
          nrow = 2, ncol = 1,
          common.legend = TRUE,
          legend = "none") 


