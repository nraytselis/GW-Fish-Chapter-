library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)
library(dplyr)

#####Summary#####
#This code uses the mesocosm data from the Abate pesticide experiment to make predictions for the fish experiment. Specifically, it plots the 
#the densities of copepod stages over time and across a fish density gradient. 


#ODE 
FishODE =function(t, y, parameters) { 
  
  N=y[1]; J=y[2]; A=y[3]
  
  with(as.list(parameters),{  
    
    Preds = ifelse(t<49,0,Preds) #fish added to tanks on week 4 
    
    Pred_A = f*(Preds/VOL)/(1 + f*h*(A+f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL) 
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    dNdt = b_M*A/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A)) - (m_N_c+d_N_c)*N - cann*(VOL/15)*A*N - Pred_N*N 
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J - Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A - Pred_A*A
    
    result = c(dNdt,dJdt,dAdt) 
    
    return(list(result)) 
  } 
  )  
}  

#ODE for model where there is only predation on nauplii
FishODEOnlyNPred =function(t, y, parameters) { 
  
  N=y[1]; J=y[2]; A=y[3]
  
  with(as.list(parameters),{  
    
    Preds = ifelse(t<49,0,Preds) #fish added to tanks on week 4 
    
    Pred_A = f*(Preds/VOL)/(1 + f*h*(A+f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL) 
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    dNdt = b_M*A/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A)) - (m_N_c+d_N_c)*N - cann*(VOL/15)*A*N - Pred_N*N 
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J #- Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A #- Pred_A*A
    
    result = c(dNdt,dJdt,dAdt) 
    
    return(list(result)) 
  } 
  )  
}  

#ODE for model where there is only predation on juveniles
FishODEOnlyJPred =function(t, y, parameters) { 
  
  N=y[1]; J=y[2]; A=y[3]
  
  with(as.list(parameters),{  
    
    Preds = ifelse(t<49,0,Preds) #fish added to tanks on week 4 
    
    Pred_A = f*(Preds/VOL)/(1 + f*h*(A+f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_J = f*f_J*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL)
    Pred_N = f*f_N*(Preds/VOL)/(1 + f*h*(A + f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL) 
    
    d_A_c = d_A*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_J_c = d_J*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    d_N_c = d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) #density dependence in deaths
    
    m_N_c = m_N*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    m_J_c = m_J*exp(-comp_m/VOL*(c_N*N + c_J*J + A)) #density dependence in maturation
    
    dNdt = b_M*A/2*exp(-comp_b/VOL*(c_N*N + c_J*J + A)) - (m_N_c+d_N_c)*N - cann*(VOL/15)*A*N #- Pred_N*N 
    
    dJdt = m_N_c*N - (m_J_c+d_J_c)*J - Pred_J*J
    
    dAdt = m_J_c*J - d_A_c*A #- Pred_A*A
    
    result = c(dNdt,dJdt,dAdt) 
    
    return(list(result)) 
  } 
  )  
}  

#import 32 chains from the Abate experiment 

setwd("~/Desktop/Rscripts/Data")
fullA = readRDS("Rebound_parameters_full_disp250kff_5.RDA")
fullB = readRDS("Rebound_parameters_full_disp250kff_3.RDA")
fullC = readRDS("Rebound_parameters_full_disp250kff_05.RDA")
fullD = readRDS("Rebound_parameters_full_disp250kff_01.RDA")
full <- c(fullA, fullB, fullC, fullD)

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
samps = get_best_fit(full)
pars = samps[[1]]
Initial_conditions = c(N = 1200, J=800 , A=100)
timespan = 98



######Simulations Over Time#####

#Run simulations under different handling time and attack rate scenarios. 

#A) relative h and a (f) both are same for all stages 
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=1,h_N=1,i_P=0)
parameters = c(pars,other_parameters)
sim1 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#B) h is same for all stages, a (f) is largest for for adults (visual preds) and smallest for nauplii
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=1,h_N=1,i_P=0)
parameters = c(pars,other_parameters)
sim2 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#C) h is largest for for adults and smallest for nauplii, a is same for all stages
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim3 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)
#D) h is largest for for adults and a (f) is largest for for adults (visual preds) and smallest for nauplii
other_parameters = c(Preds = 0.3, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim4 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)

#E) Control (no fish) 
other_parameters = c(Preds = 0, VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.001,h_J=0.1,h_N=0.01,i_P=0)
parameters = c(pars,other_parameters)
sim5 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)


colnames(sim1) = c("time","N1","J1","A1") #non-selective
colnames(sim2) = c("time","N2","J2","A2") #attack rates selective
colnames(sim3) = c("time","N3","J3","A3") #handling times selective 
colnames(sim4) = c("time","N4","J4","A4") #selective attack and handling
colnames(sim5) = c("time","N5","J5","A5") #control no fish
allscenarios = bind_cols(sim1,sim2,sim3,sim4,sim5)
names(allscenarios)[1] <- "Time"

Naup = allscenarios %>% select(c(1,2,6,10,14,18))
colnames(Naup) = c("Time","Non-selective","Selective_attack","Selective_handling","Selective_both","controlNofish")
Juv = allscenarios %>% select(c(1,3,7,11,15,19))
colnames(Juv) = c("Time","Non-selective","Selective_attack","Selective_handling","Selective_both","controlNofish")
Ad = allscenarios %>% select(c(1,4,8,12,16,20))
colnames(Ad) = c("Time","Non-selective","Selective_attack","Selective_handling","Selective_both","controlNofish")
NaupLong = Naup %>% pivot_longer(cols=c(2:6),names_to = "predation_scenario",values_to = "value" )
JuvLong = Juv %>% pivot_longer(cols=c(2:6),names_to = "predation_scenario",values_to = "value" )
AdLong = Ad %>% pivot_longer(cols=c(2:6),names_to = "predation_scenario",values_to = "value" )

palette <- c(
  "controlNofish"   = "grey",
  "Selective_attack" = "black", #noH
  "Non-selective"  ="red", #noHF
  "Selective_handling" = "red" , #noF
  "Selective_both" = "black" #normal 
)



a1 = ggplot(data=NaupLong,aes(x=Time,y=value,group=predation_scenario,color=predation_scenario,linetype = predation_scenario)) +
  geom_vline(xintercept = 49,color="blue",linewidth=1.5) + geom_line(linewidth=1.5) +theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + 
  guides(color = guide_legend(title = "Fish Scenario")) +
  scale_y_continuous(limits = c(0, NA)) + scale_linetype_manual(values = c("controlNofish" = "solid", 
                                   "Selective_attack" = "dashed", #noH
                                   "Non-selective" = "dashed", #noHF
                                   "Selective_handling" = "solid", #noF
                                   "Selective_both" = "solid")) #normal 
b1 = ggplot(data=JuvLong,aes(x=Time,y=value,group=predation_scenario,color=predation_scenario,linetype = predation_scenario)) +
  geom_vline(xintercept = 49,color="blue",linewidth=1.5) + geom_line(linewidth=1.5) + theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +
  scale_y_continuous(limits = c(0, NA)) + scale_linetype_manual(values = c("controlNofish" = "solid", 
                                                                             "Selective_attack" = "dashed", 
                                                                             "Non-selective" = "dashed", 
                                                                             "Selective_handling" = "solid",
                                                                             "Selective_both" = "solid"))
c1 = ggplot(data=AdLong,aes(x=Time,y=value,group=predation_scenario,color=predation_scenario,linetype = predation_scenario)) +
  geom_vline(xintercept = 49,color="blue",linewidth=1.5) + geom_line(linewidth=1.5) + theme_classic(base_size = 20)+
  labs(x="Time (days)",y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +
  scale_y_continuous(limits = c(0, NA)) + scale_linetype_manual(values = c("controlNofish" = "solid", 
                                   "Selective_attack" = "dashed", 
                                   "Non-selective" = "dashed", 
                                   "Selective_handling" = "solid",
                                   "Selective_both" = "solid"))
library(ggpubr)
ggarrange(a1,b1,c1,
          nrow = 1, ncol = 3,
          common.legend = TRUE,
          legend = "right")


####fish density gradient plots####
#predator gradient
preds_gradient <- seq(0, 0.4, by = 0.01)

#scenarios
scenario_params <- list(
  "Non-selective" = c(VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=1,h_N=1,i_P=0),
  Selective_attack = c(VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=1,h_N=1,i_P=0),
  Selective_handling = c(VOL = 1,f=1,f_J=1,f_N=1,h=0.006,h_J=0.1,h_N=0.01,i_P=0),
  Selective_both = c(VOL = 1,f=1,f_J=0.1,f_N=0.01,h=0.006,h_J=0.1,h_N=0.01,i_P=0)
)


sim_results <- list()

for (s in names(scenario_params)) {
  for (p in preds_gradient) {
    
    # update parameters
    parameters <- c(pars,
                    Preds = p,
                    scenario_params[[s]])
    
    # run simulation
    sim <- ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = parameters,
      method = "lsoda",
      func = FishODE
    )
    
    # store results
    sim_df <- as.data.frame(sim) %>%
      mutate(Preds = p,
             Scenario = s)
    
    sim_results[[paste(s, p, sep = "_")]] <- sim_df
  }
}

# bind into one dataframe
all_sims <- bind_rows(sim_results, .id = "run_id")


endpoints <- all_sims %>%
  group_by(Scenario, Preds) %>%
  slice_tail(n = 1) %>%   # last time step
  ungroup()


PredsN = ggplot(endpoints, aes(x = Preds, y = N, color = Scenario,linetype = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +
  scale_linetype_manual(values = c("controlNofish" = "dotted", 
                                     "Selective_attack" = "dashed", 
                                     "Non-selective" = "dashed", 
                                     "Selective_handling" = "solid",
                                     "Selective_both" = "solid"))
PredsJ = ggplot(endpoints, aes(x = Preds, y = J, color = Scenario,linetype = Scenario )) + 
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +
  scale_linetype_manual(values = c("controlNofish" = "dotted", 
                                     "Selective_attack" = "dashed", 
                                     "Non-selective" = "dashed", 
                                     "Selective_handling" = "solid",
                                     "Selective_both" = "solid"))
PredsA = ggplot(endpoints, aes(x = Preds, y = A, color = Scenario, linetype = Scenario )) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) + labs(x=expression('Fish Density, L' ^ -1), y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +
  scale_linetype_manual(values = c("controlNofish" = "dotted", 
                                     "Selective_attack" = "dashed", 
                                     "Non-selective" = "dashed", 
                                     "Selective_handling" = "solid",
                                     "Selective_both" = "solid"))

#####figure in main text####
library(ggpubr)
ggarrange(a1,b1,c1, PredsN,PredsJ,PredsA,
          labels = c("A","B","C","D","E","F"),
          nrow = 2, ncol = 3,
          common.legend = TRUE,
          legend = "none")

###figure for UVM

palette2 <- c(
  "controlNofish" = "grey",
  "Non-selective" = "red",
  "Selective_both" = "black"
)

NaupLongfiltered = NaupLong %>% filter(predation_scenario == "controlNofish" | predation_scenario == "Selective_both" | predation_scenario == "Non-selective")
JuvLongfiltered = JuvLong %>% filter(predation_scenario == "controlNofish" | predation_scenario == "Selective_both" | predation_scenario == "Non-selective")
AdLongfiltered = AdLong %>% filter(predation_scenario == "controlNofish" | predation_scenario == "Selective_both" | predation_scenario == "Non-selective")

a2 = ggplot(data=NaupLongfiltered, aes(x=Time, y=value, group=predation_scenario, color=predation_scenario)) + 
  geom_vline(xintercept = 49, color="blue", linewidth=1.5) + 
  geom_line(linewidth=1.5) +
  theme_classic(base_size = 20) + 
  labs(x="Time (days)", y=expression(Nauplii~Density~L^{-1})) + 
  # Use the named vector to map colors to data values
  scale_color_manual(values = palette2, 
                     name = "Scenario",
                     labels = c("Control (No Fish)", "No Preference", "Preference for Adults")) 


b2 = ggplot(data=JuvLongfiltered,aes(x=Time,y=value,group=predation_scenario,color=predation_scenario)) +
  geom_vline(xintercept = 49,color="blue",linewidth=1.5) + geom_line(linewidth=1.5) +theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Juvenile Density, L' ^ -1)) + 
  scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette2, name = "scenario") + 
  guides(color = guide_legend(title = "Fish Scenario")) +
  scale_y_continuous(limits = c(0, NA)) 
c2 = ggplot(data=AdLongfiltered,aes(x=Time,y=value,group=predation_scenario,color=predation_scenario)) +
  geom_vline(xintercept = 49,color="blue",linewidth=1.5) + geom_line(linewidth=1.5) +theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Adult Density, L' ^ -1)) + 
  scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette2, name = "scenario") + 
  guides(color = guide_legend(title = "Fish Scenario")) +
  scale_y_continuous(limits = c(0, NA)) 

endpointsfiltered = endpoints %>% filter(Scenario == "Selective_both" | Scenario == "Non-selective")

palette3 <- c(
  "Non-selective" = "red",
  "Selective_both" = "black"
)

PredsN2 = ggplot(endpointsfiltered, aes(x = Preds, y = N, group = Scenario, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Nauplii Density, L' ^ -1)) +
  scale_color_manual(values = palette3, name = "Scenario", labels = c("No Preference", "Preference for Adults")) 



PredsJ2 = ggplot(endpointsfiltered, aes(x = Preds, y = J,group = Scenario, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Juvenile Density, L' ^ -1)) +
  scale_color_manual(values = palette3, name = "Scenario", labels = c("No Preference", "Preference for Adults")) 



PredsA2 = ggplot(endpointsfiltered, aes(x = Preds, y = A,group = Scenario, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Adult Density, L' ^ -1)) +
  scale_color_manual(values = palette3, name = "Scenario", labels = c("No Preference", "Preference for Adults")) 




ggarrange(a2,b2,c2, PredsN2,PredsJ2,PredsA2,
          labels = c("A","B","C","D","E","F"),
          nrow = 2, ncol = 3,
          common.legend = TRUE,
          legend = "none")
