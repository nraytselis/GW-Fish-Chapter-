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
    
    Preds = ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
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
    
    Preds = ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
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
    
    Preds = ifelse(t<28,0,Preds) #fish added to tanks on week 4 
    
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
timespan = 365 



######Simulations Over Time#####

#Run simulations under different handling time and attack rate scenarios. 

#A) h and a (f) both are same for all stages 
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

#F) exclusive predation adults
other_parameters = c(Preds = 0.3,VOL = 1, f = 1, f_J = 0, f_N = 0,
h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
parameters = c(pars,other_parameters)
sim6 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)

#G) no preference
other_parameters = c(Preds = 0.3,VOL = 1, f = 1, f_J = 1, f_N = 1,
                     h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
parameters = c(pars,other_parameters)
sim7 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)

#H) preference for adults
other_parameters = c(Preds = 0.3,VOL = 1, f=1,f_J=0.1,f_N=0.01,
                     h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
parameters = c(pars,other_parameters)
sim8 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODE)


#I) exclusive pred N
other_parameters = c(Preds = 0.3,VOL = 1, f=1,f_J=0.1,f_N=0.01,
                     h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
parameters = c(pars,other_parameters)
sim9 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODEOnlyNPred) 

#J) exclusive pred J
other_parameters = c(Preds = 0.3,VOL = 1, f=1,f_J=0.1,f_N=0.01,
                     h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
parameters = c(pars,other_parameters)
sim10 = ode(y = Initial_conditions, times=1:timespan, parms=parameters, 
           method="lsoda", func=FishODEOnlyJPred) 

#combine all simulations into a single dataframe 
colnames(sim1) = c("time","N1","J1","A1")
colnames(sim2) = c("time","N2","J2","A2")
colnames(sim3) = c("time","N3","J3","A3")
colnames(sim4) = c("time","N4","J4","A4")
colnames(sim5) = c("time","N5","J5","A5")
colnames(sim6) = c("time","N6","J6","A6")
colnames(sim7) = c("time","N7","J7","A7")
colnames(sim8) = c("time","N8","J8","A8")
colnames(sim9) = c("time","N9","J9","A9")
colnames(sim10) = c("time","N10","J10","A10")
allscenarios = bind_cols(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8,sim9,sim10)
names(allscenarios)[1] <- "Time"


#Plot scenarios for 1) preference for adults, 2) no preference, and 3) control (sims 5,7,8)
Naup = allscenarios %>% select(c(1,18,26,30))
Juv = allscenarios %>% select(c(1,19,27,31))
Ad = allscenarios %>% select(c(1,20,28,32))
NaupLong = Naup %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )
JuvLong = Juv %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )
AdLong = Ad %>% pivot_longer(cols=c(2:4),names_to = "scenario",values_to = "value" )

NaupLong$scenario <- gsub("N5", "control", NaupLong$scenario)
NaupLong$scenario <- gsub("N7", "nopref", NaupLong$scenario)
NaupLong$scenario <- gsub("N8", "prefadults", NaupLong$scenario)

JuvLong$scenario <- gsub("J5", "control", JuvLong$scenario)
JuvLong$scenario <- gsub("J7", "nopref", JuvLong$scenario)
JuvLong$scenario <- gsub("J8", "prefadults", JuvLong$scenario)

AdLong$scenario <- gsub("A5", "control", AdLong$scenario)
AdLong$scenario <- gsub("A7", "nopref", AdLong$scenario)
AdLong$scenario <- gsub("A8", "prefadults", AdLong$scenario)


#plotting 

palette <- c(
  control   = "darkgoldenrod2",
  prefadults = "deeppink2",
  nopref  ="#77AADD"
)

a1 = ggplot(data=NaupLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) +theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) + scale_y_log10()

b1 = ggplot(data=JuvLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) + scale_y_log10()

c1 = ggplot(data=AdLong,aes(x=Time,y=value,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + theme_classic(base_size = 20)+
  labs(x="Time (days)",y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",nopref="No Preference",prefadults="Preference for Adults"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) + scale_y_log10()

library(ggpubr)
ggarrange(a1,b1,c1,
          nrow = 1, ncol = 3,
          common.legend = TRUE,
          legend = "right")


#plot 1) exclusive for adults, 2) no preference, and 3) control
Naup = allscenarios %>% select(c(1,18,22,26,30,34,38))
Juv = allscenarios %>% select(c(1,19,23,27,31,35,39))
Ad = allscenarios %>% select(c(1,20,24,28,32,36,40))
NaupLong = Naup %>% pivot_longer(cols=c(2:7),names_to = "scenario",values_to = "value" )
JuvLong = Juv %>% pivot_longer(cols=c(2:7),names_to = "scenario",values_to = "value" )
AdLong = Ad %>% pivot_longer(cols=c(2:7),names_to = "scenario",values_to = "value" )


#relabel scenarios to be more informative 
NaupLong$scenario <- gsub("N6", "exclusiveadults", NaupLong$scenario)
NaupLong$scenario <- gsub("N7", "nopref", NaupLong$scenario)
NaupLong$scenario <- gsub("N8", "prefadults", NaupLong$scenario)
NaupLong$scenario <- gsub("N9", "exclusiveNaup", NaupLong$scenario)
NaupLong$scenario <- gsub("N10", "exclusiveJuv", NaupLong$scenario)
NaupLong$scenario <- gsub("N5", "control", NaupLong$scenario)
JuvLong$scenario <- gsub("J6", "exclusiveadults", JuvLong$scenario)
JuvLong$scenario <- gsub("J7", "nopref", JuvLong$scenario)
JuvLong$scenario <- gsub("J8", "prefadults", JuvLong$scenario)
JuvLong$scenario <- gsub("J9", "exclusiveNaup", JuvLong$scenario)
JuvLong$scenario <- gsub("J10", "exclusiveJuv", JuvLong$scenario)
JuvLong$scenario <- gsub("J5", "control", JuvLong$scenario)
AdLong$scenario <- gsub("A6", "exclusiveadults", AdLong$scenario)
AdLong$scenario <- gsub("A7", "nopref", AdLong$scenario)
AdLong$scenario <- gsub("A8", "prefadults", AdLong$scenario)
AdLong$scenario <- gsub("A9", "exclusiveNaup", AdLong$scenario)
AdLong$scenario <- gsub("A10", "exclusiveJuv", AdLong$scenario)
AdLong$scenario <- gsub("A5", "control", AdLong$scenario)

# offset lines so you can see the control 
offset_factor <- c(
  control = 1.00,
  exclusiveadults = 1.02,
  nopref = 0.98,
  prefadults = 1.04,
  exclusiveNaup = 0.96,
  exclusiveJuv = 1.06
)

# Apply the offsets to each dataset
NaupLong <- NaupLong %>%
  mutate(value_offset = value * offset_factor[scenario])

JuvLong <- JuvLong %>%
  mutate(value_offset = value * offset_factor[scenario])

AdLong <- AdLong %>%
  mutate(value_offset = value * offset_factor[scenario])



palette = c(exclusiveadults = "#EE8866", nopref = "#77AADD",prefadults = "thistle",exclusiveNaup = "deeppink3", exclusiveJuv = "goldenrod", control = "#BBCC33")

a = ggplot(data=NaupLong,aes(x=Time,y=value_offset,group=scenario,color=scenario)) + geom_line(linewidth=1.5) +theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) + scale_y_log10()

b = ggplot(data=JuvLong,aes(x=Time,y=value_offset,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + theme_classic(base_size = 20) +
  labs(x="Time (days)",y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                 exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +  scale_y_log10()

c = ggplot(data=AdLong,aes(x=Time,y=value_offset,group=scenario,color=scenario)) + geom_line(linewidth=1.5) + theme_classic(base_size = 20)+
  labs(x="Time (days)",y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(control = "Control (No Fish)",exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                            exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) +  scale_y_log10()


######Fish Density Gradient#####

#Now run the model across a fish density gradient 

library(ggpubr)
ggarrange(a,b,c,
          nrow = 1, ncol = 3,
          common.legend = TRUE,
          legend = "right")

#predator gradient
preds_gradient <- seq(0, 0.4, by = 0.01)

#scenarios
scenario_params <- list(
  nopref = c(VOL = 1, f = 1, f_J = 1, f_N = 1,
             h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0),
  
  prefadults = c(VOL = 1, f=1,f_J=0.1,f_N=0.01,
                 h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
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

PredsN = ggplot(endpoints, aes(x = Preds, y = N, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                      exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario"))

PredsJ = ggplot(endpoints, aes(x = Preds, y = J, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                       exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) 

PredsA = ggplot(endpoints, aes(x = Preds, y = A, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) + labs(x=expression('Fish Density, L' ^ -1), y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                                                    exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) 


ggarrange(a1,b1,c1,PredsN,PredsJ,PredsA,
          labels = c("A","B","C","D","E","F"),
          nrow = 2, ncol = 3,
          common.legend = TRUE,
          legend = "right")


#####Stage Specific ODEs#####
####run sims with nauplii and juvenile only ODEs
scenario_params_N <- list(
  exclusiveNaup = c(VOL = 1, f=1,f_J=0.1,f_N=0.01, #doesn't matter what f is here because no predA and predJ in model
                 h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
)

sim_results_N <- list()

for (s in names(scenario_params_N)) {
  for (p in preds_gradient) {
    
    # update parameters
    parameters <- c(pars,
                    Preds = p,
                    scenario_params_N[[s]])
    
    # run simulation
    sim <- ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = parameters,
      method = "lsoda",
      func = FishODEOnlyNPred
    )
    
    # store results
    sim_df <- as.data.frame(sim) %>%
      mutate(Preds = p,
             Scenario = s)
    
    sim_results_N[[paste(s, p, sep = "_")]] <- sim_df
  }
}

all_sims_N <- bind_rows(sim_results_N, .id = "run_id")


###sim for J
scenario_params_J <- list(
  exclusiveJuv = c(VOL = 1, f=1,f_J=0.1,f_N=0.01, #doesn't matter what f is here because no predA and predJ in model
                  h = 0.006, h_J = 0.1, h_N = 0.01, i_P = 0)
)

sim_results_J <- list()

for (s in names(scenario_params_J)) {
  for (p in preds_gradient) {
    
    # update parameters
    parameters <- c(pars,
                    Preds = p,
                    scenario_params_J[[s]])
    
    # run simulation
    sim <- ode(
      y = Initial_conditions,
      times = 1:timespan,
      parms = parameters,
      method = "lsoda",
      func = FishODEOnlyJPred
    )
    
    # store results
    sim_df <- as.data.frame(sim) %>%
      mutate(Preds = p,
             Scenario = s)
    
    sim_results_J[[paste(s, p, sep = "_")]] <- sim_df
  }
}

all_sims_J <- bind_rows(sim_results_J, .id = "run_id")


all_sims <- bind_rows(all_sims,all_sims_N,all_sims_J)

library(dplyr)

#take the last time point from each simulation
endpoints <- all_sims %>%
  group_by(Scenario, Preds) %>%
  slice_tail(n = 1) %>%   # last time step
  ungroup()

PredsN = ggplot(endpoints, aes(x = Preds, y = N, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Nauplii Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                    exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario"))

PredsJ = ggplot(endpoints, aes(x = Preds, y = J, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) +
  labs(x=expression('Fish Density, L' ^ -1), y=expression('Juvenile Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                      exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) 

PredsA = ggplot(endpoints, aes(x = Preds, y = A, color = Scenario)) +
  geom_line(size = 1.5) +
  theme_classic(base_size = 20) + labs(x=expression('Fish Density, L' ^ -1), y=expression('Adult Density, L' ^ -1)) + scale_color_manual(labels = c(exclusiveadults ="Exclusive for Adults",nopref="No Preference",prefadults="Preference for Adults",exclusiveNaup="Exclusive for Nauplii",
                                                                                                                      exclusiveJuv="Exclusive for Juveniles"),values = palette, name = "scenario") + guides(color = guide_legend(title = "Fish Scenario")) 


ggarrange(a,b,c,PredsN,PredsJ,PredsA,
          nrow = 2, ncol = 3,
          common.legend = TRUE,
          legend = "right")

