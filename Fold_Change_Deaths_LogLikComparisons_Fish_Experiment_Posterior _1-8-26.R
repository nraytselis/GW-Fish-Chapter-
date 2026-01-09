
setwd("~/Desktop/Rscripts/Data/Fish")
Fish_Time_Series = read.csv("Fish_experiment_data.csv")

#make sure 11 fish treatment is cleaned up 
Fish_Time_Series$Preds[Fish_Time_Series$Preds == 11] <- 12

#sort by Day 
Fish_Time_Series <- Fish_Time_Series[order(Fish_Time_Series$Tank, Fish_Time_Series$Day), ]
Fish_Time_Series$Day   <- as.numeric(Fish_Time_Series$Day)

#make fish-day and preds numeric
Fish_Time_Series$Live_Fish <- as.numeric(Fish_Time_Series$Live_Fish)
Fish_Time_Series$Preds   <- as.numeric(Fish_Time_Series$Preds)

#rename tank column unit
colnames(Fish_Time_Series)[3] = "unit"

#bring in fish treatment info
Fishtreatments = read_csv("FishTreatments.csv") 

#merge fish raw data with treatment info 
Fish_Time_Series <- left_join(Fish_Time_Series, Fishtreatments, by = "unit")
Fish_Time_Series$FishDensity[Fish_Time_Series$FishDensity == 0.275] <- 0.3

#calculate copepods per tank (one replicated)
Fish_Time_Series = Fish_Time_Series %>% mutate(Adults_per_Tank = 0.5*(AF1/counted_volume1+AF2/counted_volume2)*20*40,
                                               Juveniles_per_tank = 0.5*(JOA1/counted_volume1+JOA2/counted_volume2)*20*VOL,
                                               Nauplii_per_Tank = 0.5*(N1/counted_volume1+N2/counted_volume2)*20*VOL,
                                               total = 0.5*((AF1+JOA1+N1)/counted_volume1+(AF2+JOA2+N2)/counted_volume2)*20*VOL) 

#filter missing rows
Fish_Time_Series = Fish_Time_Series[is.finite(Fish_Time_Series$Adults_per_Tank),]

Fish_Time_Series$Densitytotal = Fish_Time_Series$total/Fish_Time_Series$VOL
Fish_Time_Series$DensityN = Fish_Time_Series$Nauplii_per_Tank/Fish_Time_Series$VOL
Fish_Time_Series$DensityJ = Fish_Time_Series$Juveniles_per_tank/Fish_Time_Series$VOL
Fish_Time_Series$DensityA = Fish_Time_Series$Adults_per_Tank/Fish_Time_Series$VOL

#ratios of stages
mean(Fish_Time_Series$Nauplii_per_Tank/Fish_Time_Series$total) #16/29
mean(Fish_Time_Series$Juveniles_per_tank/Fish_Time_Series$total)#12/29
mean(Fish_Time_Series$Adults_per_Tank/Fish_Time_Series$total) #1/29


Fish_Time_Series_control = Fish_Time_Series %>% filter (Preds == 0)


Fish_Time_Series_control$Densitytotal = Fish_Time_Series_control$total/Fish_Time_Series_control$VOL
Fish_Time_Series_control$DensityN = Fish_Time_Series_control$Nauplii_per_Tank/Fish_Time_Series_control$VOL
Fish_Time_Series_control$DensityJ = Fish_Time_Series_control$Juveniles_per_tank/Fish_Time_Series_control$VOL
Fish_Time_Series_control$DensityA = Fish_Time_Series_control$Adults_per_Tank/Fish_Time_Series_control$VOL

#ratios of stages in controls (use these ratios)
mean(Fish_Time_Series_control$Nauplii_per_Tank/Fish_Time_Series_control$total) #8/29
mean(Fish_Time_Series_control$Juveniles_per_tank/Fish_Time_Series_control$total)#7/16
mean(Fish_Time_Series_control$Adults_per_Tank/Fish_Time_Series_control$total) #1/16







chainsA <- readRDS("Joint_GW_full_A2_2.RDA")
chainsB <- readRDS("Joint_GW_full_B2_2.RDA")
chainsC <- readRDS("Joint_GW_full_C2_2.RDA")

#remove chains that did not converge
chainsAconverged = chainsA[-c(1:3)]
chainsBconverged = chainsB[-c(2)]
chainsCconverged = chainsC[-c(1,6,7)]

#combine converged chains
fishchains = c(chainsAconverged,chainsBconverged,chainsCconverged)


#extract samples from all chains in the list 
get_samples <- function(chain.list) {
  L <- length(chain.list)
  samples_list <- vector("list", L)  # empty list to store samples
  
  for (i in 1:L) {
    samples <- chain.list[[i]]$samples  # extract 'samples' from each chain
    samples_list[[i]] <- coda::mcmc(samples)  # convert to mcmc object
  }
  
  return(samples_list)
}

all_chains <- get_samples(fishchains)



#thin all chains
all_chains_coda <- convert.to.coda(all_chains)
all_chains_coda_thin <- window(all_chains_coda,start=60000,end=100000,thin=1000)

#make a combined list with all parameter sets from each chain
pars <- as.mcmc.list(all_chains_coda_thin)

#loop through all parameter sets
pars_matrix <- as.matrix(pars)  

pars_matrix = data.frame(pars_matrix)

#back transform the log-transformed parameters
pars_matrix$comp_d <- exp(pars_matrix$comp_d)
pars_matrix$comp_b <- exp(pars_matrix$comp_b)
pars_matrix$comp_m <- exp(pars_matrix$comp_m)
pars_matrix$cann <- exp(pars_matrix$cann)

pars_matrix = as.matrix(pars_matrix)


n_iter <- nrow(pars_matrix)

#set up ranges
Total_Abundance_range <- 10 * (0:10000) #calculate from raw data 
n_Total_Abundance<- length(Total_Abundance_range)

#adult_range <- 10 * (0:500)  # 501 values
#n_adults <- length(adult_range)

#results matrices
matrix_Total_Abundance_deaths <- matrix(NA, nrow = n_iter, ncol = n_Total_Abundance)


Per_capita_deaths = function(d_N,comp_d,c_N,c_J,f,h,f_J,h_J,f_N,h_N,i_P,Total_Abundance,cann,VOL=40, Preds = 0){
  N = 8/16*Total_Abundance
  J = 7/16*Total_Abundance
  A = 1/16*Total_Abundance
  Pred_N = f*(Preds/VOL)/(1 + f*h*(A+f_J*h_J*J+f_N*h_N*N)/VOL  + i_P*max(Preds-1, 0)/VOL) #per capita
  d_N*exp(comp_d / VOL * (c_N * N + c_J * J + A)) + (cann*A)/VOL + Pred_N
}

#deaths
for (i in 1:n_iter) {
  p <- pars_matrix[i, ]
  
  matrix_Total_Abundance_deaths[i, ] <- Per_capita_deaths(
    d_N    = p["d_N"],
    comp_d = p["comp_d"],
    c_N    = p["c_N"],
    c_J    = p["c_J"],
    f = p["f"],
    h = p["h"],
    f_J = p["f_J"],
    h_J = p["h_J"],
    f_N = p["f_N"],
    h_N = p["h_N"],
    i_P = p["i_P"],
    cann   = p["cann"],
    Total_Abundance = Total_Abundance_range
  )
}

Hi = apply(matrix_Total_Abundance_deaths, 2, quantile, 0.975)
Per_capita_death_rate = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.5)
Lo = apply(matrix_Total_Abundance_deaths, 2, quantile, probs=0.025)
plotting_data_TotalAbundance = data.frame("Total_Abundance" = Total_Abundance_range, Hi, Per_capita_death_rate, Lo)
deaths_Total = ggplot(data=plotting_data_TotalAbundance, aes(x=Total_Abundance, y=Per_capita_death_rate)) + 
  geom_line() + geom_ribbon(aes(ymin=Lo, ymax=Hi, alpha=0.2),fill = "gold2") + labs(x="Total Abundance", y="Per Capita Death Rate") + theme_classic(base_size = 13) + theme(legend.position="none") 


deaths_Total

max = max(plotting_data_TotalAbundance$Per_capita_death_rate)
min = min(plotting_data_TotalAbundance$Per_capita_death_rate)

foldchangeD = (max-min)/min

foldchangeD

#parse apart likelihood from converged chains 
ChainLL = matrix(unlist(fishchains[[1]]$extra.values),ncol=3,nrow=250000,byrow = TRUE) 
colnames(ChainLL) = c("Prior","Rebound","Fish")
plot(ChainLL[,1], type = "l")
plot(ChainLL[,2], type = "l")
plot(ChainLL[,3], type = "l")


ChainLLdf = data.frame(matrix(unlist(fishchains[[1]]$extra.values),ncol=3,nrow=250000,byrow = TRUE)) 
ChainLLdf[,"iteration"] = 1:250000

for(i in 2:length(fishchains)){
  temp = data.frame(matrix(unlist(fishchains[[i]]$extra.values),ncol=3,nrow=250000,byrow = TRUE))
  temp[,"iteration"] = 1:250000
  ChainLLdf = rbind(ChainLLdf,temp)  
}

colnames(ChainLLdf)[1:3] = c("Prior","Rebound","Fish")

ChainLLdfnoBI = ChainLLdf %>% filter(iteration > 50000)

ggplot(ChainLLdfnoBI, aes(x=Rebound,y=Fish)) + geom_bin2d()

mean(ChainLLdfnoBI$Rebound) #-12441.75
max(ChainLLdfnoBI$Rebound) #-12425.71



