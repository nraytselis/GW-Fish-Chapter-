library(panelPomp)
library(adaptMCMC)
library(deSolve) 
library(ggplot2) 
library(tidyr)
library(ggsci)
library(deSolve)
library(mgcv)
library(itsadug)
library(cowplot)
library(wesanderson)
library(stringr)
library(dplyr)
library(tidyverse)
library(coda)

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

#bring in full_df from cluster
full_df = readRDS("Fish_predictions.RDA")  

#bring in fish treatment info
Fishtreatments = read_csv("FishTreatments.csv") 

#merge dataframes with treatment info and sim
merged_df <- left_join(full_df, Fishtreatments, by = "unit")

#replace 0.275 with 0.3 (12 fish)
merged_df$FishDensity[merged_df$FishDensity == 0.275] <- 0.3
merged_df$FishDensity[merged_df$FishDensity == 11] <- 12

#convert to densities (incorrect because Vol is not always 40)
#merged_df = merged_df %>% mutate(AF1Density = AF1/40, JOA1Density = JOA1/40, N1Density = N1/40, totalDensity = total/40) 


#summary stats (per tank instead of density)
merged_df_summary = merged_df %>% group_by(FishDensity,time) %>% summarise(mean_A = mean(AF1), 
                                                                           mean_J = mean(JOA1), 
                                                                           mean_N = mean(N1),
                                                                           mean_total = mean(total),
                                                                           sd_A = sd(AF1),
                                                                           sd_J = sd(JOA1),
                                                                           sd_N = sd(N1),
                                                                           sd_total = sd(total),
                                                                           n_A = n(),
                                                                           n_J = n(),
                                                                           n_N = n(),
                                                                           n_total = n(),
                                                                           se_A = sd_A / sqrt(n_A),
                                                                           se_J = sd_J / sqrt(n_J),
                                                                           se_N = sd_N / sqrt(n_N),
                                                                           se_total = sd_total / sqrt(n_total),
                                                                           lower_ci_A = quantile(AF1, prob=0.025),
                                                                           upper_ci_A = quantile(AF1, prob=0.975),
                                                                           lower_ci_J = quantile(JOA1, prob=0.025),
                                                                           upper_ci_J = quantile(JOA1, prob=0.975),
                                                                           lower_ci_N =quantile(N1, prob=0.025),
                                                                           upper_ci_N = quantile(N1, prob=0.975),
                                                                           lower_ci_total =quantile(total, prob=0.025),
                                                                           upper_ci_total =quantile(total, prob=0.975)) 

#merge fish raw data with treatment info 
Fish_Time_Series <- left_join(Fish_Time_Series, Fishtreatments, by = "unit")
Fish_Time_Series$FishDensity[Fish_Time_Series$FishDensity == 0.275] <- 0.3

#calculate copepods per tank (one replicated)
Fish_Time_Series = Fish_Time_Series %>% mutate(Adults_per_Tank = 0.5*(AF1/counted_volume1+AF2/counted_volume2)*20*40,
                                               Juveniles_per_tank = 0.5*(JOA1/counted_volume1+JOA2/counted_volume2)*20*VOL,
                                               Nauplii_per_Tank = 0.5*(N1/counted_volume1+N2/counted_volume2)*20*VOL,
                                               total = 0.5*((AF1+JOA1+N1)/counted_volume1+(AF2+JOA2+N2)/counted_volume2)*20*VOL) 


#write.csv(Fish_Time_Series, file = "Fish_Time_Series_Densities_11-5-25.csv", row.names = FALSE)

#filter missing rows
Fish_Time_Series = Fish_Time_Series[is.finite(Fish_Time_Series$Adults_per_Tank),]

#covert to densities (this is incorrect because Vol is not always 40)
#Fish_Time_Series = Fish_Time_Series %>% mutate(AdultDensity = Adults_per_Tank/40,JuvenileDensity = Juveniles_per_tank/40,
# NaupliiDensity = Nauplii_per_Tank/40, totalDensity = total/40)

#calculate summary statistics for raw data (abundances)
Fish_Time_Series_summary = Fish_Time_Series %>% group_by(Day,FishDensity) %>% summarise(mean_A = mean(Adults_per_Tank), 
                                                                                              mean_J = mean(Juveniles_per_tank), 
                                                                                              mean_N = mean(Nauplii_per_Tank),
                                                                                              mean_total = mean(total),
                                                                                              n_A = n(),
                                                                                              n_J = n(),
                                                                                              n_N = n(),
                                                                                              n_total = n(),
                                                                                              se_N = sd(Nauplii_per_Tank) / sqrt(n_N),
                                                                                              se_J = sd(Juveniles_per_tank)/ sqrt(n_J),
                                                                                              se_A = sd(Adults_per_Tank)/sqrt(n_A),
                                                                                              se_total = sd(total) / sqrt(n_total)) 


total = ggplot(merged_df_summary, aes(x = time, y = mean_total)) +
  geom_ribbon(aes(ymin=lower_ci_total, ymax=upper_ci_total), 
              fill = "grey", alpha = 0.5) +
  geom_line(color = "black", linewidth=1) +
  theme_classic() + 
  geom_line(data = Fish_Time_Series_summary, 
            aes(x = Day, y = mean_total), 
            color = "#3C7DB1", linewidth=1) +   
  geom_errorbar(data = Fish_Time_Series_summary, color = "#3C7DB1",
                aes(x = Day, ymin = mean_total - se_total, ymax = mean_total + se_total)) +
  geom_point(data = Fish_Time_Series_summary, 
             aes(x = Day, y = mean_total), 
             color = "#3C7DB1", size = 2) + 
  facet_wrap(~ FishDensity,
             nrow = 3, ncol = 3) +   
  labs(x = "Time (days)", y = expression('Total Copepod Density, L' ^ -1)) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none") 

total




