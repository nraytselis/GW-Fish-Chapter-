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

#calculate proportions
Fish_Time_Series = Fish_Time_Series %>% mutate(PropAdults=Adults_per_Tank/(Adults_per_Tank+Juveniles_per_tank+Nauplii_per_Tank))
Fish_Time_Series = Fish_Time_Series %>% mutate(PropJuv=Juveniles_per_tank/(Adults_per_Tank+Juveniles_per_tank+Nauplii_per_Tank))
Fish_Time_Series = Fish_Time_Series %>% mutate(PropNaup=Nauplii_per_Tank/(Adults_per_Tank+Juveniles_per_tank+Nauplii_per_Tank))

Fish_Time_Series_summary = Fish_Time_Series %>% group_by(Preds,Day,FishDensity) %>% summarise(mean_A = mean(PropAdults), 
                                                                                              mean_J = mean(PropJuv), 
                                                                                              mean_N = mean(PropNaup))
Fish_Time_Series_summary_long <- Fish_Time_Series_summary %>%
  pivot_longer(cols = starts_with("mean"), names_to = "Stage", values_to = "Proportion")


palette = c(mean_A="deeppink3",mean_J="#77AADD",mean_N="goldenrod")
ggplot(Fish_Time_Series_summary_long, aes(x = Day, y = Proportion, fill = Stage)) +
  geom_area(position = "stack") +
  labs(
    x = "Time",
    y = "Proportion",
    fill = "Variable"
  ) + facet_wrap(~ FishDensity) + 
  theme_classic(base_size=20) +
  scale_fill_manual(
    labels = c(mean_A = "Adults", mean_J = "Juveniles", mean_N = "Nauplii"),
    values = palette,
    name = "Stage"
  ) 





