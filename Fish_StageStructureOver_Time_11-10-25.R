#Packages
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

setwd("~/Desktop/Rscripts/Data/Fish")

####Summary####
#This code plots the proportion of individuals of each stage class throughout the fish mesocosm experiment

Fish_Time_Series_og = read.csv("Fish_experiment_data.csv")

#format this dataframe to look like fish_time_series
Fish_Time_Series2 = read.csv("Fish_Data_Week2on_10-22-25.csv") 

# 
# Fish_Time_Series2 <- Fish_Time_Series2 %>%
#   mutate(NumericID = as.numeric(gsub("[A-Z]", "", Sample))) %>% # Extract numeric part of ID
#   group_by(Day, NumericID,Tank) %>%
#   summarise(AF1 = mean(AF), JOA1 = mean(JOA), N1 = mean(N), counted_vol = mean(counted_vol), Preds = mean(Preds)) %>%
#   ungroup()


Fish_Time_Series_og = Fish_Time_Series_og[, c(2,3,7)]
#convert 28 day to 35 because they have the same number of livefish
Fish_Time_Series_og$Day[Fish_Time_Series_og$Day == 28] <- 35

Fish_Time_Series2$number <- as.numeric(gsub("[A-Z]", "", Fish_Time_Series2$Sample))
Fish_Time_Series2$letter <- gsub("\\d", "", Fish_Time_Series2$Sample)

Fish_Time_Series2A = Fish_Time_Series2 %>% filter(letter == "A")

Fish_Time_Series2B = Fish_Time_Series2 %>% filter(letter == "B")

colnames(Fish_Time_Series2A) = c("Week","Day1","Tank", "Sample","AF1","JOA1","N1","counted_volume1","Preds","number","letter")

colnames(Fish_Time_Series2B) = c("Week","Day1","Tank", "Sample","AF2","JOA2","N2","counted_volume2","Preds","number","letter")

Fish_Time_Series2AB = left_join(Fish_Time_Series2A,Fish_Time_Series2B,by=c("Week","Day1","Tank","Preds"))

Fish_Time_Series2AB = Fish_Time_Series2AB[, -c(4,10,11,12,17,18)]

names(Fish_Time_Series2AB)[names(Fish_Time_Series2AB) == 'Day1'] <- 'Day'

#Fish_Time_Series2ABlivefish <- merge(
#  Fish_Time_Series2AB, 
#  Fish_Time_Series_og,
#  by = c("Tank", "Day"),   # merge only where both Tank & Day match
#  all.x = TRUE             # keep all rows from the first df
#)

Fish_Time_Series2ABlivefish <- Fish_Time_Series_og %>%
  left_join(Fish_Time_Series2AB, by = c("Tank", "Day"))


#make fish from day 35 equal to zero instead of NA
#Fish_Time_Series2ABlivefish$Live_Fish[Fish_Time_Series2ABlivefish$Day == "35"] <- 0

#make preds zero for earlier timepoints
#Fish_Time_Series2ABlivefish$Preds[Fish_Time_Series2ABlivefish$Day <= "42"] <- 0


Fish_Time_Series = Fish_Time_Series2ABlivefish

#make sure 11 fish treatment is cleaned up 
Fish_Time_Series$Preds[Fish_Time_Series$Preds == 11] <- 12

Fish_Time_Series = Fish_Time_Series %>% mutate(FishDensity = (Preds/40),Adults_per_Tank = 0.5*(AF1/counted_volume1+AF2/counted_volume2)*20*40,
                                               Juveniles_per_tank = 0.5*(JOA1/counted_volume1+JOA2/counted_volume2)*20*40,
                                               Nauplii_per_Tank = 0.5*(N1/counted_volume1+N2/counted_volume2)*20*40,
                                               total = 0.5*((AF1+JOA1+N1)/counted_volume1+(AF2+JOA2+N2)/counted_volume2)*20*40) 



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


