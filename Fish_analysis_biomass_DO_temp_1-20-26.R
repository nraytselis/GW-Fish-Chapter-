setwd("~/Desktop/Rscripts/Data/Fish")

FishData = read.csv("Fish_Biomass_DO_Temp_1-20-26.csv") 

FishData <- FishData %>%
  mutate(time = rep(1:30, length.out = n()))

library(mgcv)

library(glmmTMB)
FishData = na.omit(FishData)

FishData$Treatment <- as.character(FishData$Treatment)
FishData$Treatment[FishData$Treatment == "11"] <- "12"
FishData$Treatment <- factor(FishData$Treatment)

#make sure treatment is continuous
FishData$Number_of_Fish <- as.numeric(as.character(FishData$Number_of_Fish))
FishData$Treatment      <- as.numeric(as.character(FishData$Treatment))


#exclude controls
FishData = FishData %>% filter(Treatment!= "0")

#fish survival
FishData$PercentSurvival = 100*(FishData$Number_of_Fish/FishData$Treatment)
FishDataSurvival = FishData %>% filter(Date == "5/3/24")
FishDataSurvival = na.omit(FishDataSurvival)
sum(FishDataSurvival$PercentSurvival)/30

#exclude other treatments with zero fish (due to loss of fish during experiment)
FishData = FishData %>% filter(Number_of_Fish != "0")


#center covariates
FishData$DO_c   <- scale(FishData$DO, center = TRUE, scale = FALSE)
FishData$time_c <- scale(FishData$time, center = TRUE, scale = FALSE)
FishData$Treatment_c   <- scale(FishData$Treatment, center = TRUE, scale = FALSE)

#The main effect of DO_c = effect of DO at average time

#The main effect of time_c = effect of time at average DO


FishModel <- glmmTMB(
  Mean_Fish_Biomass ~ DO * Treatment * time + (1 | Tank),
  family = Gamma(link = "log"),
  data = FishData
)


summary(FishModel)

FishModelc <- glmmTMB(
  Mean_Fish_Biomass ~ DO_c * Treatment_c * time_c + (1 | Tank),
  family = Gamma(link = "log"),
  data = FishData
)

summary(FishModelc)


FishDataSummary = FishData %>% group_by(Date,Treatment) %>% summarise(mean=mean(Mean_Fish_Biomass),
                                                                      sd = sd(Mean_Fish_Biomass),
                                                                      n = n(),
                                                                      se_Mean_Fish_Biomass = sd(Mean_Fish_Biomass) / sqrt(n))

library(lubridate)
library(RColorBrewer)
FishDataSummary$Date <- mdy(FishDataSummary$Date)

my_palette <- brewer.pal(n = 5, name = "PuBu")[2:6]

ggplot(data=FishDataSummary,aes(x=Date,y=mean,group=as.factor(Treatment),color=as.factor(Treatment))) + geom_line(linewidth=1) +
  geom_errorbar(data = FishDataSummary,
                aes(x = Date, ymin = mean - se_Mean_Fish_Biomass, ymax = mean + se_Mean_Fish_Biomass)) + labs(x = "Day", y = "Mean Fish Biomass (g)") + theme_classic() + scale_color_manual(labels = c("0.025","0.05","0.1","0.175","0.3"),values = my_palette) + theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18)
  ) + guides(color = guide_legend(title = expression('Fish Density, L' ^ -1),  ))




