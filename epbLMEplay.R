#Load libraries
require(tidyverse)
require(ggpubr)
require(cowplot)
require(scales)
require(rstatix)
require(broom)
require(mgcv)
require(lme4)
require(lmerTest)
require(MuMIn)

#Load data
data <-read.csv(file="SteadyState.data.csv", header=TRUE, as.is=TRUE)

#Convert Ct value to copies/reaction in qPCR reaction
#Slope and intercept from qPCR standard curve for O. mykiss
#Equation: 10^(Ct value - intercept)/slope = copies/reaction 
slope <- -3.402646
intercept <- 38.63023
data <- data %>% 
  dplyr::mutate(copies.reaction = 10^((Ct_value - intercept)/slope))

#Convert copies per qPCR reaction to copies per liter of water filtered
#X (copies/L) =	(Y·b)/(c·a)
#X(copies/L)	The eDNA copy number in the original water sample
#a	volume of water filtered (in liters)
#b	the re-suspended volume of purified DNA (fixed b at 100 µL)
#c	the amount of purified DNA added to each qPCR-reaction volume (c: 6 µL/qPCR reaction)
#Y	the measured copy number in the qPCR-reaction volume, based on the standard curve (Y, copies/PCR reaction)
data <- data %>% 
  dplyr::mutate(copies.liter = (copies.reaction * Elution_volume)/Template * Volume_filtered) %>% 
  dplyr::mutate(log.copies.liter = log(copies.liter)) %>% 
  dplyr::mutate(individualWaterSample = paste0(Time,Water_sample))


avg.data <- data %>% 
  group_by(Fish_num) %>% 
  group_by(Replicate, .add=TRUE) %>% 
  group_by(Time, .add=TRUE) %>% 
  group_by(Water_sample, .add=TRUE) %>% 
  dplyr::filter(!is.na(Ct_value), .add=TRUE) %>%
  summarise(avg.copies.liter = mean(copies.liter), 
            sd.copies.liter = sd(copies.liter), 
            n.qPCR.reps = n(), 
            se.copies.liter = sd(copies.liter)/sqrt(n()),
            Group=mean(Group)) %>% 
  dplyr::mutate(log.copies.liter = log(avg.copies.liter)) %>% 
  dplyr::mutate(individualWaterSample = Water_sample)

data_05_1 <- data %>% 
  dplyr::filter(Fish_num == 5 & Replicate == 1) 
data_05_2 <- data %>% 
  dplyr::filter(Fish_num == 5 & Replicate == 2) 
data_25_1 <- data %>% 
  dplyr::filter(Fish_num == 25 & Replicate == 1) 
data_25_2 <- data %>% 
  dplyr::filter(Fish_num == 25 & Replicate == 2) 

avg.data_05_1 <- avg.data %>%
  dplyr::filter(Fish_num == 5 & Replicate == 1)
avg.data_05_2 <- avg.data %>%
  dplyr::filter(Fish_num == 5 & Replicate == 2)
avg.data_25_1 <- avg.data %>%
  dplyr::filter(Fish_num == 25 & Replicate == 1)
avg.data_25_2 <- avg.data %>%
  dplyr::filter(Fish_num == 25 & Replicate == 2)


#####
# plots for LME rand int; qPCR replicates not averaged

lmeOfDecay <- function(data) {
  
  nFish <- mean(data$Fish_num) %>% 
    as.numeric()
  
  replicate <- mean(data$Replicate) %>% 
    as.numeric()
  
  mod.a <- lmer(log.copies.liter ~ Time + (1|individualWaterSample),
                data = data)
  plotInuse <- ggplot(data) +
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample))) +
    theme(legend.position='none') +
    ggtitle(paste0(nFish,':',replicate,'  k = ',
                   round(summary(mod.a)$coefficients[2],5))) +
    coord_cartesian(ylim=c(0,13))
  
  fInt <- fixef(mod.a)[1]
  fSlp <- fixef(mod.a)[2]
  rInt <- ranef(mod.a) %>% unlist() %>% as.numeric()
  for (i in 1:length(rInt)) {
    plotInuse <- plotInuse + 
      geom_abline(intercept = fInt + rInt[i], slope = fSlp,alpha=0.4)
  }
  plotInuse <- plotInuse + 
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample))) 
  
  # hist(rInt)
  # print(summary(mod.a))
    
  return(plotInuse)
}

lmePlots <- ggarrange(lmeOfDecay(data_05_1),
                      lmeOfDecay(data_05_2),
                      lmeOfDecay(data_25_1),
                      lmeOfDecay(data_25_2) ,
                      nrow=2, ncol=2)



lmePlots


#####
# plots for LME rand int & slope; qPCR replicates not averaged

lme2OfDecay <- function(data) {


  nFish <- mean(data$Fish_num) %>%
    as.numeric()

  replicate <- mean(data$Replicate) %>%
    as.numeric()

  mod.a <- lmer(log.copies.liter ~  (1+Time|individualWaterSample),
                data = data, REML = FALSE,
                control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))
  plotInuse <- ggplot(data) +
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample))) +
    theme(legend.position='none') +
    ggtitle(paste0(nFish,':',replicate,'  k = ',
                   round(summary(mod.a)$coefficients[2],5))) +
    coord_cartesian(ylim=c(0,13))

  fInt <- fixef(mod.a)[1]
  fSlp <- fixef(mod.a)[2]
  rInt <- ranef(mod.a)$individualWaterSample$'(Intercept)' %>% unlist() %>% as.numeric()
  rSlp <- ranef(mod.a)$individualWaterSample$'Time' %>% unlist() %>% as.numeric()
  # rSlp <- ranef(mod.a) %>% unlist() %>% as.numeric()
  # print(ranef(mod.a))
  # print(rSlp)
  for (i in 1:length(rSlp)) {
    plotInuse <- plotInuse +
      geom_abline(intercept = fInt+rInt[i], slope = fSlp+rSlp[i],alpha=0.4)
    # geom_abline(intercept = fInt, slope = fSlp+rSlp[i],alpha=0.4)
  }
  plotInuse <- plotInuse +
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample)))

  # print(plotInuse)

  # plot(rInt)
  # print(summary(mod.a))

  return(plotInuse)
}

lme2Plots <- ggarrange(lme2OfDecay(data_05_1),
                      lme2OfDecay(data_05_2),
                      lme2OfDecay(data_25_1),
                      lme2OfDecay(data_25_2) ,
                      nrow=2, ncol=2)



lme2Plots


#####
# plots for lm across different time horizons; qPCR replicates averaged


restrictedTimeRangeDecay <- function(data) {
  
  nFish <- mean(data$Fish_num) %>% 
    as.numeric()
  
  replicate <- mean(data$Replicate) %>% 
    as.numeric()
  
  plotInuse <- ggplot(data) +
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample))) +
    theme(legend.position='none') +
    ggtitle(paste0(nFish,':',replicate))
  
  times <- unique(data$Time)
  
  for (t in seq(3,length(times))) {
  
    dataInuse <- data %>% 
      dplyr::filter(Time <= times[t])
    print(dim(dataInuse))
    
    mod <- lm(log.copies.liter ~ Time, dataInuse)
    print(mod$coefficients)
    int <- mod$coefficients[1] %>% as.numeric()
    slp <- mod$coefficients[2] %>% as.numeric()
    print(c(int,slp))
    plotData <- data.frame(Time=unique(dataInuse$Time))
    plotData <- plotData %>% 
      dplyr::mutate(log.copies.liter = int + slp * Time)

    plotInuse <- plotInuse + 
      geom_abline(intercept = int, slope = slp,alpha=0.2) +
      geom_line(data = plotData, aes(x=Time,y=log.copies.liter),linewidth=1) +
      coord_cartesian(ylim=c(0,13))
    
  }
  plotInuse <- plotInuse + 
    geom_point(aes(x=Time,y=log.copies.liter,color=as.factor(Water_sample))) 
  
  return(plotInuse)
}

restRangePlots <- ggarrange(restrictedTimeRangeDecay(avg.data_05_1),
                            restrictedTimeRangeDecay(avg.data_05_2),
                            restrictedTimeRangeDecay(avg.data_25_1),
                            restrictedTimeRangeDecay(avg.data_25_2) ,
                       nrow=2, ncol=2)



restRangePlots
