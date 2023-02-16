################################################################################################
################################################################################################
# 2022 Steady State Experiments ##
# Cal Poly Humboldt
# Batch reactor experiments to determine environmental DNA shedding and decay rates in steelhead
# Andrew Kinziger, Gavin Bandy, Edwin Millard, Braden Herman, Jason Shaffer, Andre Buchheister, and Eric Bjorkstedt
################################################################################################
################################################################################################

#Set working directory
setwd("/Users/apk7001/Desktop/SteadyState")

#Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
require(tidyverse)
require(scales)
library(rstatix)
library(broom)

#Load data
data <-read.csv(file="SteadyState.data.csv", header=TRUE)
data$Ct_value <- as.numeric(data$Ct_value)

#Convert Ct value to copies/reaction in qPCR reaction
#Slope and intercept from qPCR standard curve for O. mykiss
#Equation: 10^(Ct value - intercept)/slope = copies/reaction 
slope <- -3.402646
intercept <- 38.63023
data$copies.reaction <- 10^((data$Ct_value - intercept)/slope)

#Convert copies per qPCR reaction to copies per liter of water filtered
#X (copies/L) =	(Y·b)/(c·a)
#X(copies/L)	The eDNA copy number in the original water sample
#a	volume of water filtered (in liters)
#b	the re-suspended volume of purified DNA (fixed b at 100 µL)
#c	the amount of purified DNA added to each qPCR-reaction volume (c: 6 µL/qPCR reaction)
#Y	the measured copy number in the qPCR-reaction volume, based on the standard curve (Y, copies/PCR reaction)
data$copies.liter <- (data$copies.reaction * data$Elution_volume)/(data$Template * data$Volume_filtered)

#Check of copies per qPCR reaction by using another method
#CN/L=CN.PCR x (Ve/Vt) * (1/Vs), Ve is elution volume (uL), Vt is template volume (uL), and Vs is volume of water filtered in L
cor(data$copies.liter,data$copies.reaction * (data$Elution_volume/data$Template) * (1/data$Volume_filtered))
#Here we get and expect a perfect correlation =1

#Add column for ln(copies.liter).  The default for "log" in R is natural log
data$ln.copies.liter <- log(data$copies.liter)


########################################################################
########################################################################
########################################################################
#Decay and shedding rates
#Estimated following:
#Sansom, B. J., & Sassoubre, L. M. (2017). Environmental DNA (eDNA) shedding and decay rates to model freshwater mussel eDNA transport in a river. Environmental Science & Technology, 51(24), 14244– 14253. https://doi.org/10.1021/acs.est.7b05199
#Sassoubre, L. M., Yamahara, K. M., Gardner, L. D., Block, B. A., & Boehm, A. B. (2016). Quantification of environmental DNA (eDNA) shedding and decay rates for three marine fish. Environmental Science & Technology, 50(19), 10456– 10464.
########################################################################
########################################################################
########################################################################


####################################
#Decay Rate
####################################
#Data filters

#Remove lines with values below the LOD
#nrow(data)
#data <- subset(data,copies.reaction >= 4.6 )
#nrow(data)


#Get the average of qPCR replicates
avg.data <- data %>% group_by(Fish_num) %>% group_by(Replicate, .add=TRUE) %>% group_by(Time, .add=TRUE) %>% 
  group_by(Water_sample, .add=TRUE) %>% filter(!is.na(Ct_value), .add=TRUE) %>%
  summarise(avg.copies.liter = mean(copies.liter), sd.copies.liter = sd(copies.liter), n.qPCR.reps = n(), se.copies.liter = sd(copies.liter)/sqrt(n()),Group=mean(Group))

#Filter data to include from a specified time period
#nrow(data)
#avg.data <- subset(avg.data,Time < 72 )
#avg.data <- subset(avg.data,Time >= 12 )
#nrow(data)

#Format output file
header<- c("tank", "decay", "decay.se","t.halflife")
cat(header, file = "decayrates.txt", sep = " ")

#Initialize 
F <- c(5,25) #Number of fish in tank
R <- c(1,2)  #Replicate number

for (f in F) { #open loop for fish number
  for (r in R) { #open loop for replicate number
    #Subset data to include the relevant tank
    data.tank<-avg.data[(avg.data$Fish_num==f & avg.data$Replicate==r) ,]

   plot<-ggplot(data=data.tank, aes(x=Time, y=avg.copies.liter)) +
     geom_smooth(method="lm") +
     geom_point() + 
     stat_regline_equation(label.x=30, label.y=10)+
     stat_cor(aes(label=..rr.label..), label.x=30, label.y=9.5)+
     labs(title=paste0('Tank.',f,'.',r), x ="Hours", y = "Copies per liter") +
     scale_y_continuous(trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma)
   
   assign(paste0('Plot.',f,'.',r), plot, envir=.GlobalEnv)
    
    #Produce a linear model of time and ln(concentration)
    model<-lm( log(data.tank$avg.copies.liter) ~ data.tank$Time )
    summary(model)
    #plot(model) #Plots to evaluate model assumptions

    #The decay rate or k is the slope in units copies per hour
    k<-(model$coefficients[2])
    k # decay rate
    se.k<-sqrt(diag(vcov(model)))[2]
    se.k #standard error of decay rate
    t.halflife<- (0.693/k) #half-life 
    t.halflife
    
    #Write to file
    output <- c(paste0("Tank",f,".",r), -k, se.k,-t.halflife)
    cat("", file = "decayrates.txt", sep = "\n ", append = TRUE)
    cat(output, file = "decayrates.txt", sep = " ", append = TRUE)  # Apply cat & append
   
  } #close loop for tank number
} #close loop for replicate number

decay.rates <- read.table('decayrates.txt',header = TRUE)
decay.rates

#Plots of the decay results
plot_grid(Plot.5.1, Plot.5.2, Plot.25.1, Plot.25.2, align = "hv")


####################################
#Decay Plots
####################################

#reduce margins
margins_left <- c(0.4,-0.2,-0.4,0.2)
margins_right <- c(0.4,0.2,-0.4,0.2)


#5 Fish, Tank 1 - Subset and plot
data.5.1<-avg.data[(avg.data$Fish_num==5 & avg.data$Replicate==1) ,]
plot5.1<-ggplot(data=data.5.1, aes(x=data.5.1$Time, y=data.5.1$avg.copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=10)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=9.5)+
  labs(title=paste0('5 Fish, Tank 1'), x ="Hours", y = "Copies per liter") +
  scale_y_continuous(limits = c(1, 10^5),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(margins_left, "cm"))

#5 Fish, Tank 2 - Subset and plot
data.5.2<-avg.data[(avg.data$Fish_num==5 & avg.data$Replicate==2) ,]
plot5.2<-ggplot(data=data.5.2, aes(x=data.5.2$Time, y=data.5.2$avg.copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=10)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=9.5)+
  labs(title=paste0('5 Fish, Tank 2'), x ="Hours", y = "") +
  scale_y_continuous(limits = c(1, 10^5),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(margins_right, "cm"))

#25 Fish, Tank 1 - Subset and plot
data.25.1<-avg.data[(avg.data$Fish_num==25 & avg.data$Replicate==1) ,]
plot25.1<-ggplot(data=data.25.1, aes(x=data.25.1$Time, y=data.25.1$avg.copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=12)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=11)+
  labs(title=paste0('25 Fish, Tank 1'), x ="Hours", y = "Copies per liter") +
  scale_y_continuous(limits = c(1, 10^6),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(plot.margin=unit(margins_left, "cm"))

#25 Fish, Tank 2 - Subset and plot
data.25.2<-avg.data[(avg.data$Fish_num==25 & avg.data$Replicate==2) ,]
plot25.2<-ggplot(data=data.25.2, aes(x=data.25.2$Time, y=data.25.2$avg.copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=12)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=11)+
  labs(title=paste0('25 Fish, Tank 2'), x ="Hours", y = "") +
  scale_y_continuous(limits = c(1, 10^6), trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(margins_right, "cm"))

plot_grid(plot5.1, plot5.2,plot25.1,plot25.2,align = "hv")


####################################
#ANCOVA to compare slopes/decay rates
####################################

#Prepare the data

#Get the average of qPCR replicates
avg.data <- data %>% group_by(Fish_num) %>% group_by(Replicate, .add=TRUE) %>% group_by(Time, .add=TRUE) %>% 
  group_by(Water_sample, .add=TRUE) %>% filter(!is.na(Ct_value), .add=TRUE) %>%
  summarise(avg.copies.liter = mean(copies.liter), sd.copies.liter = sd(copies.liter), n.qPCR.reps = n(), se.copies.liter = sd(copies.liter)/sqrt(n()),Group=mean(Group))

avg.data$Group <- as.factor(avg.data$Group)
avg.data$ln.avg.copies.liter <- log(avg.data$avg.copies.liter)

#Subsetting the data to include from a specified time period
#nrow(data)
#avg.data <- subset(avg.data,Time < 72 )
#avg.data <- subset(avg.data,Time >= 9 )
#nrow(data)

#Check linearity assumption by making a plot
ggscatter(avg.data, x = "Time", y = "ln.avg.copies.liter", color = "Group", add = "reg.line")+
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Group))

#Test for homogeneity of regression slopes
avg.data <- as.data.frame(avg.data)
avg.data %>% anova_test(ln.avg.copies.liter ~ Group*Time)


#EXTRA CHECKS ON THE DATA

#Normality of residuals
# Fit the model, the covariate goes first
model <- lm(ln.avg.copies.liter ~ Time + Group, data = avg.data)

# Inspect the model diagnostic metrics
model.metrics <- augment(model) 
model.metrics

# Assess normality of residuals using shapiro wilk test
shapiro_test(model.metrics$.resid)

#Homogeneity of variances
#ANCOVA assumes that the variance of the residuals is equal for all groups. This can be checked using the Levene’s test:
model.metrics %>% levene_test(.resid ~ Group)

#Outliers?
model.metrics %>% 
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

#Computation
mod1 <- aov(ln.avg.copies.liter ~ Time*Group, data = avg.data)
summary(mod1)

####################################
#END Decay Rate
####################################


####################################
#Shed Rate
####################################

#Data filters (Note, must include time 0 to execute calculations.)

#Remove lines with values below the LOD
#nrow(data)
#data <- subset(data,copies.reaction >= 4.6 )
#nrow(data)

#Initialize the data
#Get the average of qPCR replicates
avg.data <- data %>% group_by(Fish_num) %>% group_by(Replicate, .add=TRUE) %>% group_by(Time, .add=TRUE) %>% 
  group_by(Water_sample, .add=TRUE) %>% filter(!is.na(Ct_value), .add=TRUE) %>%
  summarise(avg.copies.liter = mean(copies.liter), sd.copies.liter = sd(copies.liter), n.qPCR.reps = n(), se.copies.liter = sd(copies.liter)/sqrt(n()),Group=mean(Group))


#Format output file
header<- c("c.zero", "c.zero.se","s", "s.error", "s.fish", "s.fish.error", "s.g", "s.g.error")
cat(header, file = "shedrates.txt", sep = " ")

#Read in weight data.  Total weight (g) of fish in each tank.
tank.weights<-read.table('tank.weights.txt', header = TRUE)

#Initialize counters
F <- c(5,25) #Fish number
R <- c(1,2) #Replicate

for (f in F) { #open loop for fish number
  for (r in R) { #open loop for replicate number
    #Subset data to include the relevant tank
    data.tank<-avg.data[(avg.data$Fish_num==f & avg.data$Replicate==r) ,]
    
    #dteady state concentration, or C0, units: copies per liter
    c.zero  <- mean(data.tank[(data.tank$Time==0),]$avg.copies.liter)
    c.zero.se <- c.zero/sqrt(length(data.tank[(data.tank$Time==0),]$avg.copies.liter)-1)

    #decay rate, k, units: per hour.  Taken from previous calculations (see above).
    k <- decay.rates[(decay.rates$tank==paste0("Tank",f,".",r)),2]
    k.se <- decay.rates[(decay.rates$tank==paste0("Tank",f,".",r)),3]
    
    #shedding rate, s = k c v, for the tank (copies per hour)
    s <- k * c.zero * 300 #the volume of the experimental tanks was fixed at 300 liters
    s.error <- (sqrt((k.se/k)^2+(c.zero.se/c.zero)^2))*s #error propagated from decay and steady state
    
    #shedding rate per fish
    s.fish <- s/f
    s.fish.error <- (sqrt((k.se/k)^2+(c.zero.se/c.zero)^2))*(s/f)
    
    #shedding rate per gram
    s.g <-s/tank.weights[(tank.weights$tank==paste0("Tank",f,".",r)),2] #shedding rate per gram
    s.g.error <- (sqrt((k.se/k)^2+(c.zero.se/c.zero)^2))*(s/tank.weights[(tank.weights$tank==paste0("Tank",f,".",r)),2])
    
    #Write to file
    output <- c(paste0("Tank",f,".",r), c.zero, c.zero.se, s, s.error, s.fish, s.fish.error, s.g, s.g.error)
    cat("", file = "shedrates.txt", sep = "\n ", append = TRUE)
    cat(output, file = "shedrates.txt", sep = " ", append = TRUE)  # Apply cat & append
    
  } #close loop for tank number
} #close loop for replicate number

#Read and print shed rates
shed.rates <- read.table('shedrates.txt', header = TRUE)
shed.rates

#combine and print both shed and decay rates
cbind(decay.rates, shed.rates)
      
####################################
#END Shed Rate
####################################







####################################
####################################
#Exploratory plots
####################################
####################################

#Temperature vs time
data$Group <- as.factor(data$Group)
ggscatter(data, x = "Time", y = "Temperature", color = "Group")


#Decay Plots
#These plots include concentrations estimated from each replicate qPCR as an independent observation.  Overall, this produces similar results
#to the above plots where qPCR replicates are averaged.


#reduce margins
margins_left <- c(0.4,-0.2,-0.4,0.2)
margins_right <- c(0.4,0.2,-0.4,0.2)


#5 Fish, Tank 1 - Subset and plot
data.5.1<-data[(data$Fish_num==5 & data$Replicate==1) ,]
plot5.1<-ggplot(data=data.5.1, aes(x=Time, y=copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=10)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=9.5)+
  labs(title=paste0('5 Fish, Tank 1'), x ="Hours", y = "Copies per liter") +
  scale_y_continuous(limits = c(1, 10^5),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(margins_left, "cm"))

#5 Fish, Tank 2 - Subset and plot
data.5.2<-data[(data$Fish_num==5 & data$Replicate==2) ,]
plot5.2<-ggplot(data=data.5.2, aes(x=Time, y=copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=10)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=9.5)+
  labs(title=paste0('5 Fish, Tank 2'), x ="Hours", y = "") +
  scale_y_continuous(limits = c(1, 10^5),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x =element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(margins_right, "cm"))

#25 Fish, Tank 1 - Subset and plot
data.25.1<-data[(data$Fish_num==25 & data$Replicate==1) ,]
plot25.1<-ggplot(data=data.25.1, aes(x=Time, y=copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=12)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=11)+
  labs(title=paste0('25 Fish, Tank 1'), x ="Hours", y = "Copies per liter") +
  scale_y_continuous(limits = c(1, 10^6),trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(plot.margin=unit(margins_left, "cm"))

#25 Fish, Tank 2 - Subset and plot
data.25.2<-data[(data$Fish_num==25 & data$Replicate==2) ,]
plot25.2<-ggplot(data=data.25.2, aes(x=Time, y=copies.liter)) +
  geom_smooth(method="lm") +
  geom_point() + 
  stat_regline_equation(label.x=10, label.y=12)+
  stat_cor(aes(label=..rr.label..), label.x=10, label.y=11)+
  labs(title=paste0('25 Fish, Tank 2'), x ="Hours", y = "") +
  scale_y_continuous(limits = c(1, 10^6), trans = log_trans(), breaks = 10^seq(-1,5), labels = scales::comma) +
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(margins_right, "cm"))

plot_grid(plot5.1, plot5.2,plot25.1,plot25.2,align = "hv")