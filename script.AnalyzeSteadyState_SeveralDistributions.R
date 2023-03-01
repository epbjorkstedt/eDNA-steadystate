########################################################################
########################################################################
########################################################################
#Decay rates using more complex models
#Package ‘nlsMicrobio’
#Data sets and nonlinear regression models dedicated to predictive microbiology.
########################################################################
########################################################################
########################################################################

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
require(nlsMicrobio)
require(AICcmodavg)
require(dplyr)
require(ggplot2)
require(tidybayes)

#Set working directory
# setwd("/Users/apk7001/Desktop/Decay")



#Load data
data <-read.csv(file="./Decay.Data.csv", header=TRUE, as.is=TRUE)

#Convert Ct value to copies/reaction in qPCR reaction
#Slope and intercept from qPCR standard curve for O. mykiss
#Equation: 10^(Ct value - intercept)/slope = copies/reaction
slope <- -3.402646
intercept <- 38.63023
data <- data %>%
  dplyr::mutate(copies.reaction = 10^((Ct - intercept)/slope))

#Convert copies per qPCR reaction to copies per liter of water filtered
#X (copies/L) = (Y·b)/(c·a)
#X(copies/L) The eDNA copy number in the original water sample
#a volume of water filtered (in liters)
#b the re-suspended volume of purified DNA (fixed b at 100 µL)
#c the amount of purified DNA added to each qPCR-reaction volume (c: 6 µL/qPCR reaction)
#Y the measured copy number in the qPCR-reaction volume, based on the standard curve (Y, copies/PCR reaction)
data <- data %>%
  dplyr::mutate(copies.liter = (copies.reaction * Elution_Volume)/Template * Volume_filtered_L) %>%
  dplyr::mutate(ln.copies.liter = log(copies.liter)) %>%
  dplyr::mutate(log10.copies.liter = log10(copies.liter))  %>%
  dplyr::mutate(individualWaterSample = paste0(Target,Temperature,Time_Hours,Water_sample))

#Average by qPCR replicate
avg.data <- data %>%
  group_by(Target) %>%
  group_by(Temperature, .add=TRUE) %>%
  group_by(Time_Hours, .add=TRUE) %>%
  group_by(Water_sample, .add=TRUE) %>%
  dplyr::filter(!is.na(Ct), .add=TRUE) %>%
  summarise(avg.copies.liter = mean(copies.liter),
            sd.copies.liter = sd(copies.liter),
            n.qPCR.reps = n(),
            se.copies.liter = sd(copies.liter)/sqrt(n())) %>%
  dplyr::mutate(ln.copies.liter = log(avg.copies.liter)) %>%
  dplyr::mutate(log10.copies.liter = log10(avg.copies.liter)) %>%
  dplyr::mutate(individualWaterSample = Water_sample)

#Pull for ATR at 12 C
data_ATR_12 <- data %>%
  dplyr::filter(Target == "Atr" & Temperature == 12)

#Pull for ATR at 12 C, averaged by qPCR replicate.
avg.data_ATR_12 <- avg.data %>%
  dplyr::filter(Target == "Atr" & Temperature == 12)


######################################
#ATR 12
survival_ATR_12 <- as.data.frame(cbind(avg.data_ATR_12$Time_Hours, avg.data_ATR_12$log10.copies.liter))
survival_ATR_12 <- survival_ATR_12 %>%
  rename("t" = "V1",
         "LOG10N" = "V2")


survivalcurve <- survival_ATR_12

#geeraerd
geeraerd <- nls(geeraerd, survivalcurve,
                list(Sl = 5, kmax = 1.5, LOG10N0 = 7, LOG10Nres = 1))
overview(geeraerd)
plotfit(geeraerd, smooth = TRUE)

#geeraerd_without_Nres
geeraerd_without_Nres <- nls(geeraerd_without_Nres, survivalcurve,
                             list(Sl = 10, kmax = 1.5, LOG10N0 = 7.5))
overview(geeraerd_without_Nres)
plotfit(geeraerd_without_Nres, smooth = TRUE)

#geeraerd_without_Sl
geeraerd_without_Sl <- nls(geeraerd_without_Sl, survivalcurve,
                           list(kmax = 4, LOG10N0 = 7.5, LOG10Nres = 1))
overview(geeraerd_without_Sl)
plotfit(geeraerd_without_Sl, smooth = TRUE)

#mafart
mafart <- nls(mafart, survivalcurve,
              list(p = 1.5, delta = 8, LOG10N0 = 7.5))
overview(mafart)
plotfit(mafart, smooth = TRUE)

#albert
albert <- nls(albert,survivalcurve,
              list(p = 1.2, delta = 4, LOG10N0 = 7, LOG10Nres = 1))
overview(albert)
plotfit(albert, smooth = TRUE)


#trilinear - failed to fit, not included in model selection
trilinear <- nls(trilinear, survivalcurve, list(Sl = 5, kmax = 1.5, LOG10N0 = 7, LOG10Nres = 1))
overview(trilinear)
plotfit(trilinear, smooth = TRUE)

#bilinear_without_Nres - failed to fit, not included in model selection
bilinear_without_Nres <- nls(bilinear_without_Nres, survivalcurve, list(Sl = 10, kmax = 1.5, LOG10N0 = 7.5))
overview(bilinear_without_Nres)
plotfit(bilinear_without_Nres, smooth = TRUE)

#bilinear_without_Sl - failed to fit, not included in model selection
bilinear_without_Sl <- nls(bilinear_without_Sl, survivalcurve,
                           list(kmax = 4, LOG10N0 = 7.5, LOG10Nres = 1))
overview(bilinear_without_Sl)
plotfit(bilinear_without_Sl, smooth = TRUE)

#Bigelow - log linear model fit in nls.  Null model
bigelow <- nls(LOG10N ~ (LOG10N0-kmax*t)/log(10), survivalcurve, start= list(kmax = .05, LOG10N0 = 7.5))
overview(bigelow)
plotfit(bigelow, smooth = TRUE)

#####AICc model selection
#define list of models
models <- list(geeraerd, geeraerd_without_Nres, geeraerd_without_Sl, mafart, albert, bigelow)

#specify model names
mod.names <- c('geeraerd', 'geeraerd_without_Nres', 'geeraerd_without_Sl','mafart', 'albert', 'bigelow')

#calculate AIC of each model
aictab(cand.set = models, modnames = mod.names)