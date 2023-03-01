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
require(cmocean)

nIter <- 100
modInterval <- floor(nIter/20)

eDNA_base <- data.frame(time = 0,
                   size_class = seq(from=-2, to=-0.2, by=0.2)) %>% 
  dplyr::mutate(size_class = 10^size_class) %>%     # size classes
  dplyr::mutate(eDNAconc = 5*size_class^0.667) %>%  # eDNA load per particle
  dplyr::mutate(k = (1/size_class)^1/4) %>%         # decay rates
  dplyr::mutate(k = k*0.05) %>%                     # decay rates
  dplyr::mutate(a = 0.15*size_class^(2/3) ) %>%     # degradation rates surface:volume scaling
  dplyr::mutate(N = (1/size_class)^1/4)             # original particle concentrations

eDNA <- eDNA_experiment <- eDNA_base


for (time in seq(nIter + 3)) {
  eDNA <- eDNA %>% 
    dplyr::mutate(N = N*exp(-k)) %>% 
    dplyr::mutate(time = time + 1)
  eDNA_experiment <- eDNA_experiment %>% 
    bind_rows(eDNA)
}


eDNA_filtered <- eDNA_experiment %>% 
  dplyr::mutate(total_eDNA = N*eDNAconc) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarize(total_eDNA = sum(total_eDNA)) %>% 
  dplyr::mutate(log10_total_eDNA = log10(total_eDNA)) %>% 
  dplyr::ungroup()

plotInuse <- ggplot(data = eDNA_filtered, aes(x=time, y=total_eDNA)) + 
  scale_y_log10()

for (t in seq(from=3, to=max(eDNA_filtered$time), by=10)) {
  
  dataInuse <- eDNA_filtered %>% 
    dplyr::filter(time <= t)
  
  mod <- lm(log10_total_eDNA ~ time, dataInuse)
  int <- mod$coefficients[1] %>% as.numeric() 
  slp <- mod$coefficients[2] %>% as.numeric()
  plotData <- data.frame(time=unique(dataInuse$time))
  plotData <- plotData %>% 
    dplyr::mutate(log10_total_eDNA = int + slp * time) %>% 
    dplyr::mutate(total_eDNA = 10^log10_total_eDNA)
  plotInuse <- plotInuse + 
    geom_abline(intercept = int, slope = slp,alpha=0.2) +
    geom_line(data = plotData, aes(x=time,y=total_eDNA),linewidth=0.65)
  
  
}

plotInuse <- plotInuse + 
  geom_point(data = eDNA_filtered, aes(x=time,y=total_eDNA),color='blue',size=2,alpha=0.7) 

plotInuseSimpleDecay <- plotInuse


eDNA_distributions <- eDNA_experiment %>% 
  dplyr::filter(time%%3 == 0)

sizeClass_plot <- ggplot(data = eDNA_distributions) +
  geom_line(aes(x=size_class,y=N,group=time,color=time)) +
  scale_x_log10()

logSizeClass_plot <- sizeClass_plot +
  scale_y_log10() 

sizeClass_plot_SimpleDecay <- sizeClass_plot
logSizeClass_plot_SimpleDecay <- logSizeClass_plot

### include degradation to next smaller size class



eDNA <- eDNA_experiment <- eDNA_base

decayMatrix <- diag(exp(-eDNA_base$k))
degraMatrix <- diag(1-eDNA_base$a)
for (i in 1:(dim(degraMatrix)[1]-1)) {
  degraMatrix[i,i+1] <- eDNA_base$a[i+1]
}





for (time in seq(nIter + 3)) {
  eDNA <- eDNA %>% 
    dplyr::mutate(N = decayMatrix%*%N) %>% 
    dplyr::mutate(N = degraMatrix%*%N) %>% 
    dplyr::mutate(time = time + 1)
  eDNA_experiment <- eDNA_experiment %>% 
    bind_rows(eDNA)
}


eDNA_filtered <- eDNA_experiment %>% 
  dplyr::mutate(total_eDNA = N*eDNAconc) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarize(total_eDNA = sum(total_eDNA)) %>% 
  dplyr::mutate(log10_total_eDNA = log10(total_eDNA)) %>% 
  dplyr::ungroup()

# plotInuse <- ggplot(data = eDNA_filtered, aes(x=time, y=total_eDNA)) + 
#   scale_y_log10()

for (t in seq(from=3, to=max(eDNA_filtered$time), by=10)) {
  
  dataInuse <- eDNA_filtered %>% 
    dplyr::filter(time <= t)
  
  mod <- lm(log10_total_eDNA ~ time, dataInuse)
  int <- mod$coefficients[1] %>% as.numeric() 
  slp <- mod$coefficients[2] %>% as.numeric()
  plotData <- data.frame(time=unique(dataInuse$time))
  plotData <- plotData %>% 
    dplyr::mutate(log10_total_eDNA = int + slp * time) %>% 
    dplyr::mutate(total_eDNA = 10^log10_total_eDNA)
  plotInuse <- plotInuse + 
    geom_abline(intercept = int, slope = slp,alpha=0.2, color='pink') +
    geom_line(data = plotData, aes(x=time,y=total_eDNA),linewidth=0.65, color='red')
  
  
}

plotInuse <- plotInuse + 
  geom_point(data = eDNA_filtered, aes(x=time,y=total_eDNA),color='green',size=2,alpha=0.7) 

plotInuse


eDNA_distributions <- eDNA_experiment %>% 
  dplyr::filter(time%%modInterval == 0)

sizeClass_plot_degra <- ggplot(data = eDNA_distributions) +
  geom_line(aes(x=size_class,y=N,group=time,color=time)) +
  scale_x_log10()

logSizeClass_plot_degra <- sizeClass_plot +
  scale_y_log10() 

sizeClass_barplot_degra <- ggplot(data = eDNA_distributions,aes(fill=size_class,y=N,x=time)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_cmocean(name='haline')

### include degradation to next smaller size class



eDNA <- eDNA_experiment <- eDNA_base

decayMatrix <- diag(exp(-eDNA_base$k))
degraMatrix <- diag(1-eDNA_base$a)

decdegMatrix <- decayMatrix * degraMatrix




for (time in seq(nIter + 3)) {
  eDNA <- eDNA %>% 
    dplyr::mutate(N = decdegMatrix%*%N) %>% 
    # dplyr::mutate(N = degraMatrix%*%N) %>% 
    dplyr::mutate(time = time + 1)
  eDNA_experiment <- eDNA_experiment %>% 
    bind_rows(eDNA)
}


eDNA_filtered <- eDNA_experiment %>% 
  dplyr::mutate(total_eDNA = N*eDNAconc) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarize(total_eDNA = sum(total_eDNA)) %>% 
  dplyr::mutate(log10_total_eDNA = log10(total_eDNA)) %>% 
  dplyr::ungroup()

# plotInuse <- ggplot(data = eDNA_filtered, aes(x=time, y=total_eDNA)) + 
#   scale_y_log10()

for (t in seq(from=3, to=max(eDNA_filtered$time), by=10)) {
  
  dataInuse <- eDNA_filtered %>% 
    dplyr::filter(time <= t)
  
  mod <- lm(log10_total_eDNA ~ time, dataInuse)
  int <- mod$coefficients[1] %>% as.numeric() 
  slp <- mod$coefficients[2] %>% as.numeric()
  plotData <- data.frame(time=unique(dataInuse$time))
  plotData <- plotData %>% 
    dplyr::mutate(log10_total_eDNA = int + slp * time) %>% 
    dplyr::mutate(total_eDNA = 10^log10_total_eDNA)
  plotInuse <- plotInuse + 
    geom_abline(intercept = int, slope = slp,alpha=0.2, color='white') +
    geom_line(data = plotData, aes(x=time,y=total_eDNA),linewidth=0.65, color='yellow')
  
  
}


plotInuse <- plotInuse + 
  geom_point(data = eDNA_filtered, aes(x=time,y=total_eDNA),color='orange',size=2,alpha=0.7) 

plotInuse


eDNA_distributions <- eDNA_experiment %>% 
  dplyr::filter(time%%modInterval == 0)

sizeClass_plot_degra <- ggplot(data = eDNA_distributions) +
  geom_line(aes(x=size_class,y=N,group=time,color=time)) +
  scale_x_log10()

logSizeClass_plot_decdeg <- sizeClass_plot +
  scale_y_log10() 

sizeClass_barplot_decdeg <- ggplot(data = eDNA_distributions,aes(fill=size_class,y=N,x=time)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_cmocean(name='haline')
