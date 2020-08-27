#-------------------------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This script analyses the data output from the PVA model run using `Run_analysis.R`. The script
# also produces the visualisations used within the manuscript written for this PVA. The analysis
# is broken into two parts. First only the predation level effects are explored, using the no
# human-growth (i.e. no human induced mortality) runs. The second part of the script analyses the
# the effect of the different human growths, primarily on the baseline predation level runs to 
# see if there is any sort of tipping effect.

# Author: C.E. Simpkins

# Date: 25/04/2020
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

#~~SETUP~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

## Load all the required packages
library("progress") # get a standard progress bar
library("esvis") # calculate hedges_g
library("ggpubr") # help create well organised plots
library("Hmisc") # summary statistics
library("psych") # more summary statistics
library("tidyverse") # To help with data manipulation and plotting
library("wesanderson") # colours
library("plotrix") # improved plotting
library("data.table") # improved speed for rbind
library("foreach") # To create a for loop that can be parallelized
library("doParallel") # to run loops in parallel
library("snow") # To register the backend for the parallel loop

##~~Import data~~##
#-------------------------------------------------------------------------------------------------#
## Import all the data created in `Run_analysis.R` and perform some simple data cleaning to make
## the other analyses easier

## Scenario data
scen_df <- readRDS("Output/Pred_simulations.RData")

## Sensitivity analysis data
SA_baseline <- readRDS("Output/Sensitivity_Analysis/baseline.RData") %>% 
  mutate(parameter = "base")


SA_clutch_plus_25 <- readRDS("Output/Sensitivity_Analysis/clutch_size_plus25.RData") %>% 
  mutate(parameter = "clutch_size_plus25")

SA_clutch_plus_10 <- readRDS("Output/Sensitivity_Analysis/clutch_size_plus10.RData") %>% 
  mutate(parameter = "clutch_size_plus10")

SA_clutch_minus_10 <- readRDS("Output/Sensitivity_Analysis/clutch_size_minus10.RData") %>% 
  mutate(parameter = "clutch_size_minus10")

SA_clutch_minus_25 <- readRDS("Output/Sensitivity_Analysis/clutch_size_minus25.RData") %>% 
  mutate(parameter = "clutch_size_minus25")


SA_breeding_plus_25 <- readRDS("Output/Sensitivity_Analysis/pop_breeding_plus25.RData") %>% 
  mutate(parameter = "breeding_plus25")

SA_breeding_plus_10 <- readRDS("Output/Sensitivity_Analysis/pop_breeding_plus10.RData") %>% 
  mutate(parameter = "breeding_plus10")

SA_breeding_minus_10 <- readRDS("Output/Sensitivity_Analysis/pop_breeding_minus10.RData") %>% 
  mutate(parameter = "breeding_minus10")

SA_breeding_minus_25 <- readRDS("Output/Sensitivity_Analysis/pop_breeding_minus25.RData") %>% 
  mutate(parameter = "breeding_minus25")


SA_pop_plus_25 <- readRDS("Output/Sensitivity_Analysis/pop_size_plus25.RData") %>% 
  mutate(parameter = "pop_plus25")

SA_pop_plus_10 <- readRDS("Output/Sensitivity_Analysis/pop_size_plus10.RData") %>% 
  mutate(parameter = "pop_plus10")

SA_pop_minus_10 <- readRDS("Output/Sensitivity_Analysis/pop_size_minus10.RData") %>% 
  mutate(parameter = "pop_minus10")

SA_pop_minus_25 <- readRDS("Output/Sensitivity_Analysis/pop_size_minus25.RData") %>% 
  mutate(parameter = "pop_minus25")

## Bind all the sensitivity analyses into one data frame
SA_df <- dplyr::bind_rows(list(SA_baseline,
                               SA_clutch_plus_25, SA_clutch_plus_10, SA_clutch_minus_10, SA_clutch_minus_25,
                               SA_breeding_plus_25, SA_breeding_plus_10, SA_breeding_minus_10, SA_breeding_minus_25,
                               SA_pop_plus_25, SA_pop_plus_10, SA_pop_minus_10, SA_pop_minus_25))

## Clean up the global environment
rm(SA_baseline,
   SA_clutch_plus_25, SA_clutch_plus_10, SA_clutch_minus_10, SA_clutch_minus_25,
   SA_breeding_plus_25, SA_breeding_plus_10, SA_breeding_minus_10, SA_breeding_minus_25,
   SA_pop_plus_25, SA_pop_plus_10, SA_pop_minus_10, SA_pop_minus_25)

##~~Data cleaning~~##
#-------------------------------------------------------------------------------------------------#
## There are some idiosyncrasies in the raw data which may impact any analyses we perform on 
## them. So in this section we just clean these up

## There are a few places where the population or stage classes go below zero for a single step.
## This is a result of the non-negative catch only resetting once a negative has occurred. So
## while this wouldn't impact the results of the model it, may impact any summary statistics we
## produce and therefore we just convert any negative population to 0
scen_df[scen_df < 0] <- 0
SA_df[SA_df < 0] <- 0

## We chose to use numeric IDs to identify the different predation levels and growth types. These
## numerics may get a little confusing and annoying going forward so we change them into actual
## character descriptors here
scen_df$growth_type <- recode(scen_df$growth_type, 
                              `1` = "none", `2` = "static", `4` = "linear", `6` = "logistic")

scen_df$predation_level <- recode(scen_df$predation_level, 
                                  `1` = "low", `2` = "baseline", `3` = "high")

## As the first half of the analysis explores the impacts of predation alone, we create a new
## filtered data frame to use in that section of the script
pred_df <- scen_df %>% 
  filter(growth_type == "none", time >= 1)

## As the second half of the analysis explores the relative impact of HIM to look for tipping
## points or similar, and since the high predation level shows such a dramatic rate of extinction
## even with no HIM, we create data frame with high predation level removed. To keep the messaging
## of the paper focussed we decided to not examine logistic growth opting for a more simple story 
## of comparing no HIM, a static HIM, and one type of increasing HIM (linear)
HIM_df <- scen_df %>% 
  filter(predation_level != "high", time >= 1, growth_type != "logistic")


#~~PREDATOR LEVEL ANALYSIS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#
## This section of the script explores the impacts of predation level, when there is no human-
## induced mortality. 

##~~Time series of population~~##
#-------------------------------------------------------------------------------------------------#

## Create a data frame of summary statistics grouped by predation level and time-step (so 250 
## summary points per predation level) to use in the time-series
pred_ts_sum <- pred_df %>% 
  group_by(predation_level, as.integer(time)) %>% 
  summarise(mean_pop = mean(Population, na.rm = TRUE),
            std_pop = sd(Population, na.rm = TRUE),
            stderr_pop = std_pop/ sqrt(length(Population)),
            median_pop = median(Population, na.rm = TRUE),
            pop_75 = list(quantile(Population, 0.75))[[1]][[1]],
            pop_25 = list(quantile(Population, 0.25))[[1]][[1]],
            pop_05 = list(quantile(Population, 0.05))[[1]][[1]],
            pop_95 = list(quantile(Population, 0.95))[[1]][[1]],
            pop_33 = list(quantile(Population, 0.33))[[1]][[1]],
            pop_66 = list(quantile(Population, 0.66))[[1]][[1]]) %>% 
  rename(time = "as.integer(time)")

## Create a time series plots showing the change in median population at each time step for each
## predation level. A ribbon showing the area between the 25th and 75th percentiles is included to 
## indicate variance

palette <- c("#334FDA", "#DF0011", "#21A900")

pred_ts <- ggplot(pred_ts_sum, 
                  aes(y = median_pop, x = time, group = as.factor(predation_level))) +
  geom_ribbon(aes(ymin = pop_33, 
                  ymax = pop_66, 
                  alpha = as.factor(predation_level)), fill = "lightgray") +
  scale_alpha_manual(values = rep(1, 3)) +
  #geom_line(aes(linetype = as.factor(predation_level)), size = 1) +
  #scale_linetype_manual(values = c("dashed", "dotted", "solid")) +
  geom_line(aes(colour = as.factor(predation_level)), size = 1) +
  scale_color_manual(values = palette) +
  ylab("number of individuals") + 
  theme_classic() + theme(legend.position = "none")

pred_ts

##~~Calculate and visualise average population~~##
#-------------------------------------------------------------------------------------------------#

## Again we create a set of summary statistics to help with plotting. Here we first calculate a 
## weighted mean accounting for the fact that any run that went extinct may be artificially 
## inflated due to a far reduced number of ticks in the model run. We then take the median of the
## weighted mean to get the summary statistics for each predation level.
pred_pop_sum <- pred_df %>% 
  group_by(Run, predation_level) %>% 
  summarise(weighted_mean_pop = sum(Population) / 250) %>% 
  group_by(predation_level) %>% 
  summarise(y0 = min(weighted_mean_pop),
            y25 = quantile(weighted_mean_pop, 0.01),
            y50 = quantile(weighted_mean_pop, 0.5),
            y75 = quantile(weighted_mean_pop, 0.99),
            y100 = max(weighted_mean_pop),
            average = mean(weighted_mean_pop),
            stDev = sd(weighted_mean_pop))

pred_pop_plot <- ggplot(pred_pop_sum, aes(x = as.factor(predation_level), group = predation_level)) +
  geom_boxplot(
    aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
    stat = "identity") + 
  ylab("number of individuals") + xlab("scenario") +
  theme_classic()

pred_pop_plot

## Calculate the effect size differences in average population between the predation levels
temp <- pred_df %>% 
  group_by(Run, predation_level) %>% 
  summarise(weighted_mean_pop = sum(Population) / 250)

effect_pop_p <- esvis::hedg_g(weighted_mean_pop ~ predation_level, data = temp)
effect_pop_p$hedg_g <- abs(effect_pop_p$hedg_g)

rm(temp)

##~~Calculate the total number of quasi-extinctions~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how many total extinctions occurred for firstly overall
## (i.e. all no human induced mortality runs), and then for each predation level. We also calculate
## the total number of full extinctions (i.e. pop = 0) to allow for a comparison

repeats <- length(unique(pred_df$Run))

extinct_pred_total <- pred_df %>%
  filter(time == 250) %>% 
  summarise(quasi_extinctions = repeats - sum(not_quasi_extinct),
            quasi_ex_percentage = (quasi_extinctions / repeats) * 100,
            extinctions = repeats - sum(not_extinct),
            ex_percentage = (extinctions / repeats) * 100)

#! Note that the amount of quasi-extinctions are identical to the amount of extinctions, likely
## indicating that if there is sufficiently bad conditions to cause a quasi-extinction in 250 years
## it will result in a full extinction quickly.

# Now we examine how many extinctions occurred under each predation_level
levels <- 3 # there are three predation levels

extinct_pred <- pred_df %>% 
  filter(time == 250) %>% 
  group_by(predation_level) %>%
  summarise(extinctions = (repeats / levels) - sum(not_quasi_extinct),
            percentage = (extinctions / (repeats / levels)) * 100)


##~~Calculate the average time until quasi-extinctions~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how long it took any run which went extinct to go
## extinct

time_extinct_pred <- pred_df %>% 
  group_by(Run, predation_level) %>%
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250) %>% 
  group_by(predation_level) %>% 
  summarise(av_time_extinct = mean(time_till_extinct), sd_time_extinct = sd(time_till_extinct))
  


##~~Calculate the average growth rate~~##
#-------------------------------------------------------------------------------------------------#

## Firstly, we calculate the growth rate at each time step for each run using a simple for loop
## This is done in parallel as it takes a little while otherwise.
registerDoParallel(detectCores() - 1)

r_pred <- as.data.frame(NULL)
r_pred <- foreach(i=1:length(unique(pred_df$Run)),
                  .packages = "data.table",
                  .combine = "rbind") %dopar% {
                    zero_t_pop <- 500 # Note starting number of individuals is hard coded here
                    temp <- pred_df[pred_df$Run == i,]
                    temp$R <- (temp$Population / zero_t_pop) ^ (1 / temp$time)
                    return(temp)
                  }

## Again we create a set of summary statistics to help with plotting. Here we first calculate a 
## weighted mean accounting for the fact that any run that went extinct may be artificially 
## inflated due to a far reduced number of ticks in the model run. We then take the median of the
## weighted mean to get the summary statistics for each predation level.
growth_pred_sum <- ungroup(r_pred) %>% 
  group_by(Run, predation_level) %>% 
  summarise(mean_growth = mean(R, na.rm = TRUE), 
            weighted_mean_growth = sum(R) / 250)%>% 
  group_by(predation_level) %>% 
  summarise(y0 = min(weighted_mean_growth),
            y25 = quantile(weighted_mean_growth, 0.01),
            y50 = quantile(weighted_mean_growth, 0.5),
            y75 = quantile(weighted_mean_growth, 0.99),
            y100 = max(weighted_mean_growth),
            mean_av = mean(weighted_mean_growth),
            stDev = sd(weighted_mean_growth))

## Calculate the effect size differences in growth rate between the predation levels
effect_growth_p <- esvis::hedg_g(R ~ predation_level, data = r_pred)
effect_growth_p$hedg_g <- abs(effect_growth_p$hedg_g)


#~~HUMAN-INDUCED MORTALITY ANALYSIS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#
## This section of the script explores the impacts of human-induced mortailty and changes in its
## magnitude resulting from changes in the human population. We exclude the high predation level
## runs as these have been shown to already be highly impacted by mortality and obscure differences
## which may occur at the other levels.

##~~Time series of population~~##
#-------------------------------------------------------------------------------------------------#
## Create a data frame of summary statistics grouped by predation level, growth type, and time-step 
## (so 250 summary points per predation level and growth type) to use in the time-series
HIM_ts_sum <- HIM_df %>% 
  group_by(predation_level, growth_type, as.integer(time)) %>% 
  summarise(mean_pop = mean(Population, na.rm = TRUE),
            std_pop = sd(Population, na.rm = TRUE),
            stderr_pop = std_pop/ sqrt(length(Population)),
            median_pop = median(Population, na.rm = TRUE),
            pop_75 = list(quantile(Population, 0.75))[[1]][[1]],
            pop_25 = list(quantile(Population, 0.25))[[1]][[1]],
            pop_05 = list(quantile(Population, 0.05))[[1]][[1]],
            pop_95 = list(quantile(Population, 0.95))[[1]][[1]],
            pop_33 = list(quantile(Population, 0.33))[[1]][[1]],
            pop_66 = list(quantile(Population, 0.66))[[1]][[1]],
            pop_45 = list(quantile(Population, 0.45))[[1]][[1]],
            pop_55 = list(quantile(Population, 0.55))[[1]][[1]]) %>% 
  rename(time = "as.integer(time)")

## Create a time series plots showing the change in median population at each time step for each
## growth type. Two separate plots are made, one for low predation level and one for baseline
## predation, these are then presented in a two panel plot with equal axis. A ribbon showing the
## area between the 25th and 75th percentiles is included to indicate variance
## First low predation level
HIM_low_ts <- ggplot(HIM_ts_sum[HIM_ts_sum$predation_level == "low",], 
                     aes(y = median_pop, x = time, group = as.factor(growth_type))) +
  geom_ribbon(aes(ymin = pop_33, 
                  ymax = pop_66, 
                  alpha = as.factor(growth_type)), fill = "gray") +
  scale_alpha_manual(values = rep(0.4, 3)) +
  geom_line(aes(colour = as.factor(growth_type)), size = 0.8) +
  scale_colour_manual(values = c("#DF0011", "#21A900", "#334FDA")) +
  ylab("number of individuals") + ylim(0, 575) + # Note axis does not begin at 0
  theme_classic() + theme(legend.position = "none")

HIM_low_ts

## Second baseline predation level
HIM_base_ts <- ggplot(HIM_ts_sum[HIM_ts_sum$predation_level == "baseline",], 
                      aes(y = median_pop, x = time, group = as.factor(growth_type))) +
  geom_ribbon(aes(ymin = pop_33, 
                  ymax = pop_66, 
                  alpha = as.factor(growth_type)), fill = "gray") +
  scale_alpha_manual(values = rep(0.4, 3)) +
  geom_line(aes(colour = as.factor(growth_type)), size = 0.8) +
  scale_colour_manual(values = c("#DF0011", "#21A900", "#334FDA")) +
  ylab("number of individuals") + ylim(0, 575) + # Note axis does not begin at 0
  theme_classic() + theme(legend.position = "none")

HIM_base_ts

ggarrange(HIM_low_ts + rremove("legend"), 
          HIM_base_ts + rremove("ylab") + rremove("y.text") + rremove("legend"),
          labels = c("A", "B"),
          ncol = 2, nrow = 1, hjust = -32)


##~~Calculate and visualise average population~~##
#-------------------------------------------------------------------------------------------------#

## Again we create a set of summary statistics to help with plotting. Here we first calculate a 
## weighted mean accounting for the fact that any run that went exticnt may be artificially 
## inflated due to a far reduced number of ticks in the model run. 
HIM_pop_sum <- HIM_df %>% 
  group_by(Run, predation_level, growth_type) %>% 
  summarise(weighted_mean_pop = sum(Population) / 250)

HIM_pop_sum$growth_type <- factor(HIM_pop_sum$growth_type,
                                  levels = c("none", "static", "linear"))

HIM_pop_plot <- ggplot(HIM_pop_sum, 
                       aes(x = predation_level, y = weighted_mean_pop, colour = growth_type,
                           fill = growth_type)) +
  geom_violin(scale = "count", trim = FALSE) + 
  geom_boxplot(width = 0.5, position=position_dodge(0.9), aes(fill = "white"),
               outlier.size = 0.5) +
  scale_colour_manual(values = rep("black", 3)) +
  scale_fill_manual(values = c("#DF0011", "#21A900", "#334FDA", "white")) +
  ylim(0, 580) + ylab("mean number of individuals") + xlab("predation level") +
  theme_classic() + theme(legend.position = "none")

HIM_pop_plot

HIM_pop_summary <- HIM_pop_sum %>% 
  group_by(growth_type, predation_level) %>% 
  summarise(average = mean(weighted_mean_pop),
            stDev = sd(weighted_mean_pop))

## Calculate the effect size differences in average population between the predation levels
temp <- HIM_df %>% 
  group_by(Run, predation_level, growth_type) %>% 
  summarise(weighted_mean_pop = sum(Population) / 250)

effect_pop_HIM <- esvis::hedg_g(ungroup(temp), weighted_mean_pop ~ growth_type + predation_level)
effect_pop_HIM$hedg_g <- abs(effect_pop_HIM$hedg_g)

rm(temp)

##~~Calculate the total number of quasi-extinctions~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how many total extinctions occurred for each growth 
## level.

repeats <- length(unique(HIM_df$Run))
levels <- 6 # there are three growth types and two predation_levels

extinct_HIM <- HIM_df %>% 
  filter(time == 250) %>% 
  group_by(predation_level, growth_type) %>%
  summarise(extinctions = (repeats / levels) - sum(not_quasi_extinct),
            percentage = (extinctions / (repeats / levels)) * 100)

##~~Calculate the average time until quasi-extinctions~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how long it took any run which went extinct to go
## extinct

time_extinct_HIM <- HIM_df %>% 
  group_by(Run, predation_level, growth_type) %>%
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250) %>% 
  group_by(predation_level, growth_type) %>% 
  summarise(av_time_extinct = mean(time_till_extinct), sd_time_extinct = sd(time_till_extinct))

## Calculate the effect size differences in time until extinction between growth types for each of
## the predation levels
temp <- HIM_df %>% 
  group_by(Run, predation_level, growth_type) %>% 
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250)

effect_et_HIM <- esvis::hedg_g(ungroup(temp), time_till_extinct ~ growth_type + predation_level)
effect_et_HIM$hedg_g <- abs(effect_et_HIM$hedg_g)

rm(temp)

##~~Calculate the average growth rate~~##
#-------------------------------------------------------------------------------------------------#

## Firstly, we calculate the growth rate at each time step for each run using a simple for loop
## This is done in parallel as it takes a little while otherwise. Same as for predation level
registerDoParallel(detectCores() - 1)

r_HIM <- as.data.frame(NULL)
r_HIM <- foreach(i=1:length(unique(HIM_df$Run)),
                 .packages = "data.table",
                 .combine = "rbind") %dopar% {
                   zero_t_pop <- 500 # Note starting number of individuals is hard coded here
                   temp <- HIM_df[HIM_df$Run == i,]
                   temp$R <- (temp$Population / zero_t_pop) ^ (1 / temp$time)
                   return(temp)
                 }

## Again we create a set of summary statistics to help with plotting. Here we first calculate a 
## weighted mean accounting for the fact that any run that went extinct may be artificially 
## inflated due to a far reduced number of ticks in the model run. We then take the median of the
## weighted mean to get the summary statistics for each predation level.
growth_HIM_sum <- ungroup(r_HIM) %>% 
  group_by(Run, predation_level, growth_type) %>% 
  summarise(mean_growth = mean(R, na.rm = TRUE), 
            weighted_mean_growth = sum(R) / 250)%>% 
  group_by(predation_level, growth_type) %>% 
  summarise(y0 = min(weighted_mean_growth),
            y25 = quantile(weighted_mean_growth, 0.01),
            y50 = quantile(weighted_mean_growth, 0.5),
            y75 = quantile(weighted_mean_growth, 0.99),
            y100 = max(weighted_mean_growth),
            mean_av = mean(weighted_mean_growth))

## Calculate the effect size differences in growth rate between the predation levels
effect_growth_HIM <- esvis::hedg_g(R ~ predation_level + growth_type, data = ungroup(r_HIM))
effect_growth_HIM$hedg_g <- abs(effect_growth_HIM$hedg_g)

## Make a visualisation
HIM_growth <- ungroup(r_HIM) %>% 
  group_by(Run, predation_level, growth_type) %>% 
  summarise(weighted_mean_growth = sum(R) / 250)

HIM_growth$growth_type <- factor(HIM_growth$growth_type,
                                 levels = c("none", "static", "linear"))

HIM_growth_plot <- ggplot(HIM_growth, 
                          aes(x = predation_level, y = weighted_mean_growth, colour = growth_type,
                              fill = growth_type)) +
  geom_violin(scale = "count", trim = FALSE) + 
  geom_boxplot(width = 0.5, position=position_dodge(0.9), aes(fill = "white"),
               outlier.size = 0.5) +
  scale_colour_manual(values = rep("black", 3)) +
  scale_fill_manual(values = c("#DF0011", "#21A900", "#334FDA", "white")) +
  ylim(0.98, 1.005) + ylab("mean growth rate") + xlab("predation level") +
  theme_classic() + theme(legend.position = "none") + 
  geom_hline(yintercept = 1, linetype = "dashed")

HIM_growth_plot


#~~SENSITIVITY ANALYSIS~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#
# Based off of population
# Get sensitivity values
SA_sum <- SA_df %>% 
  group_by(parameter) %>% 
  summarise(av_pop = mean(Population, na.rm = TRUE),
            var_pop = sd(Population, na.rm = TRUE))

sensitivity <- function(base_in, base_out, percent_change, new_out){
  new_in <- base_in * percent_change
  
  delta_in <- new_in - base_in
  delta_out <- new_out - base_out
  
  s <- abs((delta_out / base_out) / (delta_in / base_in))
  
  return(s)
}

clutch_plus25 <- sensitivity(1.5, SA_sum$av_pop[1], 1.25, SA_sum$av_pop[9])
clutch_plus10 <- sensitivity(1.5, SA_sum$av_pop[1], 1.1, SA_sum$av_pop[8])
clutch_minus10 <- sensitivity(1.5, SA_sum$av_pop[1], 0.9, SA_sum$av_pop[6])
clutch_minus25 <- sensitivity(1.5, SA_sum$av_pop[1], 0.75, SA_sum$av_pop[7])

breeding_plus25 <- sensitivity(0.6, SA_sum$av_pop[1], 1.25, SA_sum$av_pop[5])
breeding_plus10 <- sensitivity(0.6, SA_sum$av_pop[1], 1.1, SA_sum$av_pop[4])
breeding_minus10 <- sensitivity(0.6, SA_sum$av_pop[1], 0.9, SA_sum$av_pop[2])
breeding_minus25 <- sensitivity(0.6, SA_sum$av_pop[1], 0.75, SA_sum$av_pop[3])

pop_plus25 <- sensitivity(500, SA_sum$av_pop[1], 1.25, SA_sum$av_pop[13])
pop_plus10 <- sensitivity(500, SA_sum$av_pop[1], 1.1, SA_sum$av_pop[12])
pop_minus10 <- sensitivity(500, SA_sum$av_pop[1], 0.9, SA_sum$av_pop[10])
pop_minus25 <- sensitivity(500, SA_sum$av_pop[1], 0.75, SA_sum$av_pop[11])

# Effect size difference (avoid p-values as can just simulate to significance)
# Rule of thumb <=0.2 small >= 0.8 large rest medium
hedges_g_SA <- esvis::hedg_g(Population ~ parameter, data = SA_df)
hedges_g_SA$hedg_g <- abs(hedges_g_SA$hedg_g)

##~~Time until extinction for sensitivity parameters~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how long it took any run which went extinct to go
## extinct

SA_extinct <- SA_df %>% 
  group_by(Run, parameter) %>%
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250) %>% 
  group_by(parameter) %>% 
  summarise(av_time_extinct = mean(time_till_extinct), sd_time_extinct = sd(time_till_extinct))

## Calculate the effect size differences in time until extinction between growth types for each of
## the predation levels
temp <- SA_df %>% 
  group_by(Run, parameter) %>% 
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250)

effect_et_SA <- esvis::hedg_g(ungroup(temp), time_till_extinct ~ parameter)
effect_et_SA$hedg_g <- abs(effect_et_SA$hedg_g)

rm(temp)

#~~EXPLORING HIGH PREDATION FURTHER~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

HP_baseline <- readRDS("Output/Sensitivity_Analysis/High_pred_pop_base.RData") %>% 
  mutate(parameter = "base")

HP_ten <- readRDS("Output/Sensitivity_Analysis/High_pred_pop_plus10.RData") %>% 
  mutate(parameter = "ten")

HP_twenty <- readRDS("Output/Sensitivity_Analysis/High_pred_pop_plus25.RData") %>% 
  mutate(parameter = "twenty")

## Bind all the sensitivity analyses into one data frame
HP_df <- dplyr::bind_rows(list(HP_baseline, HP_ten, HP_twenty))

## Clean up the global environment
rm(HP_baseline, HP_ten, HP_twenty)


##~~Calculate the total number of quasi-extinctions for high predation scenarios~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how many total extinctions occurred for each growth 
## level.

repeats <- length(unique(HP_df$Run))
#levels <- 3 # there are three growth types and two predation_levels

count_extinct_HP <- HP_df %>% 
  filter(time == 250) %>% 
  group_by(parameter) %>%
  summarise(extinctions = repeats - sum(not_quasi_extinct),
            percentage = extinctions / repeats * 100)

##~~Time until extinction for high predation scenarios~~##
#-------------------------------------------------------------------------------------------------#
## We defined an extinction as any point where the population dropped below 50 individuals and
## in this section of the code we calculate how long it took any run which went extinct to go
## extinct

HP_extinct_time <- HP_df %>% 
  group_by(Run, parameter) %>%
  summarise(time_till_extinct = sum(not_quasi_extinct)) %>% 
  filter(time_till_extinct != 250) %>% 
  group_by(parameter) %>% 
  summarise(av_time_extinct = mean(time_till_extinct), sd_time_extinct = sd(time_till_extinct))

