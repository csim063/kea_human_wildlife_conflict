#-------------------------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This script pulls together all the other scripts used to create the PVA model and actually runs
# them. It then exports the raw results to keep a record before doing some analysis and plotting
# the results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

#~~~Setup~~~#
#-------------------------------------------------------------------------------------------------#
# Import the required libraries and source the required scripts
library(deSolve) # This is the key library that actually pieces the model together
library(tidyverse) # To help with data manipulation and plotting
library(foreach) # To create a for loop that can be parallelized
library(doParallel) # to run loops in parallel
library(doSNOW) # to run loops in parallel and get a parallel progress bar
library(wesanderson) # best plotting colour palette
library(plotrix) # for standard error calculation

source("src/Input_variables.R") # Define the constants
source("src/PVA_Function.R") # Create the model function

# Set the key dymanic values
repeats <- 1000 # How many times the model is run under each scenario 1000
iterations <- 250 # How many time steps each run has 250
predation <- c(1, 2, 3)
growth <- c(rep(1, length(predation)),
            rep(2, length(predation)),
            rep(4, length(predation)))

# Set up how many cores to run in parallel
num_cores <- (detectCores() - 1)

#~~~Run the PVA~~~#
#-------------------------------------------------------------------------------------------------#
# Run the model the required number of times
time_check <- Sys.time()

data_out <- data.frame()

cl <- snow::makeCluster(num_cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = repeats * length(growth), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

data_out <- foreach (i = 1:(repeats * length(growth)),
                     .packages = c("deSolve", "tidyverse"),
                     .combine = rbind,
                     .options.snow = opts) %dopar% {
  iterations <- iterations # number of time steps
  
  growth_type_list <- growth
  
  ifelse((i %% length(growth_type_list)) == 0,
         current_g <- length(growth_type_list),
         current_g <- (i %% length(growth_type_list)))
  
  growth_type <- growth_type_list[current_g]
  #         (1 = None; 2 = static; 3 = stochastic, 4 = linear; 5 = exponenial; 6 = logistic)
  
  
  predation_level_list <- predation
  
  ifelse((i %% length(predation_level_list)) == 0,
         current_p <- length(predation_level_list),
         current_p <- (i %% length(predation_level_list)))
  
  predation_level <- predation_level_list[current_p] 
  
  # Run the pva function (1 = low, 2 = current, 3 = high)
  out <- data.frame(deSolve::ode(y = y_initial,
                        times = seq(0, iterations, by = 1),
                        func = model,
                        parms = pars,
                        events = list(func = posfun,
                                      time = c(0:iterations))))
  
  out <- dplyr::rename(out, Population = X5, predation_level = X6, growth_type = X7, 
                       not_extinct = X8, not_quasi_extinct = X9)
  
  out <- out %>% 
    mutate(Run = i)
  
  
  return(out)
}

close(pb)
stopCluster(cl)

print(Sys.time() - time_check)

saveRDS(data_out, "Output/Pred_simulations.RData")

#~~~Calculate summary statistics~~~#
#-------------------------------------------------------------------------------------------------#
# TODO

# Set the growth_type and predation_level as factors
data_out$growth_type <- as.factor(data_out$growth_type)
data_out$predation_level <- as.factor(data_out$predation_level)

# Summarise
d_sum <- data_out %>% 
  group_by(growth_type, predation_level, time) %>% 
  summarise(av_juv = median(JUVENILES, na.rm = TRUE),
            lower_juv = median(JUVENILES, na.rm = TRUE) - std.error(JUVENILES, na.rm = TRUE),
            upper_juv = median(JUVENILES, na.rm = TRUE) + std.error(JUVENILES, na.rm = TRUE),
            av_sub = median((SUB_ADULTS_1 + SUB_ADULTS_2), na.rm = TRUE),
            lower_sub = median((SUB_ADULTS_1 + SUB_ADULTS_2), na.rm = TRUE) - 
              std.error((SUB_ADULTS_1 + SUB_ADULTS_2), na.rm = TRUE),
            upper_sub = median((SUB_ADULTS_1 + SUB_ADULTS_2), na.rm = TRUE) + 
              std.error((SUB_ADULTS_1 + SUB_ADULTS_2), na.rm = TRUE),
            av_adult = median(ADULTS, na.rm = TRUE),
            lower_adult =  median(ADULTS, na.rm = TRUE) -  std.error(ADULTS, na.rm = TRUE),
            upper_adult =  median(ADULTS, na.rm = TRUE) +  std.error(ADULTS, na.rm = TRUE),
            av_pop = median(Population, na.rm = TRUE),
            lower_pop = median(Population, na.rm = TRUE) - std.error(Population, na.rm = TRUE),
            upper_pop = median(Population, na.rm = TRUE) + std.error(Population, na.rm = TRUE))

# Remove na and non-integer time step rows as these seem to be edge cases which ended early
# Also get rid of time step 0
d_sum <- d_sum %>% 
  na.omit(.) %>% 
  filter(time %% 1 == 0) %>% 
  filter(time != 0)

# Change the names of the growth types to look nicer when plotting
d_sum$growth_type <- plyr::revalue(d_sum$growth_type, c("2" = "Static", 
                                                        "1" = "No growth",
                                                        "4" = "Linear"))

#~~~Plot~~~#
#-------------------------------------------------------------------------------------------------#
# TODO

palette <- wes_palette("Cavalcanti1", 3, type = "discrete")

theme_set(theme_classic())
p <- ggplot(data = d_sum, aes(x = time, y = av_pop, 
                              group = growth_type, color = growth_type,
                              name = "Try")) +
  geom_ribbon(aes(ymin = lower_pop, ymax = upper_pop),
              alpha = 0.2, linetype=2) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = palette, name = "human population \n growth type") +
  xlab("time steps") + ylab("number of individuals")
  

p
n



baseline <- readRDS("Output/SA_Baseline.RData")

breeding_low <- readRDS("Output/SA_Breeding_low.RData")
breeding_high <- readRDS("Output/SA_Breeding_high.RData")

clutch_high <- readRDS("Output/SA_clutch_high.RData")
clutch_low <- readRDS("Output/SA_clutch_low.RData")

baseline <- baseline %>% 
  mutate(d_variable = "Baseline")

breeding_low <- breeding_low %>% 
  mutate(d_variable = "breeding_low")

breeding_high <- breeding_high %>% 
  mutate(d_variable = "breeding_high")

clutch_low <- clutch_low %>% 
  mutate(d_variable = "clutch_low")

clutch_high <- clutch_high %>% 
  mutate(d_variable = "clutch_high")

SA_df <- rbind(baseline, breeding_high, breeding_low, clutch_high, clutch_low)
saveRDS(SA_df, "Output/SA_data.RData")

