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

# Set the key dynamic values
repeats <- 1000 # How many times the model is run under each scenario 1000
iterations <- 250 # How many time steps each run has 250
predation <- 3
growth <- 1

# Set the sensitivity values
JUVENILES <- 125 #125 #94/113/138/157
SUB_ADULTS_1 <- 125 #125 #94/113/138/157
SUB_ADULTS_2 <- 125 #125 #94/113/138/157
ADULTS <- 125 #125 #94/113/138/157

prop_breeding <- 0.6 #0.45/0.54/0.66/0.75
clutch_size_lambda <- 1.5 #1.125/1.35/1.65/1.875

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

saveRDS(data_out, "Output/Sensitivity_Analysis/High_pred_pop_base.RData")