#-------------------------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This script creates a function to run the PVA model. This allows this function to easily be 
# called in an apply, or map manner thus allowing multiple scenarios and repeats to be easily run.

# DOES NOT WORK AND MAY NOT BE NEEDED
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

run_PVA <- function(iterations = 100, growth_type, predation_level){
  # Define the dynamic variables unique to each scenario
  iterations <- iterations # number of time steps
  
  growth_type <- growth_type # Something to do with a scenario (1 = None; 2 = static; 3 = stochastic,
  #                                               4 = linear; 5 = exponenial; 6 = logistic)
  predation_level <- predation_level # 1 = low; 2 = medium; 3 = high
  
  # Run the pva function
  out <- data.frame(ode(y = y_initial, 
                        times = seq(0, iterations, by = 1), 
                        func = model, 
                        parms = pars,
                        events = list(func = eventfun, 
                                      times = iterations)))
  
  out <- rename(out, Population = X5, predation_level = X6, Growth_type = X7)
  
  return(out)
}