#-------------------------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This script creates the function used in the ode() function. This function is written with all 
# the variables matching those required by the deSolve functions. It also returns a list of the 
# individuals in each age class, as well as the total population. Included in this script is an
# event function to prevent the population falling below 0. All of the variable names match those
# created in the `input_variables.R`. More dynamic varaibles such as predation_level and
# growth_type should be defined in the analysis script calling this function.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#
model <- function(Times, State, Pars){
  with(as.list(c(State, Pars)), {
    # Calculate number of individuals dying
    # First need to know the population number
    population <- JUVENILES + SUB_ADULTS_1 + SUB_ADULTS_2 + ADULTS
    
    # Calculate HIM rate
    if (growth_type == 1 | growth_type == 2){ # none or static growth
      juvenile_HIM <- juvenile_him[growth_type]
      s_adult_HIM <- s_adult_him[growth_type]
      adult_HIM <- adult_him[growth_type]
      # Don't need to worry about egg as there is no additional HIM for this age class
    }
    
    if (growth_type == 3){ #Stochastic human growth
      juvenile_HIM <- rnorm(1, 0.005, 0.0025)
      s_adult_HIM <- rnorm(1, 0.005, 0.0025)
      adult_HIM <- rnorm(1, 0.005, 0.0025)
    }
    
    if (growth_type == 4){ #Linear human growth
      extra_HIM <- Times * 3e-04 # This value is the current HIM multiplied by Stats NZ estimate of
                                 # NZ growth rate
      
      juvenile_HIM <- HIM + extra_HIM
      s_adult_HIM <- HIM + extra_HIM
      adult_HIM <- HIM + extra_HIM
      
    }
    
    if (growth_type == 5){ #Exponenial human growth
      extra_HIM <- Times * (0.01 * HIM)
      
      juvenile_HIM <- (extra_HIM + HIM)
      s_adult_HIM <- (extra_HIM + HIM)
      adult_HIM <- (extra_HIM + HIM)
    }
    
    if (growth_type == 6){ #Logistic human growth
      extra_HIM <- Times * (0.005 * HIM * (1 - (HIM / 0.5)))
      
      juvenile_HIM <- (HIM + extra_HIM)
      s_adult_HIM <- (HIM + extra_HIM)
      adult_HIM <- (HIM + extra_HIM)
    }
    
    if (juvenile_HIM > 1 | s_adult_HIM > 1 | adult_HIM > 1){
      juvenile_HIM <- 1
      s_adult_HIM <- 1
      adult_HIM <- 1
    }
    
    # Check if additional mast predation is needed
    tmp_val <- runif(1)
    mast_pred <- ifelse(tmp_val < mast_prob, mast_predation_bump, 0)
    mast_previous_season <- ifelse(tmp_val < mast_prob, TRUE, FALSE)
    
    # Next we need the death rates
    # If statement checks that death rate isn't greater than 1
    juvenile_death_rate <- ifelse((juvenile_predation[predation_level]
                                   + juvenile_HIM
                                   + mast_pred) > 1, 
                                  1,
                                  (juvenile_predation[predation_level]
                                   + juvenile_HIM
                                   + mast_pred))
    
    s_adult_death_rate <- ifelse((s_adult_predation[predation_level]
                                  + s_adult_HIM
                                  + mast_pred) > 1, 
                                 1,
                                 (s_adult_predation[predation_level]
                                  + s_adult_HIM
                                  + mast_pred))
    
    adult_death_rate <- ifelse((adult_predation[predation_level]
                                + adult_HIM
                                + mast_pred) > 1, 
                               1,
                               (adult_predation[predation_level]
                                + adult_HIM
                                + mast_pred))
    
    # Now we need to calculate the survival rate
    egg_survival <- ifelse(1 - (egg_predation[predation_level] +
                                  egg_him[growth_type] +
                                  mast_pred) > 0,
                           1 - (egg_predation[predation_level] +
                                  egg_him[growth_type] +
                                  mast_pred),
                           0)
    
    juvenile_survival <- ((1 - juvenile_death_rate) * (1 - (population / K)))
    
    sub_adult_survival <- ((1 - s_adult_death_rate) * (1 - (population / K)))
    
    # Now calculate how many individuals actually die
    juvenile_dead <- JUVENILES * (1 - juvenile_survival)
    s_adult1_dead <- SUB_ADULTS_1 * (1 - sub_adult_survival)
    s_adult2_dead <- SUB_ADULTS_2 * (1 - sub_adult_survival)
    adult_dead <- ADULTS * adult_death_rate
    
    # Calculate number of individuals moving to next class
    to_s_adult1 <- JUVENILES * juvenile_survival
    to_s_adult2 <- SUB_ADULTS_1 * sub_adult_survival
    to_adult <- SUB_ADULTS_2 * sub_adult_survival
    
    remove_excess <- ifelse(ADULTS >= K,
                            (ADULTS - K), #To many
                            0) # < K
    
    # Calculate fecundity
    clutch_size <- rpois(1, clutch_size_lambda)
    if(clutch_size == 0) clutch_size <- 1L
    if(clutch_size > 8) clutch_size <- 8L
    numbers_breeding <- ADULTS * prop_breeding
    fecundity <- numbers_breeding * egg_survival * clutch_size
    
    # Calculate class derivative for time step
    d_juveniles <- fecundity - juvenile_dead - to_s_adult1
    d_s_adult1 <- to_s_adult1 - s_adult1_dead - to_s_adult2
    d_s_adult2 <- to_s_adult2 - s_adult2_dead - to_adult
    d_adult <- to_adult - adult_dead - remove_excess
    
    # Check if the population is extinct
    not_extinct <- population > 0
    not_quasi_extinct <- population > 50
    
    # Output values
    # The returned list must be in the same order as y_initial
    return(list(c(d_juveniles, d_s_adult1, d_s_adult2, d_adult), population, 
                predation_level, growth_type, not_extinct, not_quasi_extinct))
  })
}


#~~~Functions to prevent negative populations~~~#
###################################################################################################
# This function can give exception commands to the function. So currently it checks if any age
# class drops below one and if it does it sets that class to 0 (i.e. prevents negative values).

posfun <- function(t, y, parms){
  y[which(y < 1)] <- 0
  return(y)
}
