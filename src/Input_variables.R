#-------------------------------------------------------------------------------------------------#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This script is just used to set all the initial and constant values used in the Kea PVA, 
# independent of the primary scripts. This script is then sourced into the other scripts, ensuring
# they are not accidentally structurally changed when values are adjusted. This method of creating
# separate modules of script should also help if we choose to develop an RShiny app.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#-------------------------------------------------------------------------------------------------#

#~~~Set the number of individuals in each age class~~~#
#-------------------------------------------------------------------------------------------------#
JUVENILES <- 125
SUB_ADULTS_1 <- 125
SUB_ADULTS_2 <- 125
ADULTS <- 125
HIM <- 0.015

K <- 1000 # Carrying Capacity (i.e. maximum number of adults)

#~~~Set fecundity values~~~#
#-------------------------------------------------------------------------------------------------#
prop_breeding <- 0.6 # Proportion of females which breed each time step
clutch_size_lambda <- 1.5 # lambda value used in rpois(1,lambda) to generate the clutch size for
                        # each time step of the model run.

#~~~Set human induced mortality values~~~#
#-------------------------------------------------------------------------------------------------#
human_growth <- 0.015

# 1) No human growth values
egg_him1 <- 0
juvenile_him1 <- 0
sub_adult_him1 <- 0
adult_him1 <- 0

# 2) Static human growth values
egg_him2 <- 0
juvenile_him2 <- human_growth
sub_adult_him2 <- human_growth
adult_him2 <- human_growth

# 3) Stochastic human growth values
egg_him3 <- 0
juvenile_him3 <- rnorm(1, 0.005, 0.0025)
sub_adult_him3 <- rnorm(1, 0.005, 0.0025)
adult_him3 <- rnorm(1, 0.005, 0.0025)

# 4) Linear human growth values
egg_him4 <- 0
juvenile_him4 <- (0.0003 + human_growth)
sub_adult_him4 <- (0.0003 + human_growth)
adult_him4 <- (0.0003 + human_growth)

# 5) Exponenial human growth values
egg_him5 <- 0
juvenile_him5 <- ((0.01 * human_growth) + human_growth)
sub_adult_him5 <- ((0.01 * human_growth) + human_growth)
adult_him5 <- ((0.01 * human_growth) + human_growth)

# 6) Logistic human growth values
egg_him6 <- 0
juvenile_him6 <- (human_growth + (0.005 * human_growth * (1 - (human_growth / 0.5))))
sub_adult_him6 <- (human_growth + (0.005 * human_growth * (1 - (human_growth / 0.5))))
adult_him6 <- (human_growth + (0.005 * human_growth * (1 - (human_growth / 0.5))))

# Join them up into a list to use in the actual function
egg_him <- c(egg_him1, egg_him2, egg_him3, egg_him4, egg_him5, egg_him6)
juvenile_him <- c(juvenile_him1, juvenile_him2, juvenile_him3, 
                  juvenile_him4, juvenile_him5, juvenile_him6)
s_adult_him <- c(sub_adult_him1, sub_adult_him2, sub_adult_him3, 
                   sub_adult_him4, sub_adult_him5, sub_adult_him6)
adult_him <- c(adult_him1, adult_him2, adult_him3, adult_him4, adult_him5, adult_him6)

#~~~Set predation rates~~~#
#-------------------------------------------------------------------------------------------------#
mast_prob <- 0.25 # probability of a mast year occuring
mast_predation_bump <- 0.1 # the increase in predation rate occurring in mast periods
mast_previous_season <- FALSE # Counter to indicate whether there was a mast the previous season 

# 1) Low predation (effective predator management)
egg_pred1 <- runif(n = 1, min = 0.18, max = 0.22) #One value drawn from a uniform distribution
juvenile_pred1 <- runif(n = 1, min = 0.09, max = 0.11) 
sub_adult_pred1 <- runif(n = 1, min = 0.02, max = 0.03) 
adult_pred1 <- runif(n = 1, min = 0.045, max = 0.055) 

# 2) Medium predation
egg_pred2 <- runif(n = 1, min = 0.36, max = 0.44) 
juvenile_pred2 <- runif(n = 1, min = 0.18, max = 0.22) 
sub_adult_pred2 <- runif(n = 1, min = 0.045, max = 0.055)  
adult_pred2 <- runif(n = 1, min = 0.09, max = 0.11) 

# 3) High predation (As close to no predator management as we have data for)
egg_pred3 <- runif(n = 1, min = 0.72, max = 0.88) 
juvenile_pred3 <- runif(n = 1, min =0.36, max = 0.44) 
sub_adult_pred3 <- runif(n = 1, min = 0.18, max = 0.22)
adult_pred3 <- runif(n = 1, min = 0.27, max = 0.33) 

# Join them up into a list to use in the actual function
egg_predation <- c(egg_pred1, egg_pred2, egg_pred3)
juvenile_predation <- c(juvenile_pred1, juvenile_pred2, juvenile_pred3)
s_adult_predation <- c(sub_adult_pred1, sub_adult_pred2, sub_adult_pred3)
adult_predation <- c(adult_pred1, adult_pred2, adult_pred3)

#~~~Gather all of the variables together in appropriate lists~~~#
#-------------------------------------------------------------------------------------------------#
# *The order of lists is very important in the deSolve package so need to keep these constistent*
# Additionally, some list elements are named themselves to avoid them becoming anonymous when
# run in the ode() function.

# Initial values
y_initial <- c(JUVENILES = JUVENILES,
               SUB_ADULTS_1 = SUB_ADULTS_1, 
               SUB_ADULTS_2 = SUB_ADULTS_2,
               ADULTS = ADULTS)#,
               #HIM = HIM)

# Parameter (constants and flows)
pars <- c(K, prop_breeding, clutch_size_lambda, mast_prob, mast_previous_season,
          egg_him, juvenile_him, s_adult_him, adult_him, mast_predation_bump,
          egg_predation, juvenile_predation, s_adult_predation, adult_predation)
