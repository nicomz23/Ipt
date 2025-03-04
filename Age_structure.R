
# PACKAGES
require(deSolve)
require(reshape2)
require(ggplot2)
require(tidyverse)
require(gridExtra)
require(latticeExtra)
require(scales)

## ggplot theme
theme <- theme_set(theme_bw())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))

Nh <- 10000   # total number of hosts
Nv <- 150000   # total number of vectors

## one strain model to get to an equilibrium, and then manually set starting prevalence
## initial state values for the one strain model to run to equilibrium
## will run the model for a while without chemoprevention to see where it stabilises to use as a starting point
initial_state_values <- c(Sh = 0.9 * Nh,
                          Ih = 0.1 * Nh,
                          Sv = Nv,
                          Ev = 0,
                          Iv = 0)

## parameters
parameters <- c(a = 0.3,       # mosquito biting rate per day (~1 bite per three days)
                b_v = 0.2,   # probability of infection from an infected host to a susceptible vector
                b_h = 0.40,   # probability of infection from an infected vector to a susceptible host
                mu_v = 1/14,  # mortality/recruitment rate of the vector
                gamma = 1/8,   # recovery rate from symptomatic treated individuals ~(5-14) days
                epsilon = 1/10) # rate of transition from exposed to infected vectors ~ 10 days

times <- seq(from = 0, to = 365, by = 1)    # simulate for a length of time to find equilibrium

rm_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    Nh <- Sh + Ih     # total human population
    Nv <- Sv + Ev + Iv # total vector population
    
    ## force of infections
    lambda_h <- ((a * b_h)/Nh) * Iv
    lambda_v <- ((a * b_v)/Nh) * Ih
    
    # The differential equations
    # Host population dynamics:
    dSh <- - (lambda_h * Sh) + (gamma * Ih)
    dIh <-   (lambda_h * Sh) - (gamma * Ih)
    
    # Vector population dynamics:
    dSv <- (mu_v * Nv) - (lambda_v * Sv) - (mu_v * Sv)         
    dEv <- (lambda_v * Sv) - (epsilon * Ev) - (mu_v * Ev)
    dIv <- (epsilon * Ev) - (mu_v * Iv)
    
    # Output
    return(list(c(dSh, dIh, dSv, dEv, dIv))) 
  })
}

# Solve the equations
output_one <- as.data.frame(ode(y = initial_state_values, 
                                times = times, 
                                func = rm_model,
                                parms = parameters))


# Solve the equations
output_one <- as.data.frame(ode(y = initial_state_values, 
                                times = times, 
                                func = rm_model,
                                parms = parameters))

# Solve the equations
output_one <- as.data.frame(ode(y = initial_state_values, 
                                times = times, 
                                func = rm_model,
                                parms = parameters))

equilibrium <- output_one %>% dplyr::filter(time == 210)

equilibrium_start <- c(Sh = equilibrium$Sh,
                       Ih = equilibrium$Ih,
                       Sv = equilibrium$Sv,
                       Ev = equilibrium$Ev,
                       Iv = equilibrium$Iv)


prop_C = 0.12 # 12% of the population are children under 5
theta = 0.5  # aging children who retain some protection (partial protection)
prop_proph = 0.20 # prophylaxis coverage (baseline value)
nu = 1/(5*365) #Aging rate (1/5*365)
prev <- 0.1 ## prevalence of resistant strain
alpha_parm = 5 ## wild type and resistant weibull parms
beta_parm = 30 ## wild type strain duration of protection
beta_res = 11.7 ## resistant strain duration of protection

equilibrium_two <- c(
  # Susceptible population (subscript c is for children and m is for the adult population)
  Shc = equilibrium$Sh  *  (1-prop_proph) * prop_C,
  Shm = equilibrium$Sh  * (1 - prop_C),
  
  # Wild type strain
  Icw = equilibrium$Ih  *  (1 - prev) * (1 - prop_proph) * prop_C,
  Imw = equilibrium$Ih  *  (1-prev)  * (1- prop_C),
  
  # Resistant
  Icr = equilibrium$Ih  *  prev * (1 - prop_proph) * prop_C,
  Imr = equilibrium$Ih  * prev * (1- prop_C),
  
  # Prophylaxis (only for under 5)
  Pc  = equilibrium$Sh  * prop_proph * prop_C,
  Pm = equilibrium$Sh * prop_C *prop_proph * theta,
  Icpw = equilibrium$Ih * (1-prev) * prop_proph * prop_C,
  Icpr = equilibrium$Ih * prev * prop_proph * prop_C,
  Impw = equilibrium$Ih * prop_C * (1-prev) * prop_proph * theta,
  Impr = equilibrium$Ih * prop_C * prev * prop_proph * theta,
  
  # Mosquito/ vector population
  Sv  = equilibrium$Ev, 
  Evw = equilibrium$Ev * (1 - prev),
  Evr = equilibrium$Ev * prev,
  Ivw = equilibrium$Iv * (1 - prev),                     
  Ivr = equilibrium$Iv * prev)

## additional efficacy functions
weibull <- function(x, beta = beta_parm, alpha = alpha_parm) {
  pow <- -(x/beta)^alpha
  y <- exp(pow)
  return(y)
}

ipt_ptp <- function(time_days) {
  dose_vec <- c(seq(90, by = 365, length = 30), #12th week
                seq(114, by = 365, length = 30), # 20-24th week
                seq(194, by = 365, length = 30)) # 28-32th week
  
  
  dose_vec <- sort(dose_vec)
  ## when was their last dose of IPT
  if(time_days < min(dose_vec)) {
    protection <- 0
  } else {
    recent_dose <- dose_vec[max(which((time_days-dose_vec)>0))]
    ## days since IPT dose
    t <- time_days - recent_dose
    protection <- weibull(x = t, alpha = alpha_parm, beta = beta_parm)
  }
  return(protection)
}

ipt_ptp_res <- function(time_days) {
  dose_vec <- c(seq(90, by = 365, length = 30), #12th week for pregnant women
                seq(114, by = 365, length = 30), # 20-24th week
                seq(194, by = 365, length = 30)) # 28-32th week
  
  dose_vec <- sort(dose_vec)
  ## when was their last dose of IPT
  if(time_days < min(dose_vec)) {
    protection <- 0
  } else {
    recent_dose <- dose_vec[max(which((time_days-dose_vec)>0))]
    ## days since IPT dose
    t <- time_days - recent_dose
    protection <- weibull(x = t, alpha = alpha_parm, beta = beta_res)
  }
  return(protection)
}


rm_two_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    Nh <- Shc + Shm + Icw + Imw + Icr + Imr + Pc + Pm + Icpw + Icpr + Impw + Impr # total human population
    Nv <- Sv + Evw + Evr + Ivw + Ivr # total vector population
    
    ## force of infections are strain dependent
    lambda_hw <- ((a * b_h)/Nh) * Ivw
    lambda_hr <- ((a * b_h)/Nh) * Ivr
    lambda_vw <- ((a * b_v)/Nh) * (Icw + Icpw+ Imw + Impw)
    lambda_vr <- ((a * b_v)/Nh) * (Icr + Icpr+ Imr + Impr)
    
    ## vector birth rate
    B <- (0.2 * cos((time * 2 * pi/365) + 3*pi/2) + 1) * mu_v
    ## protection from SMC based on the weibull function
    Ww <- ipt_ptp(time_days = time)
    Wr <- ipt_ptp_res(time_days = time) ## strain dependent protection
    
    # The differential equations
    # Host population dynamics:
    # # those who don't receive prophylaxis
    
    # Without prophylaxis (Children)
    dShc <- - (lambda_hw * Shc) - (lambda_hr * Shc) + (gamma * (Icw + Icr)) - (nu*Shc)
    dIcw <- (lambda_hw * Shc) - (gamma * Icw) - (nu*Icw)
    dIcr <- (lambda_hr * Shc) - (gamma * Icr) - (nu*Icr)
    
    # Without prophylaxis (Adults)
    dShm <- - (lambda_hw * Shm) - (lambda_hr * Shm) + (gamma * (Imw + Imr)) + (nu*Shc)
    dImw <- (lambda_hw * Shm) - (gamma * Imw) + (nu*Icw)
    dImr <- (lambda_hr * Shm) - (gamma * Imr) + (nu*Icr)
    
    ## prophylaxis (Children)
    dPc <- - ((1 - Ww) * lambda_hw * Pc) - ((1 - Wr) * lambda_hr * Pc) + (gamma * (Icpw + Icpr)) - (nu*Pc)
    dIcpw <- (1 - Ww) * (lambda_hw * Pc) - (gamma * Icpw) - (nu*Icpw)
    dIcpr <- (1 - Wr) * (lambda_hr * Pc) - (gamma * Icpr) - (nu*Icpr)
    
    ## prophylaxis (Adults) aged children retain some protection
    dPm <- - ((1 - Ww) * lambda_hw * Pm) - ((1 - Wr) * lambda_hr * Pm) + (gamma * (Impw + Impr)) + (nu* theta * Pc)
    dImpw <- (1 - Ww) * (lambda_hw * Pm) - (gamma * Impw) + (nu * theta * Icpw)
    dImpr <- (1 - Wr) * (lambda_hr * Pm) - (gamma * Impr) + (nu * theta* Icpr)
    
    # Vector population dynamics:
    dSv <- (B * Nv) - ((lambda_vw + lambda_vr) * Sv) - (mu_v * Sv)
    dEvw <- (lambda_vw * Sv) - (epsilon * Evw) - (mu_v * Evw)
    dEvr <- (lambda_vr * Sv) - (epsilon * Evr) - (mu_v * Evr)
    dIvw <- (epsilon * Evw) - (mu_v * Ivw)
    dIvr <- (epsilon * Evr) - (mu_v * Ivr)
    
    
    # Output
    return(list(c(dShc, dIcw, dIcr, # no prophylaxis(Children)
                  dShm, dImw, dImr, # no prophylaxis(Adults)
                  dPc, dIcpw, dIcpr, # prophylaxis (Children)
                  dPm, dImpw, dImpr, # prophylaxis (Adults)
                  dSv, dEvw, dEvr, dIvw, dIvr # vectors 
    ))) 
  })
}

times <- seq(from = 0, to = 365*10, by = 1)    # run for long enough to see the impact of IPT

parameters_two <- c(a = 0.3,       # mosquito biting rate per day
                    b_v = 0.2,   # probability of infection from an infected host to a susceptible vector
                    b_h = 0.40,   # probability of infection from an infected vector to a susceptible host
                    theta = 0.5,  # aging children who retain some protection (partial protection)
                    mu_v = 1/14,  # mortality/recruitment rate of the vector
                    gamma = 1/8,   # recovery rate from symptomatic treated individuals ~(5-14) days
                    nu = 1/(5*365), #Aging rate for the children
                    epsilon = 1/10, # rate of transition from exposed to infected vectors ~ 10 days
                    prop_proph = 0.20, # proportion of the population receiving prophylaxis
                    prev = 0.1, # Resistance prevalence at baseline
                    alpha_parm = 5, # weibull parameter
                    beta_parm = 30, # beta wild type
                    beta_res = 11.7) # beta resistant strain

out_two <- as.data.frame(ode(y = equilibrium_two, 
                             times = times, 
                             func = rm_two_model,
                             parms = parameters_two))
head(out_two)

out_two_long <- melt(as.data.frame(out_two), id = "time")  # turn output dataset into long format


g1<- ggplot(data = filter(out_two_long, variable %in% c( "Icw","Imw","Icpw", "Impw","Icr", "Imr","Icpr", "Impr")),                                               
            aes(x = time, y = value, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (days)")+                                                   
  ylab("Number of people") +                                     
  labs(colour = "Compartment",                                           
       title = "Malaria model + Seasonality + IPT 20%") 

g1
