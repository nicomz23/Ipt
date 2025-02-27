## Two strain Ross-Macdonald model 


# PACKAGES
require(deSolve)
require(reshape2)
require(ggplot2)
require(tidyverse)
require(gridExtra)
require(latticeExtra)

## ggplot theme
theme <- theme_set(theme_bw())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))

Nh <- 10000   # total number of hosts
Nv <- 150000   # total number of vectors

## one strain model to get to an equilibrium, and then manually set starting prevalence
## initial state values for the one strain model to run to equilibrium
## will run the model for a while without chemoprevention to see where it stabilises to use as a starting point
initial_state_values <- c(Sh = 0.9 * Nh,
                          Eh = 0.1 * Nh,
                          Is = 0,   #Symptomatic treated
                          Ia = 0,    #Asymptomatic untreated
                          Sv = Nv,
                          Ev = 0,
                          Iv = 0)

## parameters
parameters <- c(a = 1/3,      # mosquito biting rate per day
                b_v = 0.2,    # probability of infection from an infected human to a susceptible vector
                b_h = 0.40,    # probability of infection from an infected vector to a susceptible human
                mu_v = 1/14,   # mortality/recruitment rate of the vector
                alpha = 1/13, # period between infection to infectiousness  ~ 12-14 days
                d = 0.3,     # proportion of infected human who develop symptoms
                gamma_a = 1/120,   # recovery rate from asymptomatic untreated humans ~ 100-200days
                gamma_s = 1/5,  # recovery from symptomatic treated humans ~ 5-7 days
                epsilon = 1/10) # rate of transition from exposed to infected vectors ~ 10 days

times <- seq(from = 0, to = 365, by = 1)    # simulate for a length of time to find equilibrium

rm_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    Nh <- Sh + Eh + Ia + Is     # total human population
    Nv <- Sv + Ev + Iv # total vector population
    
    ## force of infections
    lambda_h <- ((a * b_h)/Nh) * Iv
    lambda_v <- ((a * b_v)/Nh) * (Is + Ia)
    
    # The differential equations
    # Host population dynamics:
    dSh <- - (lambda_h * Sh) + (gamma_s * Is) + (gamma_a * Ia)
    dEh <-  (lambda_h * Sh) -  (alpha * Eh)
    dIs <-   ((d) * alpha * Eh)  - (gamma_s * Is)
    dIa <-   ((1 - d) * alpha * Eh)  - (gamma_a * Ia)
    
    
    
    # Vector population dynamics:
    dSv <- (mu_v * Nv) - (lambda_v * Sv) - (mu_v * Sv)         
    dEv <- (lambda_v * Sv) - (epsilon * Ev) - (mu_v * Ev)
    dIv <- (epsilon * Ev) - (mu_v * Iv)
    
    # Output
    return(list(c(dSh, dEh, dIs, dIa, dSv, dEv, dIv))) 
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

equilibrium <- output_one %>% dplyr::filter(time == 210)

equilibrium_start <- c(Sh = equilibrium$Sh,
                       Eh = equilibrium$Eh,
                       Is = equilibrium$Is,
                       Ia = equilibrium$Ia,
                       Sv = equilibrium$Sv,
                       Ev = equilibrium$Ev,
                       Iv = equilibrium$Iv)

prop_proph = 0.32 # prophylaxis coverage (baseline value)
prev <- 0.1 ## prevalence of resistant strain
alpha_parm = 5 ## wild type and resistant weibull parms
beta_parm = 30 ## wild type strain duration of protection
beta_res = 11.7 ## resistant strain duration of protection

equilibrium_two <- c(Sh = equilibrium$Sh  *  (1-prop_proph),
                     Ew = equilibrium$Eh  *  (1 - prev) * (1 - prop_proph),
                     Isw = equilibrium$Is  *  (1 - prev) * (1 - prop_proph),
                     Iaw = equilibrium$Ia  *  (1- prev) * (1 - prop_proph),
                     Er = equilibrium$Eh  *  prev * (1 - prop_proph),
                     Isr = equilibrium$Is * prev * (1 - prop_proph) ,
                     Iar = equilibrium$Ia * prev *  (1 - prop_proph),
                     P  = equilibrium$Sh  * prop_proph,
                     Epw = equilibrium$Eh * (1 - prev) * prop_proph,
                     Ispw = equilibrium$Is * (1 - prev) * prop_proph,
                     Iapw = equilibrium$Ia * (1 - prev) * prop_proph,
                     Epr = equilibrium$Eh *  prev * prop_proph,
                     Ispr = equilibrium$Is * prev * prop_proph,
                     Iapr = equilibrium$Ia *  prev * prop_proph,
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
    
    Nh <- Sh + Ew + Isw + Iaw + Er + Isr + Iar  + P + Epw + Ispw + Iapw + Epr + Ispr + Iapr # total human population
    Nv <- Sv + Evw + Evr + Ivw + Ivr # total vector population
    
    ## force of infections are strain dependent
    lambda_hw <- ((a * b_h)/Nh) * Ivw
    lambda_hr <- ((a * b_h)/Nh) * Ivr
    lambda_vw <- ((a * b_v)/Nh) * (Isw + Iaw+ Ispw + Iapw)
    lambda_vr <- ((a * b_v)/Nh) * (Isr+ Iar+ Ispr + Iapr)
    
    ## vector birth rate
    B <- (0.2 * cos((time * 2 * pi/365) + 3*pi/2) + 1) * mu_v
    ## protection from SMC based on the weibull function
    Ww <- ipt_ptp(time_days = time)
    Wr <- ipt_ptp_res(time_days = time) ## strain dependent protection
    
    # The differential equations
    # Host population dynamics:
  
    # # Model without chemoprevention
    
    dSh <- - (lambda_hw * Sh) - (lambda_hr * Sh) + (gamma_s * (Isw + Isr)) + (gamma_a * (Iaw + Iar))
    dEw <- (lambda_hw * Sh) - (alpha * Ew)
    dIsw <- ((d * alpha) * Ew) - (gamma_s * Isw)
    dIaw <- ((1-d) * alpha * Ew) - (gamma_a * Iaw)
    dEr <-  (lambda_hr * Sh) - (alpha * Er)
    dIsr <- ((d * alpha) * Er) - (gamma_s * Isr)
    dIar <- ((1-d) * alpha * Er) - (gamma_a * Iar)
    
    ## Model with chemoprevention
    
    dP <- - ((1 - Ww) * lambda_hw * P) - ((1 - Wr) * lambda_hr * P) + (gamma_s * (Ispw + Ispr)) + (gamma_a * (Iapw + Iapr))
    dEpw <-  (1 - Ww) * (lambda_hw * P)  - (alpha * Epw)
    dIspw <- ((d * alpha) * Epw) - (gamma_s * Ispw)
    dIapw <- ((1-d) * alpha * Epw) - (gamma_a * Iapw)
    dEpr <-  (1 - Wr) * (lambda_hr * P)  - (alpha * Epr)
    dIspr <- ((d * alpha) * Epr) - (gamma_s * Ispr)
    dIapr <- ((1-d) * alpha * Epr) - (gamma_a * Iapr)
    
    # Vector population dynamics:
    
    dSv <- (B * Nv) - ((lambda_vw + lambda_vr) * Sv) - (mu_v * Sv)
    dEvw <- (lambda_vw * Sv) - (epsilon * Evw) - (mu_v * Evw)
    dEvr <- (lambda_vr * Sv) - (epsilon * Evr) - (mu_v * Evr)
    dIvw <- (epsilon * Evw) - (mu_v * Ivw)
    dIvr <- (epsilon * Evr) - (mu_v * Ivr)
    
    
    # Output
    return(list(c(dSh, dEw, dIsw, dIaw, dEr, dIsr, dIar, # no prophylaxis
                  dP, dEpw, dIspw, dIapw, dEpr, dIspr, dIapr, # prophylaxis
                  dSv, dEvw, dEvr, dIvw, dIvr # vectors 
    ))) 
  })
}


times <- seq(from = 0, to = 365*10, by = 1)    # run for long enough to see the impact of chemoprevention

parameters_two <- c(a = 1/3,       # mosquito biting rate per day
                    b_v = 0.2,   # probability of infection from an infected host to a susceptible vector
                    b_h = 0.40,   # probability of infection from an infected vector to a susceptible host
                    mu_v = 1/14,  # mortality/recruitment rate of the vector
                    gamma_a = 1/120,   # recovery rate from asymptomatic treated humans ~ 100-200days
                    gamma_s = 1/5, #  recovery from symptomatic untreated humans ~ 12-14 days  
                    epsilon = 1/10, # rate of transition from exposed to infected vectors ~ 10 days
                    prop_proph = 0.32, # proportion of the population receiving prophylaxis
                    prev = 0.1, # Resistance prevalence at baseline
                    d  = 0.3 , # Proportion of individual who develop symptoms
                    alpha   = 1/13, # Period between infection to infectiousness  ~ 12-14 days
                    alpha_parm = 5, # weibull parameter
                    beta_parm = 30, # beta wild type
                    beta_res = 11.7) # beta resistant strain

out_two <- as.data.frame(ode(y = equilibrium_two, 
                             times = times, 
                             func = rm_two_model,
                             parms = parameters_two))
head(out_two)

out_two_long <- melt(as.data.frame(out_two), id = "time")  # turn output dataset into long format


g1<- ggplot(data = filter(out_two_long, variable %in% c( "Isw","Iaw","Ispw", "Iapw","Isr", "Iar","Ispr", "Iapr")),                                               
            aes(x = time, y = value, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (days)")+                                                   
  ylab("Number of people") +                                     
  labs(colour = "Compartment",                                           
       title = "Malaria model + Seasonality + IPT 32%") 

g1


output_two_prop <- out_two
output_two_prop$Wild_type <- ( out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw)/(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)
output_two_prop$Resistance <- (out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)/(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)

output_two_prop_long <- melt(as.data.frame(output_two_prop), id = "time")  # turn output dataset into long format


g2 <- ggplot(data = filter(output_two_prop_long, variable %in% c("Wild_type", "Resistance")),                                               
             aes(x = time/365, y = value, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("Proportion") +                                     
  labs(colour = "Compartment",                                           
       title = "Malaria model + Seasonality + IPT 32%")
g2


# plot the protection over time for each strain
protection <- data.frame(days = 1:365, wild_type = numeric(365), resistant = numeric(365))
for (i in 1:nrow(protection)) {
  protection$wild_type[i] <- ipt_ptp(time_days = protection$days[i])
  protection$resistant[i] <- ipt_ptp_res(time_days = protection$days[i])
}
protection <- protection |>
  tidyr::pivot_longer(wild_type:resistant, names_to = "strain", values_to = "protective_efficacy")

ggplot(protection, aes(x = days, y = protective_efficacy, col = strain)) + geom_line()

R=sum(out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)/sum(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)
R

S=sum(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw)/sum(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)
S

g3 <- ggplot(data = filter(output_two_prop_long, variable %in% c("Resistance")),                                               
             aes(x = time/365, y = value*100, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infection") +                                     
  labs(colour = "Compartment",                                           
       title = "Malaria model + Seasonality + IPT 32%")

g3


g4 <- ggplot(data = filter(output_two_prop_long, variable %in% c("Wild_type")),                                               
             aes(x = time/365, y = value*100, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infections") +                                     
  labs(colour = "Compartment",                                           
       title = "Malaria model + Seasonality + IPT 32%")

g4
################################################################################
# Run simulations for different prop_proph values (different coverages)
################################################################################
prop_proph_values <-c(0,0.32,0.45,0.6)
results <- lapply(prop_proph_values, function(prop_proph) {
  parameters <- c(parameters_two, prop_proph = prop_proph)
  equilibrium_two <- c(Sh = equilibrium$Sh  *  (1-prop_proph),
                       Ew = equilibrium$Eh  *  (1 - prev) * (1 - prop_proph),
                       Isw = equilibrium$Is  *  (1 - prev) * (1 - prop_proph),
                       Iaw = equilibrium$Ia  *  (1- prev) * (1 - prop_proph),
                       Er = equilibrium$Eh  *  prev * (1 - prop_proph),
                       Isr = equilibrium$Is * prev * (1 - prop_proph) ,
                       Iar = equilibrium$Ia * prev *  (1 - prop_proph),
                       P  = equilibrium$Sh  * prop_proph,
                       Epw = equilibrium$Eh * (1 - prev) * prop_proph,
                       Ispw = equilibrium$Is * (1 - prev) * prop_proph,
                       Iapw = equilibrium$Ia * (1 - prev) * prop_proph,
                       Epr = equilibrium$Eh *  prev * prop_proph,
                       Ispr = equilibrium$Is * prev * prop_proph,
                       Iapr = equilibrium$Ia *  prev * prop_proph,
                       Sv  = equilibrium$Ev, 
                       Evw = equilibrium$Ev * (1 - prev),
                       Evr = equilibrium$Ev * prev,
                       Ivw = equilibrium$Iv * (1 - prev),                     
                       Ivr = equilibrium$Iv * prev)
  
  output_spc <- as.data.frame(ode(y = equilibrium_two, times = times, func = rm_two_model, parms = parameters))
  output_spc$prop_proph <- paste0("Prophylaxis: ", round(prop_proph * 100, 1), "%")
  return(output_spc)
  
})
df <- bind_rows(results) #Combine all results
df_long <-melt(df,id = c("time", "prop_proph")) # convert data to long format

g5 <- ggplot(data = filter(df_long, variable %in% c("Isw","Iaw","Ispw", "Iapw","Isr", "Iar","Ispr", "Iapr")), 
             aes(x = time, y = value, colour = prop_proph, group = variable)) +  
  geom_line(size=1) +                                                           
  xlab("Time (years)") +                                                   
  ylab("Number of infected hosts") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on Malaria Infections")
g5


output_two1_prop <- df
output_two1_prop$Wild_type <- ( df$Isw + df$Iaw + df$Ispw + df$Iapw)/(df$Isw + df$Iaw + df$Ispw + df$Iapw + df$Isr + df$Iar+df$Ispr+df$Iapr)
output_two1_prop$Resistance <- (df$Isr + df$Iar+df$Ispr+df$Iapr)/(df$Isw + df$Iaw + df$Ispw + df$Iapw + df$Isr + df$Iar+df$Ispr+df$Iapr)

output_two1_prop_long <- melt(as.data.frame(output_two1_prop), id = c("time","prop_proph")) # turn output dataset into long format


g6 <- ggplot(data = filter(output_two1_prop_long, variable %in% c("Resistance")),                                               
             aes(x = time/365, y = value*100, colour = prop_proph, group = prop_proph)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infection") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on the Resistant Strain")

g6

g7 <- ggplot(data = filter(output_two1_prop_long, variable %in% c("Wild_type")),                                               
             aes(x = time/365, y = value*100, colour = prop_proph, group = prop_proph)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infection") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on the Wild type strain")

g7




