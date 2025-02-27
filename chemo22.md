---
title: "Chemoprevention Model"
author: "Nicholaus Mziray"
date: "2025-02-27"
output: 
  word_document: 
    keep_md: yes
---

``` r
#This model investigates the spread of malaria resistant strain in the presence and absence of chemoprevention.
```

# The model is implemented in R, so we need to load all the necessary packages.
# Load the required packages


``` r
#ggplot themes
theme <- theme_set(theme_bw())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
```

# Total population in human and vectors/ mosquitoes is defined here. 

``` r
Nh <- 10000   # total number of hosts
Nv <- 150000   # total number of vectors
```
## We start with a one strain model to get to an equilibrium, and then manually set starting prevalence.The initial state values for the one strain model to run to equilibrium.
## Then we run the model for a while without chemoprevention to see where it stabilises to use as a starting point. For the begining, we have a closed population model, the aim for now is to understand how resistant strain spreads with different levels of chemoprevention coverages.

``` r
initial_state_values <- c(Sh = 0.9 * Nh,
                          Eh = 0.1 * Nh,
                          Is = 0,   #Symptomatic treated
                          Ia = 0,    #Asymptomatic untreated
                          Sv = Nv,
                          Ev = 0,
                          Iv = 0)
```
# Parameters and time to the equilibrium point. We run a one strain model to a duration of one year, to observe where it comes to stablise.

``` r
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
```
# Define a single strain malaria model 

``` r
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
```
# Solve the equestions as in the model

``` r
output_one <- as.data.frame(ode(y = initial_state_values, 
                                times = times, 
                                func = rm_model,
                                parms = parameters))
```
# Solving the system to identify the equilibrium point, the model stablizes at around 200 days, and we take 210 as the starting point for introducing resistance.

``` r
output_one <- as.data.frame(ode(y = initial_state_values, 
                                times = times, 
                                func = rm_model,
                                parms = parameters))

equilibrium <- output_one %>% dplyr::filter(time == 210)
```
# Define the Equilibrium points

``` r
equilibrium_start <- c(Sh = equilibrium$Sh,
                       Eh = equilibrium$Eh,
                       Is = equilibrium$Is,
                       Ia = equilibrium$Ia,
                       Sv = equilibrium$Sv,
                       Ev = equilibrium$Ev,
                       Iv = equilibrium$Iv)
```
# Set the parameters for the two strain malaria model. For the chemoprevention, we have different durations of protection against wild type and resistant strains.

``` r
prop_proph = 0.32 # prophylaxis coverage (baseline value)
prev <- 0.1 ## prevalence of resistant strain
alpha_parm = 5 ## wild type and resistant weibull parms
beta_parm = 30 ## wild type strain duration of protection
beta_res = 11.7 ## resistant strain duration of protection
```
# Define the new equilibrium point where we are going to introduce resistant into the model.

``` r
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
```
# We assume that chemoprevention protection against each strain follows a weibull distribution function and it depends on time.

``` r
weibull <- function(x, beta = beta_parm, alpha = alpha_parm) {
  pow <- -(x/beta)^alpha
  y <- exp(pow)
  return(y)
}
```
# Define the dosage timing for the chemoprevention, here we use three rounds per year spaced in every 30 days. Protection against each strain depends on the days since the last dose.
# Here we define the protection for each strain

``` r
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
```

# Define the full model with and without chemoprevention, here human population is divided into 14 compartments.
# Define also the time to which the simualtion runs (we will run this model for ten years) to see how resistant spreads

``` r
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
```

# Define the parameters for the full model here

``` r
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
```
# Output the model results and convert them into a long format

``` r
out_two <- as.data.frame(ode(y = equilibrium_two, 
                             times = times, 
                             func = rm_two_model,
                             parms = parameters_two))
head(out_two)
```

```
##   time       Sh       Ew      Isw      Iaw       Er      Isr      Iar        P
## 1    0 131.1121 792.5845 91.46974 5117.945 88.06494 10.16330 568.6605 61.69980
## 2    1 131.1857 792.4957 91.46590 5117.971 88.05507 10.16288 568.6635 61.73443
## 3    2 131.6076 792.1095 91.45801 5117.986 88.01217 10.16200 568.6651 61.93299
## 4    3 132.5090 791.3339 91.43930 5117.969 87.92598 10.15992 568.6632 62.35719
## 5    4 133.8953 790.1903 91.40368 5117.900 87.79892 10.15596 568.6556 63.00955
## 6    5 135.7163 788.7455 91.34722 5117.763 87.63839 10.14969 568.6403 63.86651
##        Epw     Ispw     Iapw      Epr     Ispr     Iapr       Sv      Evw
## 1 372.9809 43.04458 2408.445 41.44233 4.782731 267.6050 27671.14 24904.02
## 2 372.9391 43.04278 2408.457 41.43768 4.782531 267.6063 30660.66 22352.97
## 3 372.7574 43.03907 2408.464 41.41749 4.782118 267.6071 33313.82 20335.82
## 4 372.3924 43.03026 2408.456 41.37693 4.781140 267.6062 35672.66 18753.71
## 5 371.8543 43.01350 2408.424 41.31714 4.779278 267.6026 37774.19 17525.20
## 6 371.1744 42.98693 2408.359 41.24159 4.776325 267.5954 39650.98 16583.26
##        Evr      Ivw      Ivr
## 1 2767.114 34862.55 3873.616
## 2 2483.663 34733.44 3859.272
## 3 2259.535 34393.99 3821.554
## 4 2083.746 33905.21 3767.246
## 5 1947.245 33315.30 3701.700
## 6 1842.585 32661.96 3629.107
```

``` r
out_two_long <- melt(as.data.frame(out_two), id = "time")  # turn output dataset into long format
```
# Calculate the proportion of infections that will become wild-type or resistant over time

``` r
output_two_prop <- out_two
output_two_prop$Wild_type <- ( out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw)/(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)
output_two_prop$Resistance <- (out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)/(out_two$Isw + out_two$Iaw + out_two$Ispw + out_two$Iapw + out_two$Isr + out_two$Iar+out_two$Ispr+out_two$Iapr)

output_two_prop_long <- melt(as.data.frame(output_two_prop), id = "time")  # turn output dataset into long format
```
# Run the simulations using different levels of chemoprevention coverages, here we consider 0%, 32%, 45% and 60% coverage for the chemoprevention

``` r
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
```
# Combine all the results and convert them into a long format.

``` r
df <- bind_rows(results) #Combine all results
df_long <-melt(df,id = c("time", "prop_proph")) # convert data to long format
```
# Again calculate the proportion of infections for each strain here and covert them into a long formart.

``` r
output_two1_prop <- df
output_two1_prop$Wild_type <- ( df$Isw + df$Iaw + df$Ispw + df$Iapw)/(df$Isw + df$Iaw + df$Ispw + df$Iapw + df$Isr + df$Iar+df$Ispr+df$Iapr)
output_two1_prop$Resistance <- (df$Isr + df$Iar+df$Ispr+df$Iapr)/(df$Isw + df$Iaw + df$Ispw + df$Iapw + df$Isr + df$Iar+df$Ispr+df$Iapr)

output_two1_prop_long <- melt(as.data.frame(output_two1_prop), id = c("time","prop_proph")) # turn output dataset into long format
```
# Plot for the resistant strain ( We plot to the how fast resistant strain spread) under different coverages for the chemoprevention. 

``` r
ggplot(data = filter(output_two1_prop_long, variable %in% c("Resistance")),                                               
             aes(x = time/365, y = value*100, colour = prop_proph, group = prop_proph)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infection") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on the Resistant Strain")
```

![](chemo22_files/figure-docx/plots-1.png)<!-- -->


#With the 60% coverage of chemoprevention, it takes only 2.5 years to have 50% of all infections resistant, it takes much longer for resistant infections to become 50% of all the infections when chemoprevention coverage is only 32%.


# Plot for the Wild type (sensitive) strain under different coverage levels. 


``` r
ggplot(data = filter(output_two1_prop_long, variable %in% c("Wild_type")),                                               
             aes(x = time/365, y = value*100, colour = prop_proph, group = prop_proph)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("% Infection") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on the Wild type strain")
```

![](chemo22_files/figure-docx/unnamed-chunk-20-1.png)<!-- -->
#Increased coverage levels of chemoprevention reduces the sensitive strain infections. It takes long time for the lower levels of chemoprevention to clear the sensitive infections.


## General plot for the simuated coverage levels (I am exporing the best way to improve this)

``` r
ggplot(data = filter(df_long, variable %in% c("Isw","Iaw","Ispw", "Iapw","Isr", "Iar","Ispr", "Iapr")), 
             aes(x = time, y = value, colour = prop_proph, group = variable)) +  
  geom_line(size=1) +                                                           
  xlab("Time (years)") +                                                   
  ylab("Number of infected hosts") +                                     
  labs(colour = "Prophylaxis Coverage",                                           
       title = "Impact of Chemoprevention on Malaria Infections")
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## ℹ Please use `linewidth` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

![](chemo22_files/figure-docx/unnamed-chunk-21-1.png)<!-- -->
# Further analysis for the model is ongoing, we are going to reduce some assumtions, and see how the model dynamics behaves. Analysis on the total number of the cases averted based on the chemprevention coverage is going to be done. The last sep for this first model is to include/ take into account the age structure, to account for the coverages of chemopreventions into different age groups (Mostly pregnant women and the children under 5).

