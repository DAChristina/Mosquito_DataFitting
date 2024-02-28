# bootstrapped confidence interval
# SOURCE: https://epurdom.github.io/Stat131A/book/curve-fitting.html

# 1. Data used mfE vs. Mosquito proportion with Larvae #########################
# Data needed:
# dat_sel$Human_mf_intensity_per20uL
# dat_sel$Mosquito_larvae_infected_proportion_fromtotaldissected
# dat_sel$mf_mean_arithmetic_pertotaldissectedmosquito
# dat_sel$mfE_mean_calculated
# dat_sel$Larvae_mean_calculated

# 1. Data preparation ##########################################################
wd_R_mosquitoes_Reference = "/home/ron/Downloads/Mosquito_DataFitting-main"
setwd(wd_R_mosquitoes_Reference)

library(tidyverse)
library(readxl)

dat <- read_excel('Mosquito_EstablishedInfection_Data.xlsx') %>% 
  # view() %>% 
  glimpse()

dat_sel <- dat %>% 
  select(Reference, Human_mf_intensity_per1mL, Human_mf_intensity_per20uL, Mosquito_species, Mosquito_totaldissected, Mosquito_mf_infected_count, Mosquito_mf_infected_proportion, mf_mean_arithmetic_pertotaldissectedmosquito, Mosquito_larvae_infected_count, Mosquito_larvae_infected_proportion_fromtotaldissected, Larvae_mean_arithmetic_perdissectedmosquito) %>% 
  filter(grepl("^[0-9.]+$", Human_mf_intensity_per20uL)) %>% # Filter out info (e.g. "Captured mosquitoes," etc.)
  filter(grepl("^[0-9.]+$", Mosquito_totaldissected)) %>% # So used data based on mfH & total dissected of mosquitoes
  filter(Mosquito_species != "An. melas") %>% # Omit An. melas
  mutate(Human_mf_intensity_per1mL = as.numeric(Human_mf_intensity_per1mL)) %>% 
  mutate(Human_mf_intensity_per20uL = as.numeric(Human_mf_intensity_per20uL)) %>% 
  mutate(Mosquito_totaldissected = as.numeric(Mosquito_totaldissected)) %>% 
  mutate(Mosquito_mf_infected_count = as.numeric(Mosquito_mf_infected_count)) %>% 
  mutate(Mosquito_mf_infected_proportion = as.numeric(Mosquito_mf_infected_proportion)) %>% 
  mutate(mf_mean_arithmetic_pertotaldissectedmosquito = as.numeric(mf_mean_arithmetic_pertotaldissectedmosquito)) %>% 
  mutate(Mosquito_larvae_infected_count = as.numeric(Mosquito_larvae_infected_count)) %>% 
  mutate(Mosquito_larvae_infected_proportion_fromtotaldissected = as.numeric(Mosquito_larvae_infected_proportion_fromtotaldissected)) %>% 
  mutate(Larvae_mean_arithmetic_perdissectedmosquito = as.numeric(Larvae_mean_arithmetic_perdissectedmosquito)) %>% 
  
  # Calculations based on Snow et al (2002 & 2006)
  mutate(Human_mf20uL_log_calculated = .1149*(log10(Human_mf_intensity_per1mL+1))^2 + .037*(log10(Human_mf_intensity_per1mL+1)) - .0309) %>% 
  mutate(Human_mf20uL_calculated = 10^Human_mf20uL_log_calculated) %>% 
  mutate(mfE_log_mean_calculated = log10(-.0061+1)+.4010*log10(Human_mf_intensity_per20uL+1)) %>% # log10(mfE(mfH)) = log10(-0.0061+1)+0.4010*log10(mfH+1)
  mutate(mfE_mean_calculated = 10^mfE_log_mean_calculated) %>%
  mutate(Larvae_mean_calculated = (.03*mfE_mean_calculated^2.3520)/(1+((mfE_mean_calculated^2.3520)/482.5680))) %>% 
  # mutate(MosquitoInfected_calculated = ) # Analyse this part by using MLE
  glimpse()

# 2. SUCCESSSS glm #############################################################
# TRIAL use glmm by family.glmm = binomial # FAILED coz' package inconsistencies
# 2.1. Use mfE vs. Proportion, ALL SPECIES #####################################
dat_anal <- dat_sel %>% 
  select(Reference, Mosquito_species, Mosquito_larvae_infected_proportion_fromtotaldissected, Mosquito_totaldissected, Mosquito_larvae_infected_count, mfE_mean_calculated) %>% 
  filter(!is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed") %>%
  mutate(Unique_ID = row_number()) %>% 
  # view() %>% 
  glimpse()

# To simplify my life a bit
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

# THE BOOTSTRAP! ###############################################################
bootstrapGLM <- function(y, x, repetitions, confidence.level = 0.95) {
  fit <- glm(y ~ x, family = binomial(link = "logit"),
             weights = AllMosq)
  stat.obs <- coef(fit)
  
  # Bootstrap function
  bootFun <- function() {
    sampled <- sample(1:length(y), size = length(y), replace = TRUE) # Replacement
    fit_sampled <- glm(y[sampled] ~ x[sampled], family = binomial(link = "logit"),
                       weights = AllMosq)
    
    coef(fit_sampled)
  }
  
  stat.boot <- replicate(repetitions, bootFun())
  
  # Format result
  nm <- deparse(substitute(x))
  row.names(stat.boot)[2] <- nm
  level <- 1 - confidence.level
  confidence.interval <- apply(stat.boot, 1, quantile, probs = c(level/2, 1 - level/2))
  
  # CI list and bootstrapped coefficients
  return(list(
    confidence.interval = cbind(lower = confidence.interval[1,], 
                                estimate = stat.obs, 
                                upper = confidence.interval[2,]),
    bootStats = stat.boot
  ))
}

# Call bootstrap with ORI data
result <- bootstrapGLM(y = dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, 
                       x = dat_anal$mfE_mean_calculated, 
                       repetitions = 1000, 
                       confidence.level = 0.95)
result

# PLOT #########################################################################
# Confidence intervals
private_conf <- result$confidence.interval
private_conf

# Plot the data
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     xlim = c(0,10), ylim = c(0,1), col = "grey50", cex = AllMosq/80)

lower <- private_conf[,1]
estim <- private_conf[,2] # the equation!
upper <- private_conf[,3]

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  1/(1+exp(-(l1+l2*x)))
}

curve(logit(x, estim[1], estim[2]), col = "grey10", add = TRUE, lwd = 2) # the equation based on bootstrap
curve(logit(x, upper[1], upper[2]), col = "grey35", add = TRUE, lty = 2) # Upper bound
curve(logit(x, lower[1], lower[2]), col = "grey35", add = TRUE, lty = 2) # Lower bound


# 2.1. Use mfE vs. Proportion, An. gambiae #####################################
# par(mfrow = c(1,2)) # for plotting together An. g & An. a in one output
dat_anal <- dat_sel %>% 
  select(Mosquito_species, Mosquito_larvae_infected_proportion_fromtotaldissected, Mosquito_totaldissected, Mosquito_larvae_infected_count, mfE_mean_calculated) %>% 
  filter(!is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed",
         Mosquito_species == "An. gambiae") %>%
  mutate(Unique_ID = row_number()) %>% 
  # view() %>% 
  glimpse()

# To simplify my life a bit
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

# THE BOOTSTRAP! ###############################################################
bootstrapGLM <- function(y, x, repetitions, confidence.level = 0.95) {
  fit <- glm(y ~ x, family = binomial(link = "logit"),
             weights = AllMosq)
  stat.obs <- coef(fit)
  
  # Bootstrap function
  bootFun <- function() {
    sampled <- sample(1:length(y), size = length(y), replace = TRUE) # Replacement
    fit_sampled <- glm(y[sampled] ~ x[sampled], family = binomial(link = "logit"),
                       weights = AllMosq)
    
    coef(fit_sampled)
  }
  
  stat.boot <- replicate(repetitions, bootFun())
  
  # Format result
  nm <- deparse(substitute(x))
  row.names(stat.boot)[2] <- nm
  level <- 1 - confidence.level
  confidence.interval <- apply(stat.boot, 1, quantile, probs = c(level/2, 1 - level/2))
  
  # CI list and bootstrapped coefficients
  return(list(
    confidence.interval = cbind(lower = confidence.interval[1,], 
                                estimate = stat.obs, 
                                upper = confidence.interval[2,]),
    bootStats = stat.boot
  ))
}

# Call bootstrap with ORI data
result <- bootstrapGLM(y = dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, 
                       x = dat_anal$mfE_mean_calculated, 
                       repetitions = 1000, 
                       confidence.level = 0.95)
result

# PLOT #########################################################################
# Confidence intervals
private_conf <- result$confidence.interval
private_conf

# Plot the data
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     xlim = c(0,10), ylim = c(0,1), col = "red", cex = AllMosq/80)

lower <- private_conf[,1]
estim <- private_conf[,2] # the equation!
upper <- private_conf[,3]

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  1/(1+exp(-(l1+l2*x)))
}

curve(logit(x, estim[1], estim[2]), col = "darkred", add = TRUE, lwd = 2) # the equation based on bootstrap
curve(logit(x, upper[1], upper[2]), col = "violetred1", add = TRUE, lty = 2) # Upper bound
curve(logit(x, lower[1], lower[2]), col = "violetred1", add = TRUE, lty = 2) # Lower bound


# 2.1. Use mfE vs. Proportion, An. arabiensis #####################################
dat_anal <- dat_sel %>% 
  select(Mosquito_species, Mosquito_larvae_infected_proportion_fromtotaldissected, Mosquito_totaldissected, Mosquito_larvae_infected_count, mfE_mean_calculated) %>% 
  filter(!is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed",
         Mosquito_species == "An. arabiensis") %>%
  mutate(Unique_ID = row_number()) %>% 
  # view() %>% 
  glimpse()

# To simplify my life a bit
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

# THE BOOTSTRAP! ###############################################################
bootstrapGLM <- function(y, x, repetitions, confidence.level = 0.95) {
  fit <- glm(y ~ x, family = binomial(link = "logit"),
             weights = AllMosq)
  stat.obs <- coef(fit)
  
  # Bootstrap function
  bootFun <- function() {
    sampled <- sample(1:length(y), size = length(y), replace = TRUE) # Replacement
    fit_sampled <- glm(y[sampled] ~ x[sampled], family = binomial(link = "logit"),
                       weights = AllMosq)
    
    coef(fit_sampled)
  }
  
  stat.boot <- replicate(repetitions, bootFun())
  
  # Format result
  nm <- deparse(substitute(x))
  row.names(stat.boot)[2] <- nm
  level <- 1 - confidence.level
  confidence.interval <- apply(stat.boot, 1, quantile, probs = c(level/2, 1 - level/2))
  
  # CI list and bootstrapped coefficients
  return(list(
    confidence.interval = cbind(lower = confidence.interval[1,], 
                                estimate = stat.obs, 
                                upper = confidence.interval[2,]),
    bootStats = stat.boot
  ))
}

# Call bootstrap with ORI data
result <- bootstrapGLM(y = dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, 
                       x = dat_anal$mfE_mean_calculated, 
                       repetitions = 1000, 
                       confidence.level = 0.95)
result

# PLOT #########################################################################
# Confidence intervals
private_conf <- result$confidence.interval
private_conf

# Plot the data
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     xlim = c(0,10), ylim = c(0,1), col = "blue", cex = AllMosq/80)

lower <- private_conf[,1]
estim <- private_conf[,2] # the equation!
upper <- private_conf[,3]

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  1/(1+exp(-(l1+l2*x)))
}

curve(logit(x, estim[1], estim[2]), col = "darkblue", add = TRUE, lwd = 2) # the equation based on bootstrap
curve(logit(x, upper[1], upper[2]), col = "steelblue", add = TRUE, lty = 2) # Upper bound
curve(logit(x, lower[1], lower[2]), col = "steelblue", add = TRUE, lty = 2) # Lower bound

par(mfrow = c(1,1))
