# 1. Data used mfE vs. Mosquito proportion with Larvae #########################
# Data needed:
# dat_sel$Human_mf_intensity_per20uL
# dat_sel$Mosquito_larvae_infected_proportion_fromtotaldissected
# dat_sel$mf_mean_arithmetic_pertotaldissectedmosquito
# dat_sel$mfE_mean_calculated
# dat_sel$Larvae_mean_calculated

# 1. Data preparation ##########################################################
wd_R_mosquitoes_Reference = "C:/Users/dac23/Downloads/Mosquito_DataFitting-main/Mosquito_DataFitting-main"
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

# Load the lme4 package
library(lme4)

# Fit GLMM
fit <- glmer(cbind(Mosquito_larvae_infected_count, Mosquito_totaldissected - Mosquito_larvae_infected_count) ~ mfE_mean_calculated + (1 | Unique_ID), # TRIAL fit by Reference & Mosquito_species
               family = binomial,
               data = dat_anal)

# For model comparisons:
fit_ALL_glmer <- fit

# View model summary
summary(fit)

# Model Result analysis visualisation??? #######################################
par(mfrow = c(2,2))
# Residual???
residuals <- resid(fit, type = "pearson") # trial deviance?
fitted_values <- fitted(fit)
plot(fitted_values, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted Values Plot")
abline(h = 0)

# qq???
qqnorm(residuals(fit))
qqline(residuals(fit))

# Cook's Distance?
plot(cooks.distance(fit), pch = 20, main = "Cook's Distance Plot",
     xlab = "Observation Number", ylab = "Cook's Distance")

par(mfrow = c(1,1))

# Residuals vs. leverage?
# recall residuals above: residuals <- resid(fit, type = "pearson")
# Warning message occurs:
# In hatvalues.merMod(fit) : the hat matrix may not make sense for GLMMs
# leverage <- hatvalues(fit)
# plot(leverage, residuals,
# xlab = "Leverage",
# ylab = "Residuals",
# main = "Residuals vs. Leverage Plot")
par(mfrow = c(1,1))


# 2.1.1. glm visualisation #####################################################
mfE_seq <- seq(min(x), max(x), length.out = 100)

# Predicted probabilities (fitted GLMM)
predicted_probs <- predict(fit, newdata = data.frame(mfE_mean_calculated = mfE_seq), type = "response", re.form = NA)

# Predicted probabilities to log-odds (logit)
log_odds <- log(predicted_probs / (1 - predicted_probs))

l1 <- coef(summary(fit))["(Intercept)", "Estimate"] # also equals to beta_0
l2 <- coef(summary(fit))["mfE_mean_calculated", "Estimate"] # also equals to beta_1

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  exp(l1+l2*x)/(1+exp(l1+l2*x))
}

# 2.1.2. Plot the data points and the fitted curves ############################
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     xlim = c(0,10), ylim = c(0,1), col = "grey50", cex = AllMosq/80)
# lines(mfE_seq, predicted_probs, col = "darkgreen", lwd = 2)  # TRIAL predicted probabilities
# lines(mfE_seq, plogis(log_odds), col = "darkred", lwd = 2)  # TRIAL fitted logit curves
curve(logit(x, l1, l2), col = "grey10", add = TRUE, lwd = 2) # the equation

# 2.1.3. Add CI using Clopper-Pearson method ###################################
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
# pred <- predict(fit, newdata = data.frame(mfE_mean_calculated = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(predicted_probs*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, l1, l2) +1.96*upper_bound -logit(x_new, l1, l2),
      col = "grey35", lty = 5, lwd = 2)
lines(x_new, logit(x_new, l1, l2) -1.96*lower_bound +logit(x_new, l1, l2),
      col = "grey35", lty = 5, lwd = 2)


# 2.1. Use mfE vs. Proportion, An. gambiae #####################################
par(mfrow = c(1,2)) # for plotting together An. g & An. a in one output
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

# Load the lme4 package
library(lme4)

# Fit GLMM
fit <- glmer(cbind(Mosquito_larvae_infected_count, Mosquito_totaldissected - Mosquito_larvae_infected_count) ~ mfE_mean_calculated + (1 | Unique_ID),
             family = binomial,
             data = dat_anal)

# For model comparisons:
fit_Angam_glmer <- fit

# View model summary
summary(fit)
# plot(fit)

# 2.1.1. glm visualisation #####################################################
mfE_seq <- seq(min(x), max(x), length.out = 100)

# Predicted probabilities (fitted GLMM)
predicted_probs <- predict(fit, newdata = data.frame(mfE_mean_calculated = mfE_seq), type = "response", re.form = NA)

# Predicted probabilities to log-odds (logit)
log_odds <- log(predicted_probs / (1 - predicted_probs))

l1 <- coef(summary(fit))["(Intercept)", "Estimate"] # also equals to beta_0
l2 <- coef(summary(fit))["mfE_mean_calculated", "Estimate"] # also equals to beta_1

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  exp(l1+l2*x)/(1+exp(l1+l2*x))
}

# 2.1.2. Plot the data points and the fitted curves ############################
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     main = expression(italic("An. gambiae")),
     xlim = c(0,10), ylim = c(0,1), col = "red", cex = AllMosq/80)
# lines(mfE_seq, predicted_probs, col = "darkgreen", lwd = 2)  # TRIAL predicted probabilities
# lines(mfE_seq, plogis(log_odds), col = "darkred", lwd = 2)  # TRIAL fitted logit curves
curve(logit(x, l1, l2), col = "darkred", add = TRUE, lwd = 2) # the equation

# 2.1.3. Add CI using Clopper-Pearson method ###################################
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
# pred <- predict(fit, newdata = data.frame(mfE_mean_calculated = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(predicted_probs*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, l1, l2) +1.96*upper_bound -logit(x_new, l1, l2),
      col = "violetred1", lty = 5, lwd = 2)
lines(x_new, logit(x_new, l1, l2) -1.96*lower_bound +logit(x_new, l1, l2),
      col = "violetred1", lty = 5, lwd = 2)


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

# Load the lme4 package
library(lme4)

# Fit GLMM
fit <- glmer(cbind(Mosquito_larvae_infected_count, Mosquito_totaldissected - Mosquito_larvae_infected_count) ~ mfE_mean_calculated + (1 | Unique_ID),
             family = binomial,
             data = dat_anal)

# For model comparisons:
fit_Anara_glmer <- fit

# View model summary
summary(fit)
# plot(fit)

# 2.1.1. glm visualisation #####################################################
mfE_seq <- seq(min(x), max(x), length.out = 100)

# Predicted probabilities (fitted GLMM)
predicted_probs <- predict(fit, newdata = data.frame(mfE_mean_calculated = mfE_seq), type = "response", re.form = NA)

# Predicted probabilities to log-odds (logit)
log_odds <- log(predicted_probs / (1 - predicted_probs))

l1 <- coef(summary(fit))["(Intercept)", "Estimate"] # also equals to beta_0
l2 <- coef(summary(fit))["mfE_mean_calculated", "Estimate"] # also equals to beta_1

logit <- function(x,l1,l2) { # Translated from the log-odds function from b0 = l1 & b1 = l2
  exp(l1+l2*x)/(1+exp(l1+l2*x))
}

# 2.1.2. Plot the data points and the fitted curves ############################
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     main = expression(italic("An. arabiensis")),
     xlim = c(0,10), ylim = c(0,1), col = "blue", cex = AllMosq/80)
# lines(mfE_seq, predicted_probs, col = "darkgreen", lwd = 2)  # TRIAL predicted probabilities
# lines(mfE_seq, plogis(log_odds), col = "darkred", lwd = 2)  # TRIAL fitted logit curves
curve(logit(x, l1, l2), col = "darkblue", add = TRUE, lwd = 2) # the equation

# 2.1.3. Add CI using Clopper-Pearson method ###################################
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
# pred <- predict(fit, newdata = data.frame(mfE_mean_calculated = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(predicted_probs*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, l1, l2) +1.96*upper_bound -logit(x_new, l1, l2),
      col = "steelblue", lty = 5, lwd = 2)
lines(x_new, logit(x_new, l1, l2) -1.96*lower_bound +logit(x_new, l1, l2),
      col = "steelblue", lty = 5, lwd = 2)

par(mfrow = c(1,1))

# Chi-squared comparisons of species model?
