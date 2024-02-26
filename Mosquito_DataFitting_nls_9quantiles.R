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

# 2. SUCCESSSS 3 PARAMETERS ####################################################
# 2.2. The logistic regression WITH ITERATION ##################################
# 2.2.1. Use mfE vs. Proportion, 9 breakpoints, ALL SPECIES ####################
dat_anal <- dat_sel %>% 
  select(Mosquito_larvae_infected_proportion_fromtotaldissected, mfE_mean_calculated, Mosquito_totaldissected) %>% 
  filter(!is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed") %>%
  # view() %>% 
  glimpse()

# mfE mean vs. proportion of mosquitoes with larvae
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

plot(dat_anal$mfE_mean_calculated, dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with Larvae",
     color = "grey20", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles
quantiles <- quantile(x, probs=seq(0,1,length=10))
quantiles
mfE_cut <- cut(x, breaks=quantiles, include.lowest=TRUE)
levels(mfE_cut)
p <- tapply(y, mfE_cut, mean)
p

quantiles_Mosq <- quantile(AllMosq, probs=seq(0,1,length=10))
quantiles_Mosq
Dissect_Mosq <- tapply(AllMosq, mfE_cut, mean)
Dissect_Mosq

mfE_cut_midpoints = rep(NA,9) # length 9
for (i in 1:9) {
  mfE_cut_midpoints[i] = (quantiles[i] + quantiles[i+1])/2
}
mfE_cut_midpoints

# Let's fit this: fit <- glm(y ~ x, family=binomial) by iterations
# Set the number of iterations
n_iterations <- 500

# Create empty vectors & matrices to store coefficients (have been checked, they produce similar result with summary)
summary_statistics <- matrix(NA, nrow = n_iterations, ncol = 4*3) # 4 summary statistics*3 parameters
l1_samples <- numeric(n_iterations)
l2_samples <- numeric(n_iterations)
l3_samples <- numeric(n_iterations)

asym_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
xmid_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
scal_summary <- matrix(NA, nrow = n_iterations, ncol = 4)

# Iterate for bootstrap resampling
for (i in 1:n_iterations) {
  # Sample with replacement from the original dataset yeah
  sampled_indices <- sample(nrow(dat_anal), replace = TRUE)
  sampled_data <- dat_anal[sampled_indices, ]
  
  fit <- nls(y ~ SSlogis(x, Asym, xmid, scal)) # Asym, xmid, scal as Asym/(1+exp((xmid-input)/scal))
  
  # Store coefficients
  l1_samples[i] <- coef(fit)["Asym"]
  l2_samples[i] <- coef(fit)["xmid"]
  l3_samples[i] <- coef(fit)["scal"]
  
  summary_fit <- summary(fit)
  
  asym_summary <- summary_fit$parameters[1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  xmid_summary <- summary_fit$parameters[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  scal_summary <- summary_fit$parameters[3, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  
  summary_statistics[i, ] <- c(asym_summary, xmid_summary, scal_summary)
}

summary_fit

# For model comparisons:
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal))
fit_ALL_glmq <- fit

# Compute mean coefficients one-by-one
mean_summary_statistics <- colMeans(summary_statistics, na.rm = TRUE)
mean_summary_statistics # A mixed vector to store Asym & smid & scal--> estimate, std. error, t value and Pr(>|t|)

mean_l1 <- mean(l1_samples)
mean_l2 <- mean(l2_samples)
mean_l3 <- mean(l3_samples)

# PLOT! ########################################################################
# Original data:
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     xlim = c(0,10), ylim = c(0,1), col = "grey50", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles (basically can use plot OR points if layered ehe)
points(mfE_cut_midpoints, p, col = "black", pch = 16,
       # xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
       xlim = c(0,10), ylim = c(0,1), cex = Dissect_Mosq/100
)

# logit function for curve
logit <- function(x,l1,l2,l3) { # Tranlated from Asym/(1+exp((xmid-input)/scal))
  l1/(1+exp((l2-x)/l3))
}
curve(logit(x, mean_l1, mean_l2, mean_l3), col = "grey10", add = TRUE, lwd = 2)

# Add CI using Clopper-Pearson method
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
pred <- predict(fit, newdata = data.frame(x = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(pred*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) +1.96*upper_bound -logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "grey35", lty = 5, lwd = 2)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) -1.96*lower_bound +logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "grey35", lty = 5, lwd = 2)

# 2.2.2. Use mfE vs. Proportion, 9 breakpoints, SPECIES SEPARATED ##############
# 2.2.2.1. Use mfE vs. Proportion, 9 breakpoints, An. gambiae ##################
dat_anal <- dat_sel %>% 
  select(Mosquito_species, Mosquito_larvae_infected_proportion_fromtotaldissected, mfE_mean_calculated, Mosquito_totaldissected) %>% 
  filter(Mosquito_species == "An. gambiae",
         !is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed") %>%
  # view() %>% 
  glimpse()

# mfE mean vs. proportion of mosquitoes with larvae
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

plot(dat_anal$mfE_mean_calculated, dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with Larvae",
     color = "grey20", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles
quantiles <- quantile(x, probs=seq(0,1,length=10))
quantiles
mfE_cut <- cut(x, breaks=quantiles, include.lowest=TRUE)
levels(mfE_cut)
p <- tapply(y, mfE_cut, mean)
p

quantiles_Mosq <- quantile(AllMosq, probs=seq(0,1,length=10))
quantiles_Mosq
Dissect_Mosq <- tapply(AllMosq, mfE_cut, mean)
Dissect_Mosq

mfE_cut_midpoints = rep(NA,9) # length 9
for (i in 1:9) {
  mfE_cut_midpoints[i] = (quantiles[i] + quantiles[i+1])/2
}
mfE_cut_midpoints

# Let's fit this: fit <- glm(y ~ x, family=binomial) by iterations
# Set the number of iterations
n_iterations <- 500

# Create empty vectors & matrices to store coefficients (have been checked, they produce similar result with summary)
summary_statistics <- matrix(NA, nrow = n_iterations, ncol = 4*3) # 4 summary statistics*3 parameters
l1_samples <- numeric(n_iterations)
l2_samples <- numeric(n_iterations)
l3_samples <- numeric(n_iterations)

asym_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
xmid_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
scal_summary <- matrix(NA, nrow = n_iterations, ncol = 4)

# Iterate for bootstrap resampling
for (i in 1:n_iterations) {
  # Sample with replacement from the original dataset yeah
  sampled_indices <- sample(nrow(dat_anal), replace = TRUE)
  sampled_data <- dat_anal[sampled_indices, ]
  
  fit <- nls(y ~ SSlogis(x, Asym, xmid, scal)) # Asym, xmid, scal as Asym/(1+exp((xmid-input)/scal))
  
  # Store coefficients
  l1_samples[i] <- coef(fit)["Asym"]
  l2_samples[i] <- coef(fit)["xmid"]
  l3_samples[i] <- coef(fit)["scal"]
  
  summary_fit <- summary(fit)
  
  asym_summary <- summary_fit$parameters[1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  xmid_summary <- summary_fit$parameters[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  scal_summary <- summary_fit$parameters[3, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  
  summary_statistics[i, ] <- c(asym_summary, xmid_summary, scal_summary)
}

summary_fit

# For model comparisons:
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal))
fit_Angam_glmq <- fit

# Compute mean coefficients one-by-one
mean_summary_statistics <- colMeans(summary_statistics, na.rm = TRUE)
mean_summary_statistics # A mixed vector to store Asym & smid & scal--> estimate, std. error, t value and Pr(>|t|)

mean_l1 <- mean(l1_samples)
mean_l2 <- mean(l2_samples)
mean_l3 <- mean(l3_samples)

# PLOT! ########################################################################
# Original data:
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     # main = expression("Proportion of "*italic("An. gambiae")*" with established larvae"),
     main = expression(italic("An. gambiae")),
     xlim = c(0,10), ylim = c(0,1), col = "grey50", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles (basically can use plot OR points if layered ehe)
points(mfE_cut_midpoints, p, col = "red", pch = 16,
       # xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
       xlim = c(0,10), ylim = c(0,1), cex = Dissect_Mosq/100
)

# logit function for curve
logit <- function(x,l1,l2,l3) { # Tranlated from Asym/(1+exp((xmid-input)/scal))
  l1/(1+exp((l2-x)/l3))
}
curve(logit(x, mean_l1, mean_l2, mean_l3), col = "darkred", add = TRUE, lwd = 2)

# Add CI using Clopper-Pearson method
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
pred <- predict(fit, newdata = data.frame(x = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(pred*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) +1.96*upper_bound -logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "violetred1", lty = 5, lwd = 2)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) -1.96*lower_bound +logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "violetred1", lty = 5, lwd = 2)

# 2.2.2.2. Use mfE vs. Proportion, 9 breakpoints, An. arabiensis ###############
dat_anal <- dat_sel %>% 
  select(Mosquito_species, Mosquito_larvae_infected_proportion_fromtotaldissected, mfE_mean_calculated, Mosquito_totaldissected) %>% 
  filter(Mosquito_species == "An. arabiensis",
         !is.na(Mosquito_larvae_infected_proportion_fromtotaldissected) & Mosquito_larvae_infected_proportion_fromtotaldissected != "not_analysed",
         !is.na(mfE_mean_calculated) & mfE_mean_calculated != "not_analysed") %>%
  # view() %>% 
  glimpse()

# mfE mean vs. proportion of mosquitoes with larvae
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
AllMosq <- dat_anal$Mosquito_totaldissected

plot(dat_anal$mfE_mean_calculated, dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with Larvae",
     color = "grey20", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles
quantiles <- quantile(x, probs=seq(0,1,length=10))
quantiles
mfE_cut <- cut(x, breaks=quantiles, include.lowest=TRUE)
levels(mfE_cut)
p <- tapply(y, mfE_cut, mean)
p

quantiles_Mosq <- quantile(AllMosq, probs=seq(0,1,length=10))
quantiles_Mosq
Dissect_Mosq <- tapply(AllMosq, mfE_cut, mean)
Dissect_Mosq

mfE_cut_midpoints = rep(NA,9) # length 9
for (i in 1:9) {
  mfE_cut_midpoints[i] = (quantiles[i] + quantiles[i+1])/2
}
mfE_cut_midpoints

# Let's fit this: fit <- glm(y ~ x, family=binomial) by iterations
# Set the number of iterations
n_iterations <- 500

# Create empty vectors & matrices to store coefficients (have been checked, they produce similar result with summary)
summary_statistics <- matrix(NA, nrow = n_iterations, ncol = 4*3) # 4 summary statistics*3 parameters
l1_samples <- numeric(n_iterations)
l2_samples <- numeric(n_iterations)
l3_samples <- numeric(n_iterations)

asym_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
xmid_summary <- matrix(NA, nrow = n_iterations, ncol = 4)
scal_summary <- matrix(NA, nrow = n_iterations, ncol = 4)

# Iterate for bootstrap resampling
for (i in 1:n_iterations) {
  # Sample with replacement from the original dataset yeah
  sampled_indices <- sample(nrow(dat_anal), replace = TRUE)
  sampled_data <- dat_anal[sampled_indices, ]
  
  fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), # Asym, xmid, scal as Asym/(1+exp((xmid-input)/scal))
             control=nls.control(maxiter = 1000)) # Control coz' I dunno why maxiter error occur with An. arabiensis data -_-)
  
  # Store coefficients
  l1_samples[i] <- coef(fit)["Asym"]
  l2_samples[i] <- coef(fit)["xmid"]
  l3_samples[i] <- coef(fit)["scal"]
  
  summary_fit <- summary(fit)
  
  asym_summary <- summary_fit$parameters[1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  xmid_summary <- summary_fit$parameters[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  scal_summary <- summary_fit$parameters[3, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
  
  summary_statistics[i, ] <- c(asym_summary, xmid_summary, scal_summary)
}

summary_fit

# For model comparisons:
fit <- nls(y ~ SSlogis(x, Asym, xmid, scal), # Asym, xmid, scal as Asym/(1+exp((xmid-input)/scal))
           control=nls.control(maxiter = 1000)) # Control coz' I dunno why maxiter error occur with An. arabiensis data -_-)
fit_Anara_glmq <- fit

# Compute mean coefficients one-by-one
mean_summary_statistics <- colMeans(summary_statistics, na.rm = TRUE)
mean_summary_statistics # A mixed vector to store Asym & smid & scal--> estimate, std. error, t value and Pr(>|t|)

mean_l1 <- mean(l1_samples)
mean_l2 <- mean(l2_samples)
mean_l3 <- mean(l3_samples)

# PLOT! ########################################################################
# Original data:
y <- dat_anal$Mosquito_larvae_infected_proportion_fromtotaldissected
x <- dat_anal$mfE_mean_calculated
plot(x,y, las = 1, pch = 1,
     xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
     # main = expression("Proportion of "*italic("An. arabiensis")*" with established larvae"),
     main = expression(italic("An. arabiensis")),
     xlim = c(0,10), ylim = c(0,1), col = "grey50", cex = dat_anal$Mosquito_totaldissected/80)

# Quantiles (basically can use plot OR points if layered ehe)
points(mfE_cut_midpoints, p, col = "blue", pch = 16,
       # xlab = "Average of ingested mf in mosquitoes", ylab = "Proportion of dissected mosquitoes with established larvae",
       xlim = c(0,10), ylim = c(0,1), cex = Dissect_Mosq/100
)

# logit function for curve
logit <- function(x,l1,l2,l3) { # Tranlated from Asym/(1+exp((xmid-input)/scal))
  l1/(1+exp((l2-x)/l3))
}
curve(logit(x, mean_l1, mean_l2, mean_l3), col = "darkblue", add = TRUE, lwd = 2)

# Add CI using Clopper-Pearson method
alpha <- 0.05  # significance level
n <- length(y)  # number of observations
x_new <- seq(0, 10, length.out = 100) #change min(x), max(x) to 0 and 10
pred <- predict(fit, newdata = data.frame(x = x_new), type = "response", se.fit = TRUE)
pred_counts <- round(pred*n) # Coz' negative errors

# Predicted probabilities to counts
lower_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[1]/n)
upper_bound <- sapply(pred_counts, function(count) binom.test(count, n, conf.level = 1-alpha)$conf.int[2]/n)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) +1.96*upper_bound -logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "steelblue", lty = 5, lwd = 2)
lines(x_new, logit(x_new, mean_l1, mean_l2, mean_l3) -1.96*lower_bound +logit(x_new, mean_l1, mean_l2, mean_l3),
      col = "steelblue", lty = 5, lwd = 2)
