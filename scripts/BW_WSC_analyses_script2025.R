
####################################################################################################################################
####### Blue whiting and Wind Stress Curl Analysis ----
####################################################################################################################################


#######  Required packages ----
library(glmmTMB)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(lmtest)
library(magrittr)
library(tibble)
library(performance)
library(stats)
library(ggpubr)
library(zoo)
library(ggrepel)
library(scales)
library(knitr)
library(kableExtra)
library(forecast)
library(cowplot)
library(reshape2)
library(stockassessment)
library(MuMIn)
library(RColorBrewer)
library(ggtext)
library(forcats)
library(ggpattern)
library(patchwork)
library(sjstats)
library(ggpmisc)
library(broom)
library(ggeffects)





#######  Set working directory ----
#setwd("")

#######  Load data: BW_WSC_data_19812023 ----
                      # year (1981-2023)
                      # R: recruitment in thousands; shifted to birth-year (1 year back)
                      # SSB: spawning stock biomass in tonnes
                      # logSI: natural log of Survival Index (SI = R/SSB)
                      # CANUM1: catch numbers of age 1; shifted back to birth-year
                      # logCANUM1: natural log of catch numbers of age 1
                      # WSClag1: WSC annual avg. (Jan-Dec) lagged 1 year ahead of spawning, e.g. WSC in 1980 shifted to year 1981
                      # WSClag1_MAM: WSC spawning season avg. (March-April-May) lagged 1 year ahead of spawning



d <- read.xlsx("Data/BW_WSC_data_19812023.xlsx")
head(d)

###########################################################################################
################################  RECRUITMENT models  #####################################  ##########################################----
###########################################################################################
#######  MODELS FIT (logR) ----
# Re-scale covariates to get numerical stability in fits ----
d$SSB_2 <- d$SSB/1e6
d$WSClag1_2 <- d$WSClag1*1e6
# Prepare data frame for model fits ----
d$logR = log(d$R)
d$logSSB = log(d$SSB)
# A1 term
d$group = factor(1)
d$time = factor(d$year)


# Fit models: RI, E, A1, RE, REA1 ----
# Note: 
# m = glmmTMB(logSI ~ SSB + WSClag1 + ar1(time+0|group) , dispformula=~0, data = exdat)
# m2 = glmmTMB(logR ~ SSB + WSClag1 + ar1(time+0|group) + offset(logSSB), dispformula=~0, data = exdat)
# AIC(m,m2)
# Note: they are the same model (REA1)

# Initiate empty list where to store the models
models_logR <- list()

## 1) Ricker model (Ri): log(SI) ~ SSB
models_logR$RI <- glmmTMB( logR ~ SSB_2 + offset(logSSB),data=d)

## 2) Environmental model (E): log(SI) ~ WSClag1 
models_logR$E <- glmmTMB( logR ~ WSClag1_2+ offset(logSSB),data=d)

## 3) First-Order Autoregressive model (A1): log(SI) ~ AR1
models_logR$A1 <- glmmTMB(logR ~ ar1(time + 0 | group) + offset(logSSB), dispformula = ~0, data = d)

## 4) Ricker model + Environmental term (RE): log(SI) ~ SSB + WSClag1
models_logR$RE <- glmmTMB( logR ~ SSB_2 + WSClag1_2 + offset(logSSB),data=d)

## 5) First-Order Autoregressive model (REA1): log(SI) ~ SSB + WSClag1 + AR1
models_logR$REA1 <- glmmTMB(logR ~ SSB_2 + WSClag1_2 + ar1(time+0|group) + offset(logSSB), dispformula = ~0,data = d)


# Models Summary ----
lapply(models_logR, summary)

# Equation REA1 ----
# Extract the fixed effects coefficients: Intercept (ln(𝑎), SSB coefficient (−𝑏), WSC coefficient (c)
fixef(models_logR$REA1)$cond
# Extract AR(1) correlation parameter (rho)
rho <- as.numeric(VarCorr(models_logR$REA1)$cond$group[1,2])
# Extract standard deviation of the AR(1) white noise component (omega_y) (Std. Dev. in model summary)
sigma_omega <- as.numeric(VarCorr(models_logR$REA1)$cond$group[1, 1])^0.5



# Baseline model (GM): logSI ~ long-term mean(logSI) (1981-2023) ----
####### Calculate long-term mean logSI from 1981 to 2023
d$meanlogSI <-  mean(d$logSI)
# Fit the GM model using the  mean
GM <- lm(logR ~ meanlogSI + offset(logSSB), data = d)
summary(GM) # Note: Intercept = meanlogSI. Slope = NA.



# MODEL PERFORMANCE METRICS: n_parms, df, AIC, AICc, LogLik, DeltaAICc, ER ----
### Note: number of parameters (fixed + random effects)
### Note: ER relative to worst model (GM)

# Compute values for GM model (GM is outside model list)
gm_logLik <- logLik(GM)[1]
gm_df <- attr(logLik(GM), "df")
gm_AIC <- AIC(GM)
gm_n_params <- 1  # Only intercept
gm_AICc <- gm_AIC + (2 * gm_df * (gm_df + 1)) / (nrow(d) - gm_df - 1)

# Initialize an empty dataframe to store model metrics   
model_metrics_logR <- data.frame(
  Model = character(),
  n_params = numeric(),
  df = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  LogLik = numeric(),
  Delta_AICc = numeric(),
  ER = numeric(),
  stringsAsFactors = FALSE
)

# Loop over models in list   
for (model_name in names(models_logR)) {
  model_fit <- models_logR[[model_name]]
  
  model_AIC <- AIC(model_fit)
  model_LogLik <- logLik(model_fit)[1]
  model_df <- attr(logLik(model_fit), "df")
  
  n_fixed <- length(fixef(model_fit)$cond)
  n_random <- length(VarCorr(model_fit)$cond)
  n_params <- n_fixed + n_random
  
  model_AICc <- model_AIC + (2 * model_df * (model_df + 1)) / (nrow(d) - model_df - 1)
  
  model_metrics_logR <- rbind(model_metrics_logR, data.frame(
    Model = model_name,
    n_params = n_params,
    df = model_df,
    AIC = model_AIC,
    AICc = model_AICc,
    LogLik = model_LogLik,
    Delta_AICc = NA,
    ER = NA
  ))
}

# Add GM model as separate row 
model_metrics_logR <- rbind(model_metrics_logR, data.frame(
  Model = "GM",
  n_params = gm_n_params,
  df = gm_df,
  AIC = gm_AIC,
  AICc = gm_AICc,
  LogLik = gm_logLik,
  Delta_AICc = NA,
  ER = NA
))

# Sort by AICc
model_metrics_logR <- model_metrics_logR %>% arrange(desc(AICc))

# Compute Delta AICc and Evidence Ratio relative to worst (highest AICc)
worst_AICc <- max(model_metrics_logR$AICc, na.rm = TRUE)
model_metrics_logR$Delta_AICc <- model_metrics_logR$AICc - min(model_metrics_logR$AICc)
model_metrics_logR$ER <- round(exp(0.5 * (worst_AICc - model_metrics_logR$AICc)))

# Round numeric columns (except ER)
model_metrics_logR <- model_metrics_logR %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(ER = as.integer(ER))

# Final model performance summary
print(model_metrics_logR)


####### PREDICTIVE SKILL ASSESSMENT of logR (6 models (GM, RI, E, A1, RE, REA1) + REA1 with WSC min/max scenarios) ----
### Prepare prediction data frame
preds_logR <- d
model_names <- names(models_logR)

# Initialize columns for all model predictions and std. error (SE)
for (model in model_names) {
  preds_logR[[paste0(model, "_predlogR")]] <- NA
  preds_logR[[paste0("SE_", model, "_predlogR")]] <- NA
}

# Rolling geometric mean baseline
preds_logR$GM_predlogR <- NA
preds_logR$SE_GM_predlogR <- NA  

# WSC scenario (max, min)
preds_logR$maxWSC_predlogR <- NA
preds_logR$minWSC_predlogR <- NA

# Determine WSC scenario values
maxWSC <- max(preds_logR$WSClag1_2, na.rm = TRUE)
minWSC <- min(preds_logR$WSClag1_2, na.rm = TRUE)

### Loop: time-series cross-validation from year 1991 onward
# Fit each model in models_logR on past data (tmp) 
# Store predictions and SEs for each model on the next year (nd) 
for (i in 11:nrow(preds_logR)) {
  
  tmp <- preds_logR[1:(i - 1), ]
  nd  <- preds_logR[i, , drop = FALSE]
  
  # Loop over all models
  for (model in model_names) {
    mod_obj <- models_logR[[model]]
    formula_mod <- as.formula(mod_obj$call$formula)
    disp_mod <- if (!is.null(mod_obj$call$dispformula)) as.formula(mod_obj$call$dispformula) else ~1
    
    fit <- glmmTMB(formula = formula_mod,         
                   dispformula = disp_mod,
                   data = tmp)
    
    pr <- predict(fit, newdata = nd, se.fit = TRUE, allow.new.levels = TRUE)
    
    preds_logR[i, paste0(model, "_predlogR")] <- pr$fit
    preds_logR[i, paste0("SE_", model, "_predlogR")] <- pr$se.fit
  }
  
  # Rolling geometric mean prediction: log(GM(R)) = mean(logR)
  # Rolling standard error of logR
  preds_logR$GM_predlogR[i] <- mean(tmp$logR, na.rm = TRUE)
  n <- sum(!is.na(tmp$logR))
  if (n > 1) {
    preds_logR$SE_GM_predlogR[i] <- sd(tmp$logR, na.rm = TRUE) / sqrt(n)
  }
  
  # Fit REA1 model on current training set for scenario predictions
  mod_REA1 <- models_logR$REA1
  formula_REA1 <- as.formula(mod_REA1$call$formula)
  disp_REA1 <- if (!is.null(mod_REA1$call$dispformula)) as.formula(mod_REA1$call$dispformula) else ~1
  
  fit_REA1 <- glmmTMB(formula = formula_REA1,
                      dispformula = disp_REA1,
                      data = tmp)
  
  # Predict logR for max and min WSC scenarios using REA1
  nd_max <- nd
  nd_max$WSClag1_2 <- maxWSC
  preds_logR$maxWSC_predlogR[i] <- predict(fit_REA1, newdata = nd_max, se.fit = TRUE, allow.new.levels = TRUE)$fit
  
  nd_min <- nd
  nd_min$WSClag1_2 <- minWSC
  preds_logR$minWSC_predlogR[i] <- predict(fit_REA1, newdata = nd_min, se.fit = TRUE, allow.new.levels = TRUE)$fit
}

###  Confidence intervals for GM-based prediction
preds_logR$GM_logR_upper <- preds_logR$GM_predlogR + 1.96 * preds_logR$SE_GM_predlogR
preds_logR$GM_logR_lower <- preds_logR$GM_predlogR - 1.96 * preds_logR$SE_GM_predlogR

# Coverage percentage (GM and REA1) ----
# Initialize coverage result data frame
coverage_logR <- data.frame(
  Model = c("REA1_predlogR", "GM_predlogR"),
  Coverage_Percentage = NA
)

# Loop over the two models
for (nn in coverage_logR$Model) {
  se_col <- ifelse(nn == "REA1_predlogR", "SE_REA1_predlogR", "SE_GM_predlogR")
  lower <- preds_logR[[nn]] - 1.96 * preds_logR[[se_col]]
  upper <- preds_logR[[nn]] + 1.96 * preds_logR[[se_col]]
  
  # Check if observed logR is within predicted CI
  within <- preds_logR$logR >= lower & preds_logR$logR <= upper
  
  # Calculate coverage percentage
  coverage_logR$Coverage_Percentage[coverage_logR$Model == nn] <- mean(within, na.rm = TRUE) * 100
}


# Extract and round
coverage_REA1   <- round(coverage_logR$Coverage_Percentage[coverage_logR$Model == "REA1_predlogR"], 0)
coverage_GM <- round(coverage_logR$Coverage_Percentage[coverage_logR$Model == "GM_predlogR"], 0)
#######  ONE-STEP AHEAD PREDICTION ERRORS of logR (observed - predicted) ----
# Define the names of the prediction columns you want to compare.
pred_cols <- c("RI_predlogR", "E_predlogR", "A1_predlogR", 
               "RE_predlogR", "REA1_predlogR", "GM_predlogR", 
               "maxWSC_predlogR", "minWSC_predlogR")

# Subset indices corresponding to 1991 onward (assuming row 11 is 1991)
rows_sub <- 11:nrow(preds_logR)

# Create a new data frame with the corresponding years
res_predlogR <- data.frame(year = preds_logR$year[rows_sub])

# Loop over each prediction column and compute the log(R) prediction errors (observed - predicted) 
for (col in pred_cols) {
  res_predlogR[[col]] <- preds_logR$logR[rows_sub] - preds_logR[[col]][rows_sub]
}


# View the resulting data frame
head(res_predlogR)



# MSE, SD, Mean of logR prediction errors ----
# Extract only prediction errors columns (excluding 'year')
res_only <- res_predlogR[, -1]

# Define model names (should match order of columns in res_only)
logR_models <- names(res_only)

# Summarize prediction errors into metrics and order them by descreasing MSE
error_predlogR_metrics <- data.frame(
  Model = logR_models,
  MSE  = sapply(res_only, function(x) mean(x^2, na.rm = TRUE)),
  SD   = sapply(res_only, function(x) sd(x, na.rm = TRUE)),
  Mean = sapply(res_only, function(x) mean(x, na.rm = TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(desc(MSE)) 

# View the result
print(error_predlogR_metrics)



# DIAGNOSTIC CHECK of logR prediction errors (autocorrelation, time trends, normality) ----
# A) Autocorrelation check of logR prediction errors 
# Select only prediction errors columns (excluding year and the two WSC scenario columns)
residual_cols <- res_predlogR[, !(names(res_predlogR) %in% c("year", "maxWSC_predlogR", "minWSC_predlogR"))]
# Reorder columns with GM_predlogR first and rename all columns
residual_cols <- residual_cols[, c("GM_predlogR", "RI_predlogR", "E_predlogR", "RE_predlogR", "A1_predlogR", "REA1_predlogR")]

# Rename the columns
colnames(residual_cols) <- c("GM", "RI", "E", "RE", "A1", "REA1")

# Check the result
head(residual_cols)
# Set up plotting layout
par(mfrow = c(3, 2))  # 3 rows, 2 columns

# Generate ACF plots
for (col in names(residual_cols)) {
  acf(residual_cols[[col]], main = paste(col), na.action = na.pass)
}


# B) Time series of logR prediction errors 
# Set up layout: 3 rows × 2 columns
par(mfrow = c(3, 2))

# Time axis (assuming res_predlogR$year exists)
years <- res_predlogR$year

# Plot time series for each model's prediction errors
for (col in names(residual_cols)) {
  plot(years, residual_cols[[col]], type = "b", pch = 16, col = "steelblue",
       xlab = "Year", ylab = "Prediction error", main = col,
       ylim = range(residual_cols, na.rm = TRUE))
  abline(h = 0, lty = 2, col = "gray40")
}


# C) Normality check of logR prediction errors 
# 1) QQ plots
par(mfrow = c(2, 3))  # Adjust layout to number of models

for (col in names(residual_cols)) {
  qqnorm(residual_cols[[col]], main = paste("QQ Plot:", col), pch = 16, col = "steelblue")
  qqline(residual_cols[[col]], col = "red", lwd = 2)
}

# 2) Density plots
par(mfrow = c(2, 3))  # Same layout

for (col in names(residual_cols)) {
  plot(density(residual_cols[[col]], na.rm = TRUE),
       main = paste("Density Plot:", col),
       xlab = "Prediction error", col = "steelblue", lwd = 2)
  rug(residual_cols[[col]], col = "gray40")
}

# 3) Shapiro test 
# Perform Shapiro-Wilk test on each model's prediction errors
shapiro_results <- sapply(residual_cols, function(x) {
  test <- shapiro.test(x)
  c(W = round(test$statistic, 3), p_value = round(test$p.value, 4))
})

# Convert matrix to data frame and format
shapiro_df <- as.data.frame(t(shapiro_results))  # Transpose and convert to data.frame
shapiro_df$Model <- rownames(shapiro_df)         # Add model names as column
rownames(shapiro_df) <- NULL                     # Remove rownames

# Rename columns correctly
colnames(shapiro_df)[1:2] <- c("W", "p_value")

# Reorder columns
shapiro_df <- shapiro_df[, c("Model", "W", "p_value")]

# Print result
print(shapiro_df)





#######  STOCK ASSESSEMENT MODEL (SAM)  ----
# Load SAM output for BW data:  R, SSB, Fbar with CIs estimated using original SAM model
load("GE/SAM_BW_2024.RData")  
# Check data
ls()
class(origfit)
str(origfit, max.level = 1)
names(origfit)
plot(origfit)  # Recruitment, SSB, Fbar with CIs 
summary(origfit)   # Gives parameter estimates and SEs
fit <- origfit

# Deterministic projections with no stochasticity (no simulation-based noise)
nosim=2

# Load myforecast() function from SAM library
source("GE/myforecast.R")
#Custom function for forecasting future SSB values using the SAM model output under fixed recruitment scenarios.
#Inputs: a) SAM model output, b) Recruitment value to impose.

# Prepare exd dataframe
# Create exd with only necessary columns for SSB forecast (REA1, GM, WSC max/min)
exd <- preds_logR %>%
  dplyr::select(
    year, WSClag1_2, R, SSB_2, group, time,
    logR, logSSB,
    REA1_predlogR,           # REA1 model
    GM_predlogR,           # Geometric mean
    maxWSC_predlogR,         # Max WSC scenario
    minWSC_predlogR)          # Min WSC scenario
exd$fcSSBmeanR = NA
exd$fcSSBwsc = NA
exd$fcSSBwscmin = NA
exd$fcSSBwscmax = NA

# Get predicted logR from the original SAM output and plot it
rsel <- which(names(fit$sdrep$value) == "logR")
logR_originalSAM <- fit$sdrep$value[rsel]
plot(fit$data$years, logR_originalSAM, type = "b")
## Note: OBS, exd years are lagged (1981:2021 correspond to 1982:2022 in SAM) 

# FORECAST SSB for 4 log(R) prediction scenarios (GM, REA1, WSC max, WSC min) ----
# For each year (yy) from 2006-2024: 
for (yy in 2006:2024) {
  # Re-fit SAM excluding years yy-1 to 2024 and 
  fit <- runwithout(origfit, year = (yy:2024))
  
  # And predict logR in year yy - 1 under different assumptions (models)
  r1 <- exp(exd$REA1_predlogR[exd$year == yy - 1])         # REA1
  r2 <- exp(exd$GM_predlogR[exd$year == yy - 1])         # Historical mean 
  r3 <- exp(exd$minWSC_predlogR[exd$year == yy - 1])       # Log WSC (min)
  r4 <- exp(exd$maxWSC_predlogR[exd$year == yy - 1])       # High WSC (max)
  
  # And forecast SSB for each R scenario (myforecast()) 
  fc1 <- myforecast(fit, fscale = c(1, 1, 1), nosim = nosim, deterministic = TRUE, recruitmentValue = r1)
  fc2 <- myforecast(fit, fscale = c(1, 1, 1), nosim = nosim, deterministic = TRUE, recruitmentValue = r2)
  fc3 <- myforecast(fit, fscale = c(1, 1, 1), nosim = nosim, deterministic = TRUE, recruitmentValue = r3)
  fc4 <- myforecast(fit, fscale = c(1, 1, 1), nosim = nosim, deterministic = TRUE, recruitmentValue = r4)
  
  # Store forecasted SSB values in exd
  exd$fcSSBwsc[exd$year == yy - 1]     <- attr(fc1, "shorttab")[3, 3]
  exd$fcSSBmeanR[exd$year == yy - 1] <- attr(fc2, "shorttab")[3, 3]
  exd$fcSSBwscmin[exd$year == yy - 1]  <- attr(fc3, "shorttab")[3, 3]
  exd$fcSSBwscmax[exd$year == yy - 1]  <- attr(fc4, "shorttab")[3, 3]
  
  # Plot predicted logR for each model 
  # Note: truncated lines because model stops at year y - 1 and re-estimates logR up to that year 
  #   Retrospective logR trajectories: what the model would have estimated if it had stopped in y-1). 
  # E.g.: for y = 2006 → fit goes only up to 2005 → line ends at 2005
  # E.g.: for y = 2024 → fit goes  up to 2023 → line ends at 2023
  rsel <- which(names(fit$sdrep$value) == "logR")
  ssbsel <- which(names(fit$sdrep$value) == "logssb")
  
  logR_originalSAM <- fit$sdrep$value[rsel]
  SSB_originalSAM  <- fit$sdrep$value[ssbsel]
  lines(fit$data$years, logR_originalSAM, col = yy %% 10)
}
## Note: last yy (2022) means fit has 1981-2021, last WSC correspond to 2022


# SSB forecast deviation (%) ----
# Calculate SSB percentage deviations from mean forecast (baseline = GM)
exd <- exd %>%
  mutate(
    pctDev     = -(fcSSBmeanR - fcSSBwsc)     / fcSSBmeanR,  # REA1
    pctDevmin  = -(fcSSBmeanR - fcSSBwscmin)  / fcSSBmeanR,  # WSCmin
    pctDevmax  = -(fcSSBmeanR - fcSSBwscmax)  / fcSSBmeanR   # WSCmax
  )

# Pivot to long format for plotting — only from 2005 onward
ssb_devs <- exd %>%
  select(year, pctDev, pctDevmin, pctDevmax) %>%
  rename(
    REA1   = pctDev,
    WSCmin = pctDevmin,
    WSCmax = pctDevmax
  ) %>%
  pivot_longer(cols = -year, names_to = "Model", values_to = "Deviation") %>%
  filter(year >= 2005) %>%
  mutate(Deviation = Deviation * 100)


# Mean ± SD, max, min deviation for REA1, WSCmin, WSCmax  ----
ssb_devs_summary <- ssb_devs %>%
  group_by(Model) %>%
  summarise(
    Mean = mean(Deviation, na.rm = TRUE),
    SD   = sd(Deviation, na.rm = TRUE),
    Min  = min(Deviation, na.rm = TRUE),
    Max  = max(Deviation, na.rm = TRUE)
  ) %>%
  arrange(Model)

print(ssb_devs_summary)






####################################################################################################
####################################################################################################




###############################################################################################
################################  CANUM1 models ############################################### ---- 
###############################################################################################
#######  MODELS FIT (logCANUM1) ----
# Fit models: A1, E, EA1 ----
# Initiate empty list to store the models
models_CANUM1 <- list()

## 1) A1 model: CANUM1_2 ~ AR1
models_CANUM1$A1 <- glmmTMB(logCANUM1 ~ ar1(time + 0 | group), dispformula = ~0, data = d)

## 2) E model: CANUM1_2 ~ WSClag1_2
models_CANUM1$E <- glmmTMB(logCANUM1 ~ WSClag1_2, data = d)

## 3) EA1 model: CANUM1_2 ~ WSClag1_2 + AR1
models_CANUM1$EA1 <- glmmTMB(logCANUM1 ~ WSClag1_2 + ar1(time + 0 | group), dispformula = ~0, data = d)



# Models Summary ----
lapply(models_CANUM1, summary)

# MODEL PERFORMANCE METRICS: n_params, df, AIC, AICc, LogLik, DeltaAICc, ER ----
### Note: number of parameters (fixed + random effects)
### Note: ER relative to worst model (A1)

# Initialize an empty dataframe
model_metrics_CANUM1 <- data.frame(
  Model = character(),
  n_params = numeric(),
  df = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  LogLik = numeric(),
  Delta_AICc = numeric(),
  ER = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model in the list and extract values
for (model_name in names(models_CANUM1)) {
  model_fit <- models_CANUM1[[model_name]]
  
  model_AIC <- AIC(model_fit)
  model_LogLik <- logLik(model_fit)[1]
  model_df <- attr(logLik(model_fit), "df")
  
  n_fixed <- length(fixef(model_fit)$cond)
  n_random <- length(VarCorr(model_fit)$cond)
  n_params <- n_fixed + n_random
  
  # Calculate AICc
  AICc_value <- model_AIC + (2 * model_df * (model_df + 1)) / (nrow(d) - model_df - 1)
  
  model_metrics_CANUM1 <- rbind(model_metrics_CANUM1, data.frame(
    Model = model_name,
    n_params = n_params,
    df = model_df,
    AIC = model_AIC,
    AICc = AICc_value,
    LogLik = model_LogLik,
    Delta_AICc = NA,
    ER = NA
  ))
}

# Sort models by AICc (descending → worst first)
model_metrics_CANUM1 <- model_metrics_CANUM1[order(-model_metrics_CANUM1$AICc), ]

# Identify worst model
worst_AICc <- max(model_metrics_CANUM1$AICc)

# Calculate ΔAICc and Evidence Ratios
model_metrics_CANUM1$Delta_AICc <- model_metrics_CANUM1$AICc - min(model_metrics_CANUM1$AICc)
model_metrics_CANUM1$ER <- round(exp(0.5 * (worst_AICc - model_metrics_CANUM1$AICc)))

# Round all numeric columns (except ER), and cast ER to integer
model_metrics_CANUM1 <- model_metrics_CANUM1 %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  mutate(ER = as.integer(ER))

# Print final table
print(model_metrics_CANUM1)


####### PREDICTIVE SKILL ASSESSMENT of logCANUM1 (A1, E, EA1) ----

# Step 1: Reserve first 10 years for model fitting
preds_CANUM1 <- data.frame(year = d$year[-c(1:10)])

# Step 2: Initialize prediction columns
for(nn in names(models_CANUM1)) { 
  preds_CANUM1[[nn]] <- NA
  preds_CANUM1[[paste0(nn, ".se")]] <- NA 
}

# Step 3: Run rolling-origin cross-validation
for(i in seq_along(models_CANUM1)) {
  model_name <- names(models_CANUM1)[i]
  cat("model", model_name)
  for(y in preds_CANUM1$year){
    
    tmp <- d[d$year < y, ]
    nd  <- d[d$year == y, ]
    
    if (inherits(models_CANUM1[[i]], "glmmTMB")) {
      tmpm <- glmmTMB(as.list(models_CANUM1[[i]]$call)$formula,
                      dispformula = as.list(models_CANUM1[[i]]$call)$dispformula,
                      data = tmp)
      pr <- predict(tmpm, newdata = nd, se.fit = TRUE, allow.new.levels = TRUE)
      
      preds_CANUM1[preds_CANUM1$year == y, model_name] <- pr$fit
      preds_CANUM1[preds_CANUM1$year == y, paste0(model_name, ".se")] <- pr$se.fit
    }
  }
  cat(" done.\n")
}


# Coverage percentage (E, EA1) ----

coverage_results_CANUM1 <- data.frame(Model = names(models_CANUM1), Coverage_Percentage = NA)

for (nn in names(models_CANUM1)) {
  lower <- preds_CANUM1[[nn]] - 1.96 * preds_CANUM1[[paste0(nn, ".se")]]
  upper <- preds_CANUM1[[nn]] + 1.96 * preds_CANUM1[[paste0(nn, ".se")]]
  
  observed <- d$logCANUM1[d$year %in% preds_CANUM1$year]
  inside <- observed >= lower & observed <= upper
  
  coverage_results_CANUM1$Coverage_Percentage[coverage_results_CANUM1$Model == nn] <- round(mean(inside, na.rm = TRUE) * 100, 0)
}


# Extract values for E and EA1 models
coverage_E   <- coverage_results_CANUM1$Coverage_Percentage[coverage_results_CANUM1$Model == "E"]
coverage_EA1 <- coverage_results_CANUM1$Coverage_Percentage[coverage_results_CANUM1$Model == "EA1"]




#######  ONE-STEP AHEAD PREDICTION ERRORS of logCANUM1 (observed - predicted) ----

# Compute the prediction errors (observed - predicted) for log(CANUM1)
errs_pred_CANUM1 <- list()
for (nn in names(models_CANUM1)) {
  errs_pred_CANUM1[[nn]] <- d$logCANUM1[d$year %in% preds_CANUM1$year] - preds_CANUM1[[nn]]
}


# View the resulting data frame
head(errs_pred_CANUM1)

# MSE, SD, Mean of logCANUM1 prediction errors ----
error_pred_CANUM1_metrics <- data.frame(
  Model = names(errs_pred_CANUM1),
  MSE  = sapply(errs_pred_CANUM1, function(x) mean(x^2, na.rm = TRUE)),
  SD   = sapply(errs_pred_CANUM1, function(x) sd(x, na.rm = TRUE)),
  Mean = sapply(errs_pred_CANUM1, function(x) mean(x, na.rm = TRUE))
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))  # Round all numeric columns to 2 decimals

# View results
print(error_pred_CANUM1_metrics)


# DIAGNOSTIC CHECK of logCANUM1 prediction errors (autocorrelation, time trends, normality) ----

# Convert the list to a data.frame
residuals_CANUM1 <- as.data.frame(errs_pred_CANUM1)
colnames(residuals_CANUM1) <- names(errs_pred_CANUM1)

# Time vector (replace if needed)
years <- res_predlogR$year

# A) Autocorrelation check of logCANUM1 prediction errors 
par(mfrow = c(3, 1))

for (col in names(residuals_CANUM1)) {
  acf(residuals_CANUM1[[col]], main = paste("ACF:", col), na.action = na.pass)
}


# B) Time series of logCANUM1 prediction errors 
par(mfrow = c(3, 1))

for (col in names(residuals_CANUM1)) {
  plot(years, residuals_CANUM1[[col]], type = "b", pch = 16, col = "steelblue",
       xlab = "Year", ylab = "Prediction error", main = col,
       ylim = range(residuals_CANUM1, na.rm = TRUE))
  abline(h = 0, lty = 2, col = "gray40")
}

# C) Normality check of logCANUM1 prediction errors 
# 1) QQ PLOTS 
par(mfrow = c(3, 1))

for (col in names(residuals_CANUM1)) {
  qqnorm(residuals_CANUM1[[col]], main = paste("QQ Plot:", col), pch = 16, col = "steelblue")
  qqline(residuals_CANUM1[[col]], col = "red", lwd = 2)
}


# 2) Density plots
par(mfrow = c(3, 1))

for (col in names(residuals_CANUM1)) {
  plot(density(residuals_CANUM1[[col]], na.rm = TRUE),
       main = paste("Density Plot:", col),
       xlab = "Prediction error", col = "steelblue", lwd = 2)
  rug(residuals_CANUM1[[col]], col = "gray40")
}

# 3) Shapiro test
shapiro_CANUM1 <- sapply(residuals_CANUM1, function(x) {
  test <- shapiro.test(x)
  c(W = round(test$statistic, 3), p_value = round(test$p.value, 4))
})

shapiro_df_CANUM1 <- as.data.frame(t(shapiro_CANUM1))
shapiro_df_CANUM1$Model <- rownames(shapiro_df_CANUM1)
rownames(shapiro_df_CANUM1) <- NULL
colnames(shapiro_df_CANUM1)[1:2] <- c("W", "p_value")
shapiro_df_CANUM1 <- shapiro_df_CANUM1[, c("Model", "W", "p_value")]

print(shapiro_df_CANUM1)





####################################################################################################
####################################################################################################

####################################################################################################
################  WSC annual avg. (Jan-Dec) vs WSC spawning season avg. (March-May) ----
####################################################################################################
#######  MODELS FIT (logR) ----
# Re-load data and store them with different name (dMAM)
dMAM <- read.xlsx("Data/BW_WSC_data_19812023.xlsx")
head(dMAM)
# Rename for simplicity: Y = annual avg. MAM = spawning season avg. ----
names(dMAM)[names(dMAM) == "WSClag1"] <- "WSClag1_Y"
# Create non-lagged WSC_Y and WSC_MAM columns by shifting the lagged columns backward and adding one more year
dMAM$WSC_Y <- c(dMAM$WSClag1_Y[-1], 4.329031e-06)
dMAM$WSC_MAM <- c(dMAM$WSClag1_MAM[-1], -2.156536e-08)

# Re-scale covariates to get numerical stability in fits ----
dMAM$SSB_2 <- dMAM$SSB/1e6
dMAM$WSClag1_Y_2 <- dMAM$WSClag1_Y*1e6
dMAM$WSClag1_MAM_2 <- dMAM$WSClag1_MAM*1e6
dMAM$WSC_Y_2 <- dMAM$WSC_Y*1e6
dMAM$WSC_MAM_2 <- dMAM$WSC_MAM*1e6

# Prepare data frame for model fits ----
dMAM$logR = log(dMAM$R)
dMAM$logSSB = log(dMAM$SSB)
# A1 term
dMAM$group = factor(1)
dMAM$time = factor(dMAM$year)

head(dMAM)



# Fit models: ANNUAL vs SPAWNING; lag0 vs lag1 ----
# Initiate empty list where to store the models
models_WSC <- list()

## 1) Environmental model (E): log(SI) ~ WSClag1 annual avg
models_WSC$ANNUAL_lag1 <- glmmTMB(logR ~ WSClag1_Y_2+ offset(logSSB),data=dMAM)
## 2) Environmental model (E): log(SI) ~ WSClag1 spawning avg
models_WSC$MAM_lag1 <- glmmTMB(logR ~ WSClag1_MAM_2+ offset(logSSB),data=dMAM)
## 3) Environmental model (E): log(SI) ~ WSClag0 annual avg
models_WSC$ANNUAL_lag0 <- glmmTMB(logR ~ WSC_Y_2+ offset(logSSB),data=dMAM)
## 4) Environmental model (E): log(SI) ~ WSClag0 spawning avg
models_WSC$MAM_lag0 <- glmmTMB(logR ~ WSC_MAM_2+ offset(logSSB),data=dMAM)

# Models Summary ----
lapply(models_WSC, summary)
# MODEL PERFORMANCE METRICS: n_parms, df, AIC, AICc, logLik, DeltaAICc, ER ----
# Initialize empty data frame
model_metrics_WSC <- data.frame(
  Model = character(),
  n_params = numeric(),
  df = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  LogLik = numeric(),
  Delta_AICc = numeric(),
  ER = numeric(),
  stringsAsFactors = FALSE
)

# Loop through models in models_WSC
for (model_name in names(models_WSC)) {
  model_fit <- models_WSC[[model_name]]
  
  model_AIC <- AIC(model_fit)
  model_LogLik <- logLik(model_fit)[1]
  model_df <- attr(logLik(model_fit), "df")
  
  n_fixed <- length(fixef(model_fit)$cond)
  n_random <- length(VarCorr(model_fit)$cond)
  n_params <- n_fixed + n_random
  
  model_AICc <- model_AIC + (2 * model_df * (model_df + 1)) / (nrow(dMAM) - model_df - 1)
  
  model_metrics_WSC <- rbind(model_metrics_WSC, data.frame(
    Model = model_name,
    n_params = n_params,
    df = model_df,
    AIC = model_AIC,
    AICc = model_AICc,
    LogLik = model_LogLik,
    Delta_AICc = NA,
    ER = NA
  ))
}

# Compute Delta AICc (relative to best) and ER (relative to worst)
best_AICc <- min(model_metrics_WSC$AICc)
worst_AICc <- max(model_metrics_WSC$AICc)

model_metrics_WSC <- model_metrics_WSC %>%
  mutate(
    Delta_AICc = AICc - best_AICc,
    ER = round(exp(0.5 * (worst_AICc - AICc)))
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(ER = as.integer(ER)) %>%
  arrange(-AICc)

# View model performance table
print(model_metrics_WSC)
