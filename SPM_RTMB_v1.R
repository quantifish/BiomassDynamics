rm(list = ls())

# libraries
library(RTMB)
library(dplyr)
library(tidyr)
library(ggplot2)

# read catch data
load("raw_data/CRA2/catch_CRA2.rda")
head(catch_CRA2)
# calculate total catch 
catch_CRA2

catch <- catch_CRA2 %>%
  mutate(TotalCatch = Commercial + Recreational) %>%
  group_by(Year) %>%
  summarise(Catch = sum(TotalCatch, na.rm = TRUE)) %>%
  ungroup()

ggplot(catch, aes(x = Year, y = Catch)) +
  geom_line(color = "steelblue", linewidth = 2) +
  geom_point(color = "black", size = 2) +
  labs(x = "Year", y = "Catch") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 


# read logbook CPUE index 
load("raw_data/CRA2/logbook_CRA2.rda")
head(logbook_CRA2)

# Extract mean CPUE and SD by year 
CPUE_logbook <- logbook_CRA2 %>%
  separate(Year, into = c("Year", "Season"), sep = "_") %>%
  mutate(Year = as.integer(Year)) %>%
  group_by(Year) %>%
  summarise(
    CPUE = mean(Mean, na.rm = TRUE),
    SD = mean(SD, na.rm = TRUE)
  ) %>%
  ungroup()
head(CPUE_logbook)
tail(CPUE_logbook)

ggplot(CPUE_logbook, aes(x = Year, y = CPUE)) +
  geom_line(color = "red", linewidth = 2) +
  geom_point(color = "black", size = 2) +
  labs(x = "Year", y = "Logbook CPUE") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) 


# read FSU CPUE index 
load("raw_data/CRA2/fsu_CRA2.rda")
head(fsu_CRA2)

# Extract mean CPUE and SD by year 
CPUE_fsu <- fsu_CRA2 %>%
  mutate(Year = as.integer(year)) %>%
  group_by(Year) %>%
  summarise(
    CPUE = mean(Mean, na.rm = TRUE),
    SD = mean(SD, na.rm = TRUE)
  ) %>%
  ungroup()
head(CPUE_fsu)
tail(CPUE_fsu)

# read CELR CPUE index 
load("raw_data/CRA2/celr_CRA2.rda")
head(celr_CRA2)

# Extract mean CPUE and SD by year 
CPUE_celr <- celr_CRA2 %>%
  separate(Year, into = c("Year", "Season"), sep = "_") %>%
  mutate(Year = as.integer(Year)) %>%
  group_by(Year) %>%
  summarise(
    CPUE = mean(Mean, na.rm = TRUE),
    SD = mean(SD, na.rm = TRUE)
  ) %>%
  ungroup()
head(CPUE_celr)
tail(CPUE_celr)

# state-space Surplus Production model for rock lobster 

## Initial guess on parameters for a Pella-Tomlinson model
K <- 1500       # Carrying capacity
r <- 0.35     # Growth rate
q <- 0.001     # catchability 
z <- 1   # With z=1, the PT model is the Shaeffer model
sigma_p <- 0.1 # process error 

N_catch <- length(catch$Year)
B <- numeric(N_catch)

N_cpue <- length(CPUE_logbook$Year)
cpue_pred <- numeric(N_cpue)

## TMB Data
data <- list(
  year = catch$Year,
  catch_obs = catch$Catch,
  
  cpue_log_obs = CPUE_logbook$CPUE,
  cpue_log_year = CPUE_logbook$Year,
  cpue_log_sd = CPUE_logbook$SD,
  
  cpue_fsu_obs = CPUE_fsu$CPUE,
  cpue_fsu_year = CPUE_fsu$Year,
  cpue_fsu_sd = CPUE_fsu$SD,
  
  cpue_celr_obs = CPUE_celr$CPUE,
  cpue_celr_year = CPUE_celr$Year,
  cpue_celr_sd = CPUE_celr$SD
)

## Initial guess on parameters for the fitted Ricker model
parameters0 <- list(
  logB = log(rep(K, N_catch) + 1e-6),
  logr = log(r),
  logK = log(K),
  logq_log = log(q),
  logq_fsu = log(q),
  logq_celr = log(q),
  logsigma_p = log(sigma_p),
  logz = log(1)
)

## Make TMB model 
fixed <- factor(NA)

fun <- function(parms) {
  getAll(data, parms, warn = FALSE)
  
  # Mark CPUEs for OSA
  cpue_log_obs <- OBS(cpue_log_obs)
  cpue_fsu_obs <- OBS(cpue_fsu_obs)
  cpue_celr_obs <- OBS(cpue_celr_obs)
  
  r <- exp(parms$logr)
  K <- exp(parms$logK)
  q_log <- exp(parms$logq_log)
  q_fsu <- exp(parms$logq_fsu)
  q_celr <- exp(parms$logq_celr)
  sigma_p <- exp(parms$logsigma_p)
  z <- exp(parms$logz)
  logB <- parms$logB
  B <- exp(logB)
  
  N <- length(data$year)
  logB_pred <- numeric(N)
  Prod <- numeric(N)
  
  logB_pred[1] <- log(K)
  Prod[1] <- r / z * logB_pred[1] * (1.0 - (logB_pred[1] / K)^z)
  
  for (y in 2:N) {
    B_prev <- B[y - 1]
    Prod[y - 1] <- r / z * B_prev * (1 - (B_prev / K)^z)
    B_t <- B_prev + Prod[y - 1] - data$catch_obs[y - 1] + 1e-6
    logB_pred[y] <- log(B_t)
  }
  
  nll_B <- -sum(dnorm(logB, mean = logB_pred, sd = sigma_p, log = TRUE))
  
  # Prediction and NLL for each index
  cpue_pred_log <- q_log * B[match(data$cpue_log_year, data$year)]
  cpue_pred_fsu <- q_fsu * B[match(data$cpue_fsu_year, data$year)]
  cpue_pred_celr <- q_celr * B[match(data$cpue_celr_year, data$year)]
  
  nll_log <- -sum(dnorm(cpue_log_obs, mean = cpue_pred_log,
                        sd = data$cpue_log_sd, log = TRUE))
  nll_fsu <- -sum(dnorm(cpue_fsu_obs, mean = cpue_pred_fsu,
                        sd = data$cpue_fsu_sd, log = TRUE))
  nll_celr <- -sum(dnorm(cpue_celr_obs, mean = cpue_pred_celr,
                         sd = data$cpue_celr_sd, log = TRUE))
  
  nll <- nll_B + nll_log + nll_fsu + nll_celr
  return(nll)
}


obj <- MakeADFun(fun,  # Function taking a parameter list (or parameter vector) as input.
                 parameters0, # Parameter list (or parameter vector) used by fun
                 random="logB", # Character vector defining the random effect parameters.
                 map=list(logz=fixed))  # List defining how to optionally collect and fix parameters

## Fit the model 
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))    
sd <- sdreport(obj)

# Recursive quantile residuals: For this to work, you need to mark the observation inside the objective function using
# the OBS function. Thereafter, residual calculation is as simple as oneStepPredict(obj). 
# However, you probably want specify a method to use.

pred  <- oneStepPredict(obj,
                        method="oneStepGeneric",discrete=TRUE,range=c(0,Inf),
                        conditional=1, ## Skip first residual
                        parallel=FALSE)
# CHEKC: Warning message: In qnorm(pred$Fx - U * pred$px) : NaNs produced
pred$residual # there NaNs here 
pred_pars <- exp(opt$par)
names(obj$env)

pred_B <- data.frame(Year=data$year, B = exp(obj$env$last.par.best[names(obj$env$last.par.best) == "logB"]))

ggplot(pred_B, aes(x = Year, y = B)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue", size = 1) + 
  labs(x = "Year", y = "B") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) 

# Predicted CPUE values
cpue_pred_log <- exp(opt$par["logq_log"]) * pred_B$B[match(CPUE_logbook$Year, pred_B$Year)]
cpue_pred_fsu <- exp(opt$par["logq_fsu"]) * pred_B$B[match(CPUE_fsu$Year, pred_B$Year)]
cpue_pred_celr <- exp(opt$par["logq_celr"]) * pred_B$B[match(CPUE_celr$Year, pred_B$Year)]

# Combine into a long data frame
cpue_df <- bind_rows(
  tibble(Year = CPUE_logbook$Year, Observed = CPUE_logbook$CPUE,
         Predicted = cpue_pred_log, Index = "Logbook"),
  tibble(Year = CPUE_fsu$Year, Observed = CPUE_fsu$CPUE,
         Predicted = cpue_pred_fsu, Index = "FSU"),
  tibble(Year = CPUE_celr$Year, Observed = CPUE_celr$CPUE,
         Predicted = cpue_pred_celr, Index = "CELR")
)

# Plot Observed vs Predicted
cpue_df$Index <- factor(cpue_df$Index, levels = c("CELR", "FSU", "Logbook"))

ggplot(cpue_df, aes(x = Year)) +
  geom_point(aes(y = Observed), color = "grey", size = 2) +
  geom_line(aes(y = Predicted), color = "red", linewidth = 1.2) +
  facet_wrap(~Index, ncol = 1) +  # 3 rows
  labs(y = "CPUE") +
  theme_minimal()

resid_fsu <- oneStepPredict(obj = obj, observation.name = "cpue_fsu_obs", method = "oneStepGeneric", trace = FALSE)$residual
resid_celr <- oneStepPredict(obj = obj, observation.name = "cpue_celr_obs", method = "oneStepGeneric", trace = FALSE)$residual
resid_log <- oneStepPredict(obj = obj, observation.name = "cpue_log_obs", method = "oneStepGeneric", trace = FALSE)$residual

df_fsu <- data.frame(Year = data$cpue_fsu_year, Residual = resid_fsu, Index = "FSU")
df_celr <- data.frame(Year = data$cpue_celr_year, Residual = resid_celr, Index = "CELR")
df_log  <- data.frame(Year = data$cpue_log_year, Residual = resid_log, Index = "LB")

resid_df <- bind_rows(df_fsu, df_celr, df_log)
resid_df$Index <- factor(resid_df$Index, levels = c("CELR", "FSU", "LB"))

ggplot(resid_df, aes(x = Year, y = Residual)) +
  geom_hline(yintercept = 0, color = "#00BA38", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = c(-2, 2), color = "#F8766D", linetype = "dashed", linewidth = 0.5) +
  geom_linerange(aes(ymin = 0, ymax = Residual), color = "#619CFF", linewidth = 0.8) +
  geom_point(color = "#619CFF", size = 2) +
  facet_wrap(~ Index, ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(y = "OSA Residuals", x = "Year")
