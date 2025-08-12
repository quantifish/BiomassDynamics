rm(list = ls())

# remotes::install_github("kaskr/TMB_contrib_R/TMBhelper")
# remotes::install_github("Cole-Monnahan-NOAA/adnuts", ref = "sparse_M")

library(tidyverse)
library(decamod)
library(RTMB)
library(TMBhelper)

theme_set(theme_bw())

# Load full model ----

load("base.rda")  #CRA 2 base model
obj_full <- obj
names(obj_full$report())

# Read catch data ----

catch <- read_csv("data/CRA2/CRA2_catch_data_input.csv") %>%
  mutate(TotalCatch = Commercial + Recreational + Illegal + Customary) %>%
  # group_by(Year) %>%
  # summarise(Catch = sum(TotalCatch, na.rm = TRUE)) %>%
  # ungroup() %>%
  select(Year, Season, Catch = TotalCatch) %>%
  arrange(Year, Season) %>%
  mutate(Period = 1:n())

ggplot(data = catch, aes(x = Year, y = Catch)) +
  geom_line(color = "steelblue", linewidth = 2) +
  geom_point(color = "black", size = 2) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.05)) +
  facet_wrap(Season ~ .) +
  labs(x = "Fishing year", y = "Catch")

# Read in CPUE data ----

load("data/CRA2/fsu_CRA2.rda")

fsu <- fsu_CRA2 %>%
  mutate(Season = ifelse(season == "AW", 1, 2)) %>%
  select(-1) %>%         
  rename(Year = year)  %>%
  mutate(Year = as.character(Year)) %>%
  mutate(q = ifelse(season == "AW", 1, 2))  %>%
  mutate(Index = "FSU") 

load("data/CRA2/celr_CRA2.rda")

celr <- celr_CRA2 %>%
  separate(Year, into = c("Year", "season"), sep = "_") %>%
  mutate(Season = ifelse(season == "AW", 1, 2)) %>%
  mutate(Year = as.character(Year)) %>%
  mutate(q = ifelse(season == "AW", 3, 4)) %>%
  mutate(Index = "CELR")

logbook <- read_csv("data/CRA2/CRA2_LB_CPUE_legal_2025.csv") %>%
  separate(Year, into = c("Year", "season"), sep = "_") %>%
  mutate(Season = ifelse(season == "AW", 1, 2)) %>%
  mutate(Year = as.character(Year)) %>%
  mutate(q = ifelse(season == "AW", 5, 6)) %>%
  mutate(Index = "Logbook")

ys <- catch %>% mutate(Year = as.character(Year)) %>% select(Year, Season, Period)

CPUE <- bind_rows(fsu, celr, logbook) %>%
  mutate(logSD = log(1 + SD / Mean)) %>%
  left_join(ys, by = join_by(Year, Season)) %>%
  select(Year, Season, Period, Median, logSD, q, Index)
tail(CPUE)
head(CPUE)
l_cpue <- c("FSU",  "CELR",  "Logbook")


# State-space Surplus Production model for rock lobster 

n_time <- length(catch$Period)

data <- list(
  catch_time = catch$Period,
  catch_season = catch$Season,
  catch_obs = catch$Catch,
  cpue_time = CPUE$Period,
  cpue_season = CPUE$Season,
  cpue_q = CPUE$q,
  cpue_obs = CPUE$Median,
  cpue_sd = CPUE$logSD
)

parameters <- list(
  log_r = log(0.35), # Growth rate 
  log_K = log(1500), # Carrying capacity
  log_z = log(1), # With z=1, the PT model is the Shaeffer model
  log_q = log(rep(0.001, 6)), # catchability
  log_cpue_pow = log(rep(1, 6)),
  #log_cpue_pro = log(rep(1e-6, 3)), # process error for CPUE
  log_sigmap = log(0.1), # process error
  log_B = log(rep(1500, n_time - 1))
)

priors <- list()
priors[["log_r"]] <- list(type = "normal", par1 = log(0.35), par2 = 1.5, index = which("log_r" == names(parameters)))
priors[["log_K"]] <- list(type = "normal", par1 = log(1000), par2 = 5, index = which("log_K" == names(parameters)))
priors[["log_z"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_z" == names(parameters)))
priors[["log_q"]] <- list(type = "normal", par1 = -5, par2 = 2, index = which("log_q" == names(parameters)))
#priors[["log_cpue_pro"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_cpue_pro" == names(parameters)))
priors[["log_sigmap"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_sigmap" == names(parameters)))
priors[["log_cpue_pow"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_cpue_pow" == names(parameters)))

priors
evaluate_priors(parameters, priors)
data$priors <- priors

fun <- function(parameters, data) {
  getAll(parameters, data, warn = FALSE)
  r <- exp(log_r) /2 # because it is by season ? 
  K <- exp(log_K)
  z <- exp(log_z)
  cpue_pow <- exp(log_cpue_pow)
  B <- c(K, exp(log_B))
  
  n_time <- length(data$catch_time)
  Prod <- r / z * B * (1 - (B / K)^z)
  Bpred <- numeric(n_time)
  Bpred[1] <- K
  for (t in 2:n_time) {
    Bpred[t] <- B[t - 1] + Prod[t - 1] - catch_obs[t - 1]
  }
  
  #nll_B <- -sum(dlnorm(B, meanlog = log(biomass), sdlog = exp(log_sigmap), log = TRUE))
  nll_B <- -sum(dnorm(log(B), mean = log(Bpred), sd = exp(log_sigmap), log = TRUE))
  REPORT(Bpred)
  
  cpue_obs <- OBS(cpue_obs)
  # cpue_pred <- exp(log_q)[cpue_q] * B[match(cpue_year, year)]^cpue_pro
  cpue_pred <- exp(log_q)[cpue_q] * B[cpue_time]^cpue_pow
  # cpue_pred <- exp(log_q)[cpue_q] * B[cpue_time]
  # cpue_pro <- exp(log_cpue_pro)[cpue_q]
  # cpue_sigma <- sqrt(cpue_sd^2 + exp(log_cpue_pro[cpue_q])^2)
  cpue_sigma <- cpue_sd
  #nll_cpue <- -dlnorm(x = cpue_obs, meanlog = log(cpue_pred), sdlog = cpue_sigma, log = TRUE)
  nll_cpue <- -1 * dlnorm(x = cpue_obs, meanlog = log(cpue_pred), sdlog = cpue_sigma, log = TRUE)
  REPORT(cpue_pred)
  
  nll <- nll_B + sum(nll_cpue) - evaluate_priors(parameters, priors)
  return(nll)
}

# Fit the model ----

map <- list(log_z = factor(NA))
# map <- list(log_z = factor(NA), log_cpue_pro = factor(c(NA, NA, NA)))
# map <- list(log_cpue_pro = factor(c(NA, NA, NA)))

obj <- MakeADFun(cmb(fun, data), parameters, random = "log_B", map = map)

Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_r", names(obj$par))] <- -10
Upr[grep("log_r", names(obj$par))] <- 10
Lwr[grep("log_K", names(obj$par))] <- log(100)
Upr[grep("log_K", names(obj$par))] <- log(10000)
Lwr[grep("log_z", names(obj$par))] <- log(0)
Upr[grep("log_z", names(obj$par))] <- log(2)
Lwr[grep("log_cpue_pow", names(obj$par))] <- log(0)
Upr[grep("log_cpue_pow", names(obj$par))] <- log(2)
Lwr[grep("log_q", names(obj$par))] <- -25
Upr[grep("log_q", names(obj$par))] <- 1
Lwr[grep("log_sigmap", names(obj$par))] <- log(0)
Upr[grep("log_sigmap", names(obj$par))] <- log(2)
#Lwr[grep("log_cpue_pro", names(obj$par))] <- log(0)
#Upr[grep("log_cpue_pro", names(obj$par))] <- log(2)

control <- list(eval.max = 100000, iter.max = 100000)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)

check_estimability(obj = obj)

sd <- sdreport(obj)
pred_pars <- exp(opt$par)
names(obj$env)

# Plot Observed vs Predicted CPUE ----

cpue_df <- data.frame(Year = CPUE$Year, Season = CPUE$Season, Observed = CPUE$Median, Predicted = obj$report()$cpue_pred, Index = CPUE$Index) %>% 
  mutate(Index = factor(Index, levels = l_cpue)) 

q <- ggplot(cpue_df, aes(x = Year)) +
  geom_point(aes(y = Observed), color = "grey", size = 2) +
  geom_line(aes(y = Predicted, group = interaction(Index, Season)), 
            color = "red", linewidth = 1.2) +
  facet_grid(Index ~ Season, scales = "free_y") +
  labs(x = "Fishing year", y = "CPUE")
q

resid <- oneStepPredict(obj = obj, observation.name = "cpue_obs", method = "oneStepGeneric", trace = FALSE)$residual
resid_df <- data.frame(Time = data$cpue_time, Season = CPUE$Season, Residual = resid, Index = CPUE$Index) %>% 
  mutate(Index = factor(Index, levels = l_cpue))

p <- ggplot(resid_df, aes(x = Time, y = Residual)) +
  geom_hline(yintercept = 0, color = "#00BA38", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = c(-2, 2), color = "#F8766D", linetype = "dashed", linewidth = 0.5) +
  geom_linerange(aes(ymin = 0, ymax = Residual), color = "#619CFF", linewidth = 0.8) +
  geom_point(color = "#619CFF", size = 2) +
  facet_wrap(Index ~ Season, ncol =2) +
  labs(x = "Fishing year", y = "OSA residuals")
p

# Plot biomass ----

names(obj_full$report())

pred_Bfull <- data.frame(Year = 1979:2024, Season = 1, B = obj_full$report()$biomass_adj_yr[,1]*2.5, Model = "Full")
pred_Bbdm <- data.frame(Year = catch$Year, Season = catch$Season, B = obj$report()$Bpred, Model = "BDM") %>%
  filter(Season == 1) %>%
  filter(Year >= 1979)
pred_B <- bind_rows(pred_Bfull, pred_Bbdm) 

# why are both biomasses in different scales? 

pp <- ggplot(data = pred_B, aes(x = Year, y = B, color = Model)) +
  geom_point() +
  geom_line() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3500)) +
  labs(x = "Fishing year", y = "Biomass")
pp
# MCMC ----

library(adnuts)

# mcmc <- sample_sparse_tmb(obj = obj, metric = "auto", iter = 4000, warmup = 3000, chains = 4, cores = 4,
#                           control = list(adapt_delta = 0.995), init = "last.par.best", globals = list(evaluate_priors = evaluate_priors))
mcmc <- sample_snuts(obj = obj, metric = "auto", chains = 4, cores = 4, init = "last.par.best", globals = list(evaluate_priors = evaluate_priors))

pairs_rtmb(fit = mcmc, pars = 1:5, order = "slow")
pairs_rtmb(fit = mcmc, pars = 1:5, order = "mismatch")
# pairs_rtmb(fit = mcmc, pars = 1:5, order = "divergent")
plot_marginals(fit = mcmc, pars = 1:16)

# B_post <- get_posterior(object=obj, posterior=mcmc, pars = "log_B", iters = NULL, option = 2, type = "df")
posteriors <- extract_samples(mcmc, as.list=TRUE)

# is there a function already build in decamod to do this? 
# combine 4 chains
post_all <- do.call(rbind, posteriors)
post_all_exp <- exp(post_all)

stats <- t(apply(post_all_exp, 2, function(x) {
  c(
    mean   = mean(x),
    median = median(x),
    q05    = quantile(x, 0.05),
    q95    = quantile(x, 0.95)
  )
}))

stats_df <- data.frame(
  parameter = sub("^log_", "", rownames(stats)),
  stats
)
rownames(stats_df) <- NULL

head(stats_df)
B_post <- stats_df %>%
  filter(grepl("^B\\[", parameter))

B_post <- data.frame(Year = catch$Year[-1], Season = catch$Season[-1], B = B_post, Model = "MCMC") %>%
  filter(Season == 1) %>%
  filter(Year >= 1979)
head(B_post) # estimated B is too large, K is too large 
