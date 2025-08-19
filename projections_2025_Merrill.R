library(decamod)
library(TMBhelper)
library(tidyverse)
library(kableExtra)
library(adnuts)
library(parallel)
library(doParallel)
library(abind)
library(forecast)
options(dplyr.summarise.inform = FALSE)

theme_set(theme_bw())

load("base.rda")
#load("runs/base_msy.rda")

n_year <- data$n_year
n_season <- data$n_season
n_region <- data$n_region
n_sex <- data$n_sex
n_fishery <- data$n_fishery
midpoint_l <- data$midpoint_l
n_iter <- max(output$Recruitment$Iteration)

#########################################
## model outputs for plotting projections
outputs <- get_posterior(object = obj, posterior = mcmc,
                         pars = c("recruitment_ytrsl", "Rdev_yr",
                                  "biomass_adj_ytrsl", "B0",
                                  "biomass_ssb_yr", "SSB0",
                                  "U_ytrf"),
                         option = 2, type = "list")

df_list <- post_to_plot(outputs)
n_iter <- max(df_list[[1]]$Iteration)

###############################
## project dynamics
n_proj <- 5

###############################
## recruitment deviates
# max p etc. come from arima keep as is
# y1 y2 = years over which we'll sample, from string of recdevs in terms of index
## should use 1980-last estimated recdev (n_year-4 but check)
## number of years to project
rdev <- project_recruitment(obj = obj, mcmc = mcmc, y1 = 2, y2 = (n_year - 4), ny = n_proj, arima = TRUE)
Rdev_yr <- rdev[[1]]

post <- adnuts::extract_samples(fit = mcmc) %>%
  mutate(Iteration = 1:nrow(.)) %>%
  pivot_longer(-Iteration, names_to = "par", values_to = "val")


post$par2 <- sapply(1:nrow(post), function(x) {
  if (grepl('\\[', post$par[x])) {
    out <- strsplit(post$par[x], "\\[")[[1]][1]
  } else {
    out <- post$par[x]
  }
  return(out)
})

post$par_num <- sapply(1:nrow(post), function(x) {
  if (grepl('\\[', post$par[x])) {
    out <- as.numeric(strsplit(strsplit(post$par[x], "\\[")[[1]][2], "\\]")[[1]])
  } else {
    out <- 1
  }
  return(out)
})


Rdev_yr_model <- post %>%
  filter(par2 == "Rdev_yr") %>%
  mutate(Year = l_year[par_num]) %>%
  select(Iteration, Year, val) %>%
  rename(value = val) %>%
  mutate(Region = l_region)

dimnames(Rdev_yr) <- list(Iteration = 1:n_iter, Year = max(l_year) + (1:n_proj), Region = l_region)
Rdev_df <- reshape2::melt(Rdev_yr)
Rdev_all <- bind_rows(Rdev_yr_model, Rdev_df)
ggplot(Rdev_all, aes(x = Year, y = value)) +
  stat_summary(fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Recruitment deviate") +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(l_year), max(l_year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~Region, scales = "free_y")

##############################################
## prepare inputs with projection dimensions
proj_inputs <- prepare_proj(data = data, obj = obj, adfit = mcmc, n_proj = n_proj)
proj_inputs$Rdev_yr <- Rdev_yr

################################
## projection: constant catch
proj_catch <- list()
catch_pcatch <- list()
u_pcatch <- list()
b_pcatch <- list()
ssb_pcatch <- list()
rec_pcatch <- list()
bt_pcatch <- list()

start <- Sys.time()
for (i in 1:n_iter) {
  proj_catch[[i]] <- with(proj_inputs, do_dynamics(log_F = NULL, 
                                                   logit_h = logit_h[i,1], 
                                                   R0 = R0[i],
                                                   recruitment_size_sl = recruitment_size_sl,
                                                   Rdev_yr = array(Rdev_yr[i,,], dim = c(n_proj, n_region)), 
                                                   Rsigma = Rsigma[i,1], 
                                                   numbers_rsl = array(numbers_rsl[i,,,], dim = c(n_region, n_sex, length(midpoint_l))), 
                                                   growth_ytrsll = array(growth_ytrsll[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_length, n_length)),
                                                   M_ytrsl = array(M_ytrsl[i,,,,,], dim = c(n_proj, n_season, n_region, n_sex, length(midpoint_l))), 
                                                   catch_ytrf = catch_ytrf, 
                                                   selectivity_ytrsfl = array(selectivity_ytrsfl[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_fishery, n_length)), 
                                                   legal_ytrsfl = legal_ytrsfl,
                                                   retained_ytrsfl = retained_ytrsfl, 
                                                   handling_mortality_y = handling_mortality_y,
                                                   weight_ytrsl = weight_ytrsl, 
                                                   maturity_ytrsl = maturity_ytrsl, 
                                                   catch_like = -1, 
                                                   rsigma_bias = 0, 
                                                   U_ytrf_in = NULL))
  
  catch <- proj_catch[[i]]$pred_catch_ytrsfl
  dimnames(catch) <- list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Fishery = l_fishery, Bin = midpoint_l)
  catch <- reshape2::melt(catch) %>%
    group_by(Region, Year, Season, Fishery) %>%
    summarise(value = sum(value)) %>%
    mutate(Iteration = i)
  catch_pcatch[[i]] <- catch
  
  u <- proj_catch[[i]]$U_ytrf
  dimnames(u) <- list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Fishery = l_fishery)
  u <- reshape2::melt(u) %>%
    group_by(Region, Year, Season, Fishery) %>%
    mutate(Iteration = i)
  u_pcatch[[i]] <- u
  
  b <- proj_catch[[i]]$biomass_adj_ytrsl
  dimnames(b) <- list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Region = l_region)
  b <- reshape2::melt(b) %>%
    mutate(Iteration = i)
  b_pcatch[[i]] <- b
  
  ssb <- proj_catch[[i]]$biomass_ssb_ytrsl
  dimnames(ssb) <- list(Year =seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l)
  ssb <- reshape2::melt(ssb) %>%
    filter(Season == "AW") %>%
    group_by(Region, Year) %>%
    summarise(value = sum(value)) %>%
    mutate(Iteration = i)
  ssb_pcatch[[i]] <- ssb
  
  rec <- proj_catch[[i]]$recruitment_ytrsl
  dimnames(rec) <- list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l)
  rec <- reshape2::melt(rec) %>%
    filter(Season == "AW") %>%
    group_by(Region, Year) %>%
    summarise(value = sum(value)) %>%
    mutate(Iteration = i)
  rec_pcatch[[i]] <- rec
  
  bt <- proj_catch[[i]]$biomass_vuln_ytrsl
  dimnames(bt) <- list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l)
  bt <- reshape2::melt(bt) %>%
    group_by(Region, Year, Season) %>%
    summarise(value = sum(value)) %>%
    mutate(Iteration = i)
  bt_pcatch[[i]] <- bt
}
end <- Sys.time() - start

catch_pcatch2 <- bind_rows(catch_pcatch)
u_pcatch2 <- bind_rows(u_pcatch)
b_pcatch2 <- bind_rows(b_pcatch)
ssb_pcatch2 <- bind_rows(ssb_pcatch)
rec_pcatch2 <- bind_rows(rec_pcatch)
bt_pcatch2 <- bind_rows(bt_pcatch)

catch_proj <- list()
catch_proj$Catch <- catch_pcatch2
catch_proj$U <- u_pcatch2
catch_proj$B <- b_pcatch2
catch_proj$SSB <- ssb_pcatch2
catch_proj$Rdev <- Rdev_df
catch_proj$Recruitment <- rec_pcatch2
catch_proj$B_t <- bt_pcatch2

#save(catch_proj, file = "runs/base_project5.rda")

# ################################
# ## projection: constant U
# proj_u <- list()
# catch_pu <- list()
# u_pu <- list()
# b_pu <- list()
# ssb_pu <- list()
# 
# start <- Sys.time()
# for (i in 1:n_iter) {
#   proj_u[[i]] <- with(proj_inputs, do_dynamics(log_F = NULL, 
#                                                    logit_h = logit_h[i,1], 
#                                                    R0 = R0[i],
#                                                    recruitment_size_sl = recruitment_size_sl,
#                                                    Rdev_yr = array(Rdev_yr[i,,], dim = c(n_proj, n_region)), 
#                                                    Rsigma = Rsigma[i,1], 
#                                                    numbers_rsl = array(numbers_rsl[i,,,], dim = c(n_region, n_sex, length(midpoint_l))), 
#                                                    growth_ytrsll = array(growth_ytrsll[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_length, n_length)),
#                                                    M_ytrsl = array(M_ytrsl[i,,,,,], dim = c(n_proj, n_season, n_region, n_sex, length(midpoint_l))), 
#                                                    catch_ytrf = catch_ytrf, 
#                                                    selectivity_ytrsfl = array(selectivity_ytrsfl[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_fishery, n_length)), 
#                                                    legal_ytrsfl = legal_ytrsfl,
#                                                    retained_ytrsfl = retained_ytrsfl, 
#                                                    handling_mortality_y = handling_mortality_y,
#                                                    weight_ytrsl = weight_ytrsl, 
#                                                    maturity_ytrsl = maturity_ytrsl, 
#                                                    catch_like = -1, 
#                                                    rsigma_bias = 0, 
#                                                    U_ytrf_in = array(U_ytrf[i,,,,], dim = c(n_proj, n_season, n_region, n_fishery))))
#   
#   catch <- proj_u[[i]]$pred_catch_ytrsfl
#   dimnames(catch) <- list(Year = 1:n_proj, Season = l_season, Region = l_region, Sex = l_sex, Fishery = l_fishery, Bin = midpoint_l)
#   catch <- reshape2::melt(catch) %>%
#     group_by(Region, Year, Season, Fishery) %>%
#     summarise(value = sum(value)) %>%
#     mutate(Iteration = i)
#   catch_pu[[i]] <- catch
#   
#   u <- proj_u[[i]]$U_ytrf
#   dimnames(u) <- list(Year = 1:n_proj, Season = l_season, Region = l_region, Fishery = l_fishery)
#   u <- reshape2::melt(u) %>%
#     group_by(Region, Year, Season, Fishery) %>%
#     mutate(Iteration = i)
#   u_pu[[i]] <- u
#   
#   b <- proj_u[[i]]$biomass_adj_ytrsl
#   dimnames(b) <- list(Year = 1:n_proj, Region = l_region)
#   b <- reshape2::melt(b) %>%
#     mutate(Iteration = i)
#   b_pu[[i]] <- b
#   
#   ssb <- proj_u[[i]]$biomass_ssb_ytrsl
#   dimnames(ssb) <- list(Year = 1:n_proj, Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l)
#   ssb <- reshape2::melt(ssb) %>%
#     filter(Season == "AW") %>%
#     group_by(Region, Year) %>%
#     summarise(value = sum(value)) %>%
#     mutate(Iteration = i)
#   ssb_pu[[i]] <- ssb
# }
# end <- Sys.time() - start
# 
# catch_pu2 <- bind_rows(catch_pu)
# u_pu2 <- bind_rows(u_pu)
# b_pu2 <- bind_rows(b_pu)
# ssb_pu2 <- bind_rows(ssb_pu)

############################
## compare projections
mcatch <- data$catch_ytrf
dimnames(mcatch) <- list(Year = l_year, Season = l_season, Region = l_region, Fishery = l_fishery)
mcatch <- reshape2::melt(mcatch)

catch_pcatch2 <- catch_pcatch2 %>%
  bind_rows(mcatch %>% filter(Year == max(Year))) %>%
  mutate(Rule = "Constant Catch")
# catch_pu2 <- catch_pu2 %>%
#   bind_rows(mcatch %>% filter(Year == max(Year))) %>%
#   mutate(Rule = "Constant U")
df <- bind_rows(mcatch %>% mutate(Rule = "Assessment"), catch_pcatch2) #, catch_pu2)
p <- ggplot(df, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Catch (tonnes)") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df$Year), max(df$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")

u_pcatch2 <- u_pcatch2 %>%
  bind_rows(df_list$U_ytrf %>% filter(Year == max(Year))) %>%
  mutate(Rule = "Constant Catch")
# u_pu2 <- u_pu2 %>%
#   bind_rows(dfm$U_ytrf %>% filter(Year == max(Year))) %>%
#   mutate(Rule = "Constant U")
df <- bind_rows(df_list$U_ytrf %>% mutate(Rule = "Assessment"), u_pcatch2)#, u_pu2)
p <- ggplot(df, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Exploitation rate") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df$Year), max(df$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")

b_pcatch2 <- b_pcatch2 %>%
  bind_rows(df_list$biomass_adj_ytrsl %>% filter(Year == max(Year))) %>%
  mutate(Rule = "Constant Catch")
# b_pu2 <- b_pu2 %>%
#   mutate(Year = max(l_year) + Year) %>%
#   bind_rows(dfm$biomass_adj_ytrsl %>% filter(Year == max(Year))) %>%
#   mutate(Rule = "Constant U")
df <- bind_rows(df_list$biomass_adj_ytrsl %>% mutate(Rule = "Assessment"), b_pcatch2)#, b_pu2)
p <- ggplot(df, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("AW adjusted vulnerable biomass (tonnes)") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df$Year), max(df$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~Region, scales = "free_y")

ssb_pcatch2 <- ssb_pcatch2 %>%
  bind_rows(df_list$biomass_ssb_yr %>% filter(Year == max(Year))) %>%
  mutate(Rule = "Constant Catch")
# ssb_pu2 <- ssb_pu2 %>%
#   mutate(Year = max(l_year) + Year) %>%
#   bind_rows(dfm$biomass_ssb_yr %>% filter(Year == max(Year))) %>%
#   mutate(Rule = "Constant U")
ssb0 <- df_list$SSB0 %>%
  rename(SSB0 = value)
df <- bind_rows(df_list$biomass_ssb_yr %>% left_join(ssb0) %>% mutate(Rule = "Assessment"), 
                ssb_pcatch2 %>% left_join(ssb0)) %>%
  mutate(RelSSB = value / SSB0)
p <- ggplot(df, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Spawning stock biomass (tonnes)") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df$Year), max(df$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~Region, scales = "free_y")

p <- ggplot(df, aes(x = Year, y = RelSSB)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Relative spawning stock biomass (tonnes)") +
  geom_hline(aes(yintercept = 0.2), color = "orange") +
  geom_hline(aes(yintercept = 0.1), color = "red") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df$Year), max(df$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~Region, scales = "free_y")
