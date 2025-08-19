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

#load("C:/Maite/CRA2/2025/models/runs/base.rda") # updated
load("base.rda")

## some outputs for the CRA 2 model have already been saved in data frames for plotting
## see B_t for vulnerable biomass by season, and B is adjusted vulnerable biomass including B0
names(output)

n_year <- data$n_year # base from 1979 to 2024
n_season <- data$n_season
n_region <- data$n_region
n_sex <- data$n_sex
n_fishery <- data$n_fishery
midpoint_l <- data$midpoint_l
n_iter <- max(output$Recruitment$Iteration)

#########################################
## model outputs for plotting projections later
## method for getting derived values from the posterior distribution
# outputs <- get_posterior(object = obj, posterior = mcmc,
#                          pars = c("recruitment_ytrsl", 
#                                   "biomass_adj_ytrsl", "B0", "cpue_pred",
#                                   "selectivity_ytrsfl", 
#                                   "U_ytrf"),
#                          option = 2, type = "list")
# save(outputs, file= "outputs.rda")

load("outputs.rda")

df_list <- post_to_plot(outputs)
n_iter <- max(df_list[[1]]$Iteration)

cpue_pred <- df_list$cpue_pred %>%
  mutate(Year_index = data$cpue_year[V1],
         Year = l_year[Year_index],
         Season = data$cpue_season[V1],
         q = data$cpue_q[V1],
         Obs = data$cpue_obs[V1],
         SD = data$cpue_sd[V1],
         Units = data$cpue_units[V1],
         Fishery = data$cpue_fishery[V1],
         Sex = data$cpue_sex[V1],
         Region = data$cpue_region[V1]) %>%
  select(-V1)

#lp <- get_posterior(object = obj, posterior = mcmc, pars = c("biomass_adj_ytrsl"), type = "list")
#lp_list <- post_to_plot(lp)

# vuln <- lp_list$biomass_adj_ytrsl %>%
#  group_by(Iteration, Year, Season, Region) %>%
#  summarise(value = sum(value))

###############################
## project dynamics
n_proj <- 30

###############################
## recruitment deviates
# max p etc. come from arima
# y1 y2 = years over which we'll sample, from string of recdevs in terms of index
# should use 1980-last estimated recdev (n_year-4 but check)
# number of years to project
# to use lognormal recruitment deviates, set arima = FALSE
rdev <- project_recruitment(obj = obj, mcmc = mcmc, y1 = 2, y2 = (n_year - 4), ny = n_proj, arima = FALSE)
Rdev_yr <- rdev[[1]]

dimnames(Rdev_yr) <- list(Iteration = 1:n_iter, Year = 1:n_proj, Region = l_region)
Rdev_df <- reshape2::melt(Rdev_yr)

ggplot(Rdev_df, aes(x = Year, y = value)) +
  stat_summary(fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(fun = "median", geom = "line") +
  # geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Recruitment deviate") +
  # scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(l_year), max(l_year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~Region, scales = "free_y")

##############################################
## prepare inputs with projection dimensions
# proj_inputs <- prepare_proj(data = data, obj = obj, adfit = mcmc, n_proj = n_proj)
# proj_inputs$Rdev_yr <- Rdev_yr
# save(proj_inputs, file="proj_inputs.rda")

load("proj_inputs.rda")
## 3 rules: one with the same catch as the last year, one with x2 of the catch of the last year and one with x3 of the catch of last year
catch_scenarios <- list(
  same   = proj_inputs$catch_ytrf,
  x2 = proj_inputs$catch_ytrf * 2,
  x3 = proj_inputs$catch_ytrf * 3
)

scenario_inputs <- lapply(catch_scenarios, function(catch_vec) {
  inputs <- proj_inputs
  inputs$catch_ytrf <- catch_vec
  inputs
})

# ###########################################################################################
# ## Helper to process metrics and checks for NULL
# process_metric <- function(array_data, dimnames_list,
#                            filter_expr = NULL, sum_by = NULL, iteration) {
#   if (is.null(array_data)) {
#     return(NULL)
#   }
#   dimnames(array_data) <- dimnames_list
#   df <- reshape2::melt(array_data)
#   if (!is.null(filter_expr)) {
#     df <- dplyr::filter(df, !!rlang::parse_expr(filter_expr))
#   }
#   if (!is.null(sum_by)) {
#     df <- df %>%
#       group_by(across(all_of(sum_by))) %>%
#       summarise(value = sum(value), .groups = "drop")
#   }
#   mutate(df, Iteration = iteration)
# }
# 
# ## Storage lists
# results_raw <- list() # raw outputs from do_dynamics
# results_df  <- list() # processed metrics
# 
# start <- Sys.time()
# for (scen in names(scenario_inputs)) {
# 
#   message("Running scenario: ", scen)
#   results_raw[[scen]] <- vector("list", n_iter)
#   results_df[[scen]] <- list(
#     catch = vector("list", n_iter),
#     u = vector("list", n_iter),
#     b = vector("list", n_iter),
#     #ssb  = vector("list", n_iter),
#     N = vector("list", n_iter),
#     rec = vector("list", n_iter),
#     bt  = vector("list", n_iter)
#   )
# 
#   for (i in seq_len(n_iter)) {
# 
#     # Run do_dynamics
#     run <- tryCatch({    # If an iteration fails, it logs a warning but the loop continues
#       with(scenario_inputs[[scen]], do_dynamics(
#         log_F = NULL,
#         logit_h = logit_h[i,1],
#         R0 = R0[i],
#         recruitment_size_sl = recruitment_size_sl,
#         Rdev_yr = array(Rdev_yr[i,,], dim = c(n_proj, n_region)),
#         Rsigma = Rsigma[i,1],
#         numbers_rsl = array(numbers_rsl[i,,,], dim = c(n_region, n_sex, length(midpoint_l))),
#         growth_ytrsll = array(growth_ytrsll[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_length, n_length)),
#         M_ytrsl = array(M_ytrsl[i,,,,,], dim = c(n_proj, n_season, n_region, n_sex, length(midpoint_l))),
#         catch_ytrf = catch_ytrf,
#         selectivity_ytrsfl = array(selectivity_ytrsfl[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_fishery, n_length)),
#         legal_ytrsfl = legal_ytrsfl,
#         retained_ytrsfl = retained_ytrsfl,
#         handling_mortality_y = handling_mortality_y,
#         weight_ytrsl = weight_ytrsl,
#         maturity_ytrsl = maturity_ytrsl,
#         catch_like = -1,
#         rsigma_bias = 0,
#         U_ytrf_in = NULL
#       ))
#     }, error = function(e) {
#       warning("Iteration ", i, " in scenario ", scen, " failed: ", e$message)
#       NULL
#     })
# 
#     results_raw[[scen]][[i]] <- run
#     if (is.null(run)) next  # skip processing if model failed
# 
#     ## Process metrics
#     results_df[[scen]]$catch[[i]] <- process_metric(
#       run$pred_catch_ytrsfl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj), Season = l_season, Region = l_region,
#                            Sex = l_sex, Fishery = l_fishery, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season", "Fishery"),
#       iteration = i
#     )
# 
#     results_df[[scen]]$u[[i]] <- process_metric(
#       run$U_ytrf,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj),
#                            Season = l_season, Region = l_region, Fishery = l_fishery),
#       sum_by = c("Region", "Year", "Season", "Fishery"),
#       iteration = i
#     )
# 
#     results_df[[scen]]$b[[i]] <- process_metric(
#       run$biomass_vuln_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
# 
#     results_df[[scen]]$N[[i]] <- process_metric(
#       run$numbers_cpue_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj),
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season", "Sex", "Bin"),
#       iteration = i
#     )
# 
#     results_df[[scen]]$rec[[i]] <- process_metric(
#       run$recruitment_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj),
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       # filter_expr = "Season == 'AW'",
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
# 
#     results_df[[scen]]$bt[[i]] <- process_metric(
#       run$biomass_vuln_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj),
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
#   }
# }
# end <- Sys.time()
# message("Total run time: ", round(end - start, 2), "minutes")
# 
# names(results_df[[1]])
# 
# # summaries to plot
# metrics <- names(results_df[[1]])
# 
# catch_proj_outputs <- list()
# 
# for (metric in metrics) {
#   # combine all scenarios and all iterations for this metric
#   catch_proj_outputs[[metric]] <- bind_rows(
#     lapply(names(results_df), function(scen) {
#       bind_rows(results_df[[scen]][[metric]]) %>%
#         mutate(Scenario = scen)  # add scenario info
#     })
#   )
# }
# save(catch_proj_outputs, file = "catch_proj_outputs.rda")
###########################################################################################
load("catch_proj_outputs.rda")
metrics <- c("catch", "u", "b", "N", "rec", "bt")  
combined_df <- bind_rows(
  lapply(metrics, function(m) {
    df <- catch_proj_outputs[[m]]
    df$metric <- m  
    df
  })
)
head(combined_df)
tail(combined_df)

############################
## PLOTS
# Observed catch
mcatch <- data$catch_ytrf
dimnames(mcatch) <- list(Year = l_year, Season = l_season, Region = l_region, Fishery = l_fishery)
mcatch <- reshape2::melt(mcatch) %>%
  mutate(Rule = "Assessment")  # observed
head(mcatch)

# Projected catch
pred_catch <- combined_df %>%
  filter(metric == "catch") %>%
  mutate(Rule = Scenario)  

# Combine observed + all projections
df_plot <- bind_rows(mcatch %>% select(Region, Year, Season, Fishery, value, Rule),
  pred_catch %>% 
    filter(Iteration==1) %>% 
    select(Region, Year, Season, Fishery, value, Rule)
)

# --- Plot ---
p <- ggplot(df_plot, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), 
               fun.min = function(x) quantile(x, 0.05), 
               fun.max = function(x) quantile(x, 0.95), 
               geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(xintercept = max(l_year), linetype = 2) +
  xlab("Fishing year") + ylab("Catch (tonnes)") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_plot$Year), max(df_plot$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")
p
ggsave("projected_catch_cte_catch.pdf",width = 10, height = 6)
### plot adjusted vulnerable biomass
names(df_list)
df_B <- df_list$biomass_adj_ytrsl %>%
  group_by(Iteration, Year, Season, Region) %>%
  summarise(value = sum(value), .groups = "drop") %>% mutate(Rule = "Assessment")

df_B <- bind_rows(df_B, 
                    combined_df %>%
                      filter(metric == "b") %>%  
                      select(Region, Year, Season, value, Iteration, Scenario) %>%
                      rename(Rule = Scenario)) 
tail(df_B)

p <- ggplot(df_B, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Adjusted vulnerable biomass (tonnes)") +
  scale_y_continuous(limits = c(0, 4000), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_B$Year), max(df_B$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~ Season, scales = "free_y")
p
ggsave("adj_vul_biomass_proj_catch.pdf", width = 10, height = 6)

# ## projections: constant U
## 3 rules: one with the same U as the last year, one with x2 of U of the last year and one with x3 of the U of last year

U_scenarios <- list(
  same = proj_inputs$U_ytrf,
  x2 = proj_inputs$U_ytrf * 2,
  x3 = proj_inputs$U_ytrf * 3
)

scenario_inputs <- lapply(U_scenarios, function(U_vec) {
  inputs <- proj_inputs
  inputs$U_ytrf <- U_vec
  inputs
})

###########################################################################################
## Storage lists
# results_Uraw <- list() # raw outputs from do_dynamics
# results_Udf  <- list() # processed metrics
# 
# start <- Sys.time()
# for (scen in names(scenario_inputs)) {
#   
#   message("Running scenario: ", scen)
#   results_Uraw[[scen]] <- vector("list", n_iter)
#   results_Udf[[scen]] <- list(
#     catch = vector("list", n_iter),
#     u = vector("list", n_iter),
#     b = vector("list", n_iter),
#     #ssb  = vector("list", n_iter),
#     N = vector("list", n_iter),
#     rec = vector("list", n_iter),
#     bt  = vector("list", n_iter)
#   )
#   
#   for (i in seq_len(n_iter)) {
#     
#     # Run do_dynamics
#     run <- tryCatch({    # If an iteration fails, it logs a warning but the loop continues
#       with(scenario_inputs[[scen]], do_dynamics(
#         log_F = NULL,
#         logit_h = logit_h[i,1],
#         R0 = R0[i],
#         recruitment_size_sl = recruitment_size_sl,
#         Rdev_yr = array(Rdev_yr[i,,], dim = c(n_proj, n_region)), 
#         Rsigma = Rsigma[i,1], 
#         numbers_rsl = array(numbers_rsl[i,,,], dim = c(n_region, n_sex, length(midpoint_l))), 
#         growth_ytrsll = array(growth_ytrsll[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_length, n_length)),
#         M_ytrsl = array(M_ytrsl[i,,,,,], dim = c(n_proj, n_season, n_region, n_sex, length(midpoint_l))), 
#         catch_ytrf = catch_ytrf, 
#         selectivity_ytrsfl = array(selectivity_ytrsfl[i,,,,,,], dim = c(n_proj, n_season, n_region, n_sex, n_fishery, n_length)), 
#         legal_ytrsfl = legal_ytrsfl,
#         retained_ytrsfl = retained_ytrsfl, 
#         handling_mortality_y = handling_mortality_y,
#         weight_ytrsl = weight_ytrsl, 
#         maturity_ytrsl = maturity_ytrsl, 
#         catch_like = -1, 
#         rsigma_bias = 0, 
#         U_ytrf_in = array(U_ytrf[i,,,,], dim = c(n_proj, n_season, n_region, n_fishery))))
#     }, error = function(e) {
#       warning("Iteration ", i, " in scenario ", scen, " failed: ", e$message)
#       NULL
#     })
#     
#     results_Uraw[[scen]][[i]] <- run
#     if (is.null(run)) next  # skip processing if model failed
#     
#     ## Process metrics
#     results_Udf[[scen]]$catch[[i]] <- process_metric(
#       run$pred_catch_ytrsfl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj), Season = l_season, Region = l_region, 
#                            Sex = l_sex, Fishery = l_fishery, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season", "Fishery"),
#       iteration = i
#     )
#     
#     results_Udf[[scen]]$u[[i]] <- process_metric(
#       run$U_ytrf,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj), 
#                            Season = l_season, Region = l_region, Fishery = l_fishery),
#       sum_by = c("Region", "Year", "Season", "Fishery"),
#       iteration = i
#     )
#     
#     results_Udf[[scen]]$b[[i]] <- process_metric(
#       run$biomass_vuln_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, by = 1, length.out = n_proj), Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
#     
#     results_Udf[[scen]]$N[[i]] <- process_metric(
#       run$numbers_cpue_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj),
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season", "Sex", "Bin"),
#       iteration = i
#     )
#     
#     results_Udf[[scen]]$rec[[i]] <- process_metric(
#       run$recruitment_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj), 
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       # filter_expr = "Season == 'AW'",
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
#     
#     results_Udf[[scen]]$bt[[i]] <- process_metric(
#       run$biomass_vuln_ytrsl,
#       dimnames_list = list(Year = seq(max(l_year)+1, length.out = n_proj), 
#                            Season = l_season, Region = l_region, Sex = l_sex, Bin = midpoint_l),
#       sum_by = c("Region", "Year", "Season"),
#       iteration = i
#     )
#   }
# }
# end <- Sys.time()
# message("Total run time: ", round(end - start, 2), "minutes")
# 
# names(results_Udf[[1]])
# 
# # summaries to plot
# metrics <- names(results_Udf[[1]])
# 
# U_proj_outputs <- list()
# 
# for (metric in metrics) {
#   # combine all scenarios and all iterations for this metric
#   U_proj_outputs[[metric]] <- bind_rows(
#     lapply(names(results_Udf), function(scen) {
#       bind_rows(results_Udf[[scen]][[metric]]) %>%
#         mutate(Scenario = scen)  # add scenario info
#     })
#   )
# }
# 
# save(U_proj_outputs, file = "U_proj_outputs.rda")
###########################################################################################
load("U_proj_outputs.rda")

combined_Udf <- bind_rows(
  lapply(metrics, function(m) {
    df <- U_proj_outputs[[m]]
    df$metric <- m  
    df
  })
)
head(combined_Udf)
tail(combined_Udf)

############################
## PLOTS
# Projected catch
pred_Ucatch <- combined_Udf %>%
  filter(metric == "catch") %>%
  mutate(Rule = Scenario)  

# Combine observed + all projections
df_Uplot <- bind_rows(mcatch %>% select(Region, Year, Season, Fishery, value, Rule),
                     pred_Ucatch %>% 
                       select(Region, Year, Season, Fishery, value, Rule)
)

# --- Plot ---
p <- ggplot(df_Uplot, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), 
               fun.min = function(x) quantile(x, 0.05), 
               fun.max = function(x) quantile(x, 0.95), 
               geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(xintercept = max(l_year), linetype = 2) +
  xlab("Fishing year") + ylab("Catch (tonnes)") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_plot$Year), max(df_plot$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")
p
ggsave("projected_catch_cte_U.pdf",width = 10, height = 6)

# Projected U
pred_U <- combined_Udf %>%
  filter(metric == "u") %>%
  mutate(Rule = Scenario)  

# Combine estimated + all projections
head(df_list$U_ytrf)
pred_U_plot <- bind_rows(df_list$U_ytrf %>% mutate(Rule="Assessment"),
          pred_U %>% 
          select(Iteration, Region, Year, Season, Fishery, value, Rule)
)
tail(pred_U_plot)
# --- Plot ---
p <- ggplot(pred_U_plot, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), 
               fun.min = function(x) quantile(x, 0.05), 
               fun.max = function(x) quantile(x, 0.95), 
               geom = "ribbon", alpha = 0.25, color = NA) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(xintercept = max(l_year), linetype = 2) +
  xlab("Fishing year") + ylab("U") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_plot$Year), max(df_plot$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")
p
ggsave("projected_U_cte.pdf",width = 10, height = 6)

### plot adjusted vulnerable biomass
df_UB <- bind_rows(df_B, 
                  combined_Udf %>%
                    filter(metric == "b") %>%  
                    select(Region, Year, Season, value, Iteration, Scenario) %>%
                    rename(Rule = Scenario)) 
tail(df_UB)

p <- ggplot(df_UB, aes(x = Year, y = value)) +
  stat_summary(aes(fill = Rule), fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95), geom = "ribbon", alpha = 0.25) +
  stat_summary(aes(color = Rule), fun = "median", geom = "line") +
  geom_vline(aes(xintercept = max(l_year)), linetype = 2) +
  xlab("Fishing year") + ylab("Adjusted vulnerable biomass (tonnes)") +
  scale_y_continuous(limits = c(0, 4000), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_B$Year), max(df_B$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~ Season, scales = "free_y")
p
ggsave("adj_vul_biomass_proj_U.pdf", width = 10, height = 6)

############## generate inputs for the Estimation Model (EM)
names(obj$simulate())
historical <- data.frame(Year = data$cpue_year, CPUE_obs = data$cpue_obs, CPUE_sd = data$cpue_sd, CPUE_sim = obj$simulate()$cpue_obs, 
                        CPUE_sigma = obj$simulate()$cpue_sigma, q= data$cpue_q)
historical <- historical %>%
  mutate(Year = 1978 + Year) %>%
  mutate(Season = rep(c("AW", "SS"), length.out = n()))

catch_tmp <- mcatch %>% 
  group_by(Year, Season)%>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
historical <- historical %>%
  left_join(catch_tmp %>% rename(Catch = value), by = c("Year", "Season"))

head(historical)

#projected catch from U scenarios 
projected <- pred_catch %>%
  group_by(Region, Year, Season, Iteration, Rule) %>%
  summarise(Catch = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(Scenario = Rule)
head(projected)

# calculate projected CPUE
# extract q3 for logbook CPUE and pow 2 
post <- adnuts::extract_samples(fit = mcmc) %>%
  mutate(Iteration = 1:nrow(.)) %>%
  pivot_longer(-Iteration, names_to = "par", values_to = "val")
# I need for Logbook CPUE: "log_q[3]" and "log_pow[2]" 
params <- post %>%
  filter(par %in% c("log_q[3]", "log_pow[2]")) %>%
  mutate(val = exp(val),                    
         par = recode(par, 
                      "log_q[3]" = "q_3",
                      "log_pow[2]" = "pow_2")) %>%
  pivot_wider(names_from = par, values_from = val)
head(params)

# CPUE ---- units=1 for logbook CPUE (proj_inputs$cpue_units)
head(U_proj_outputs$N) # this is numbers_cpue_ytrsl[iy, it, ir, ix,] 
# Legal is legal_ytrsfl
dim(proj_inputs$legal_ytrsfl)
legal <- expand.grid(
  Year = seq(max(l_year) + 1, length.out = n_proj),
  Season = l_season,
  Region = l_region,
  Sex = l_sex,
  Fishery = l_fishery,
  Bin = midpoint_l
)
# array into vector
values <- as.vector(proj_inputs$legal_ytrsfl)
# data frame
legal$value <- values
tail(legal)

# Sel is selectivity_ytrsfl
dim(proj_inputs$selectivity_ytrsfl)
Sel <- expand.grid(
  Iteration = c(1:n_iter),
  Year = seq(max(l_year) + 1, length.out = n_proj),
  Season = l_season,
  Region = l_region,
  Sex = l_sex,
  Fishery = l_fishery,
  Bin = midpoint_l
)
# array into vector
values <- as.vector(proj_inputs$selectivity_ytrsfl)
# data frame
Sel$value <- values
tail(Sel)

N <- U_proj_outputs$N
# I need N duplicated for NsL and SL fisheries 
fisheries <- unique(Sel$Fishery)
scenarios <- unique(N$Scenario)
N_expanded <- N %>% 
  crossing(Fishery = fisheries)
S_expanded <- Sel %>%  # I need it triplicated to multiply to each scenario of U
  crossing(Scenario = scenarios)

tmp <- N_expanded %>%
  inner_join(S_expanded,
             by = c("Iteration","Year","Season","Region","Sex","Bin","Fishery","Scenario"),
             suffix = c("_N","_Sel"))
head(tmp)
# here  cpue_tmp <- numbers_cpue_ytrsl[iy, it, ir, ix,] * selectivity_ytrsfl[iy, it, ir, ix, ij,]
tmp <- tmp %>%
  mutate(cpue_tmp = value_N * value_Sel)  # this is 

head(tmp)

# here  if (cpue_units[i] == 1) { cpue_tmp <- cpue_tmp * legal_ytrsfl[iy, it, ir, ix, ij,]
cpue_tmp <- tmp %>%
  left_join(legal, 
            by = c("Year", "Season", "Region", "Sex", "Fishery", "Bin")) %>%
  mutate(cpue_tmp2 = cpue_tmp * value)
head(cpue_tmp)

# here sum(cpue_tmp) 
cpue_tmp_sum <- cpue_tmp %>%
  group_by(Region, Year, Season, Iteration, Scenario) %>%
  summarise(cpue_tmp = sum(cpue_tmp2, na.rm = TRUE), .groups = "drop")
head(cpue_tmp_sum)

# here cpue_pred[i] <- q_yi[iy, iq] * sum(cpue_tmp)^exp(log_pow[iq]) + 1e-6
cpue_proj <- cpue_tmp_sum %>%
  left_join(params, by = "Iteration") %>%
  mutate(cpue_pred = q_3 * (cpue_tmp ^ pow_2) + 1e-6) 
head(cpue_proj)

projected <- projected %>%
  left_join(cpue_proj %>% select(Region, Year, Season, Iteration, Scenario, cpue_pred),
    by = c("Region", "Year", "Season", "Iteration", "Scenario")
  )
head(projected)

## Add to the data frame of adjusted vulnerable biomass (by season), B0 to compare outputs in the EM phase 
head(df_UB)
names(df_list)
head(df_list$B0)

df_UB <- df_UB %>%
  left_join(df_list$B0 %>% 
      select(Iteration, B0 = value),by = "Iteration")
head(df_UB)

head(projected)
tail(historical)
# PLOT CPUE
ggplot() +
  # Historical CPUE
  geom_line(data = historical %>% filter(q==3),aes(x = Year, y = CPUE_obs), col="red") +
  geom_point(data = historical %>% filter(q==3),aes(x = Year, y = CPUE_obs), col="red") +
  # Projected CPUE
  stat_summary(data = projected,aes(x = Year, y = cpue_pred, fill = Scenario),
    fun.min = function(x) quantile(x, 0.05),fun.max = function(x) quantile(x, 0.95),
    geom = "ribbon",alpha = 0.25) +
  stat_summary(data = projected,aes(x = Year, y = cpue_pred, color = Scenario),
    fun = median, geom = "line") +
  geom_vline(xintercept = max(historical$Year), linetype = 2) +
  xlab("Fishing year") + ylab("Logbook CPUE") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ Season, scales = "free_y")

ggsave("obs&proj_CPUE.pdf", width = 10, height = 6)
