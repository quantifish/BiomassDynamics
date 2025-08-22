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

load("C:/Maite/CRA2/2025/models/runs/base.rda") # updated
#load("base.rda")

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
# model outputs for plotting projections later
# method for getting derived values from the posterior distribution
# outputs <- get_posterior(object = obj, posterior = mcmc,
#                          pars = c("recruitment_ytrsl",
#                                   "biomass_adj_ytrsl", "B0", "cpue_pred",
#                                   "selectivity_ytrsfl",
#                                   "U_ytrf"),
#                          option = 2, type = "list")
# save(outputs, file= "outputs.rda")

load("outputs.rda")

df_list <- post_to_plot(outputs)
#n_iter <- max(df_list[[1]]$Iteration)

# cpue_pred <- df_list$cpue_pred %>%
#   mutate(Year_index = data$cpue_year[V1],
#          Year = l_year[Year_index],
#          Season = data$cpue_season[V1],
#          q = data$cpue_q[V1],
#          Obs = data$cpue_obs[V1],
#          SD = data$cpue_sd[V1],
#          Units = data$cpue_units[V1],
#          Fishery = data$cpue_fishery[V1],
#          Sex = data$cpue_sex[V1],
#          Region = data$cpue_region[V1]) %>%
#   select(-V1)
# head(cpue_pred)
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
# save(proj_inputs, file="proj_inputs.rda")

load("proj_inputs.rda")
## 3 rules: one with the same catch as the last year, one with x2 of the catch of the last year and one with x3 of the catch of last year
# catch_scenarios <- list(
#   same   = proj_inputs$catch_ytrf,
#   x2 = proj_inputs$catch_ytrf * 2,
#   x3 = proj_inputs$catch_ytrf * 3
# )
# 
# scenario_inputs <- lapply(catch_scenarios, function(catch_vec) {
#   inputs <- proj_inputs
#   inputs$catch_ytrf <- catch_vec
#   inputs
# })


# ## Helper to process metrics and checks for NULL
process_metric <- function(array_data, dimnames_list,
                           filter_expr = NULL, sum_by = NULL, iteration) {
  if (is.null(array_data)) {
    return(NULL)
  }
  dimnames(array_data) <- dimnames_list
  df <- reshape2::melt(array_data)
  if (!is.null(filter_expr)) {
    df <- dplyr::filter(df, !!rlang::parse_expr(filter_expr))
  }
  if (!is.null(sum_by)) {
    df <- df %>%
      group_by(across(all_of(sum_by))) %>%
      summarise(value = sum(value), .groups = "drop")
  }
  mutate(df, Iteration = iteration)
}
############################################################################################
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
# ggsave("projected_catch_cte_catch.pdf",width = 10, height = 6)
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
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_B$Year), max(df_B$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~ Season, scales = "free_y")
p
#ggsave("adj_vul_biomass_proj_catch.pdf", width = 10, height = 6)

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

#identical(proj_inputs, scenario_inputs[["same"]])
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

#save(U_proj_outputs, file = "U_proj_outputs.rda")
#save(results_Uraw, file = "results_U_arrays.rda")
###########################################################################################
load("U_proj_outputs.rda")
load("results_U_arrays.rda")

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
  scale_x_continuous(breaks = pretty(c(min(df_Uplot$Year), max(df_Uplot$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")
p
#ggsave("projected_catch_cte_U.pdf",width = 10, height = 6)

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
  scale_x_continuous(breaks = pretty(c(min(df_Uplot$Year), max(df_Uplot$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(Region + Season ~ Fishery, scales = "free_y")
p
#ggsave("projected_U_cte.pdf",width = 10, height = 6)

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
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks = pretty(c(min(df_B$Year), max(df_B$Year))), minor_breaks = seq(0, 1e6, 1), expand = c(0, 1)) +
  facet_wrap(~ Season, scales = "free_y")
p
# ggsave("adj_vul_biomass_proj_U.pdf", width = 10, height = 6)

############## generate inputs for the Estimation Model (EM)
input_catch  <- mcatch %>% 
  mutate(Iteration = 1) %>% rename(Catch=value)
head(input_catch)
#projected catch from U scenarios 
projected_catch <- pred_Ucatch %>%
  group_by(Region, Year, Fishery, Season, Iteration, Rule) %>%
  summarise(Catch = sum(value, na.rm = TRUE), .groups = "drop") 
head(projected_catch)

catch_input <- bind_rows(input_catch, 
                   projected_catch) 
                    
save(catch_input, file = "Input_catch.rda")

# calculate projected CPUE
# extract q3 for logbook CPUE and pow 2 
post <- adnuts::extract_samples(fit = mcmc) %>%
  mutate(Iteration = 1:nrow(.)) %>%
  pivot_longer(-Iteration, names_to = "par", values_to = "val")
# I need for Logbook CPUE: "log_q[3]" and "log_pow[2]" 
params <- post %>%
  filter(par %in% c("log_q[3]", "log_pow[2]")) %>%
  mutate(val = val,
         par = recode(par, 
                      "log_q[3]" = "log_q",
                      "log_pow[2]" = "log_pow")) %>%
  pivot_wider(names_from = par, values_from = val)
head(params)  # in log scale 

source("cpue_proj_fun.R")

# Scenario same U
cpue_proj_same <- list()
for(ii in 1:n_iter) {
    
    log_q_3 <- params$log_q[ii]
    log_pow_2 <- params$log_pow[ii]
    
    selectivity_ytrsfl <- proj_inputs$selectivity_ytrsfl[ii, , , , , , , drop = FALSE]
    selectivity_ytrsfl <- array(selectivity_ytrsfl, dim = dim(selectivity_ytrsfl)[-1])
    dim(selectivity_ytrsfl)
    
    numbers_cpue_ytrsl <- results_Uraw$same[[ii]]$numbers_cpue_ytrsl
    dim(numbers_cpue_ytrsl) # it does not have fishery 
    
    # function to calculate predicted projected CPUE
    tmp <- get_cpue_proj(data, parameters, log_q_3 = log_q_3, log_pow_2 = log_pow_2, n_proj = n_proj, numbers_cpue_ytrsl=numbers_cpue_ytrsl,
                    selectivity_ytrsfl = selectivity_ytrsfl) 
    # add iteration column
    tmp$Iteration <- ii
    
    # store 
    cpue_proj_same[[ii]] <- tmp
    
}
cpue_proj_same <- dplyr::bind_rows(cpue_proj_same, .id = "Iteration")
dim(cpue_proj_same)
head(cpue_proj_same)

# Scenario x2 U
cpue_proj_x2 <- list()
for(ii in 1:n_iter) {
  
  log_q_3 <- params$log_q[ii]
  log_pow_2 <- params$log_pow[ii]
  
  selectivity_ytrsfl <- proj_inputs$selectivity_ytrsfl[ii, , , , , , , drop = FALSE]
  selectivity_ytrsfl <- array(selectivity_ytrsfl, dim = dim(selectivity_ytrsfl)[-1])
  dim(selectivity_ytrsfl)
  
  numbers_cpue_ytrsl <- results_Uraw$x2[[ii]]$numbers_cpue_ytrsl
  dim(numbers_cpue_ytrsl) # it does not have fishery 
  
  # function to calculate predicted projected CPUE
  tmp <- get_cpue_proj(data, parameters, log_q_3 = log_q_3, log_pow_2 = log_pow_2, n_proj = n_proj, numbers_cpue_ytrsl=numbers_cpue_ytrsl,
                       selectivity_ytrsfl = selectivity_ytrsfl) 
  # add iteration column
  tmp$Iteration <- ii
  
  # store 
  cpue_proj_x2[[ii]] <- tmp
  
}
cpue_proj_x2 <- dplyr::bind_rows(cpue_proj_x2, .id = "Iteration")
dim(cpue_proj_x2)

# Scenario x3 U
cpue_proj_x3 <- list()
for(ii in 1:n_iter) {
  
  log_q_3 <- params$log_q[ii]
  log_pow_2 <- params$log_pow[ii]
  
  selectivity_ytrsfl <- proj_inputs$selectivity_ytrsfl[ii, , , , , , , drop = FALSE]
  selectivity_ytrsfl <- array(selectivity_ytrsfl, dim = dim(selectivity_ytrsfl)[-1])
  dim(selectivity_ytrsfl)
  
  numbers_cpue_ytrsl <- results_Uraw$x3[[ii]]$numbers_cpue_ytrsl
  dim(numbers_cpue_ytrsl) # it does not have fishery 
  
  # function to calculate predicted projected CPUE
  tmp <- get_cpue_proj(data, parameters, log_q_3 = log_q_3, log_pow_2 = log_pow_2, n_proj = n_proj, numbers_cpue_ytrsl=numbers_cpue_ytrsl,
                       selectivity_ytrsfl = selectivity_ytrsfl) 
  # add iteration column
  tmp$Iteration <- ii
  
  # store 
  cpue_proj_x3[[ii]] <- tmp
  
}
cpue_proj_x3 <- dplyr::bind_rows(cpue_proj_x3, .id = "Iteration")
dim(cpue_proj_x3)

cpue_proj_same <- cpue_proj_same %>% mutate(Scenario = "same")
cpue_proj_x2   <- cpue_proj_x2   %>% mutate(Scenario = "x2")
cpue_proj_x3   <- cpue_proj_x3   %>% mutate(Scenario = "x3")

cpue_proj <- bind_rows(cpue_proj_same, cpue_proj_x2, cpue_proj_x3)
head(cpue_proj)
cpue_proj <- cpue_proj %>%
  mutate(Year = Year + 2024, # replace 1 by 2025 and so on 
    Season  = l_season[Season],   # replace 1 and 2 by "AW" and "SS"
    Fishery = l_fishery[Fishery] ) # replace 1 and 2 by "SL" and "NSL"
  
head(cpue_proj)

cpue_hist<- data.frame(Year=data$cpue_year, Season = data$cpue_season, Fishery = data$cpue_fishery, 
                       cpue_obs = data$cpue_obs, cpue_sd = data$cpue_sd, q = data$cpue_q) %>%
            mutate(Year = l_year[Year], Season  = l_season[Season], Fishery = l_fishery[Fishery], Scenario = "Assessment")   

# PLOT CPUE
cpue_plot <- cpue_hist %>% filter (q==3) %>% mutate(Iteration = as.character(1)) %>%
  select(Year, Season, Fishery, cpue_obs, cpue_sd, Iteration, Scenario)

cpue_tmp <- cpue_proj %>% filter (Fishery == "SL") %>%
  select(Year, Season, Fishery, cpue_obs, cpue_sd, Iteration, Scenario)

cpue_plot <- bind_rows(cpue_plot, cpue_tmp)

# ---- Plot ----
ggplot(cpue_plot, aes(x = Year, y = cpue_obs, color = Scenario, fill = Scenario)) +
  stat_summary(fun.data = function(x) {data.frame(y = mean(x, na.rm=TRUE),
               ymin = quantile(x, 0.05, na.rm=TRUE),
               ymax = quantile(x, 0.95, na.rm=TRUE))},
  geom = "ribbon", alpha = 0.3, color = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  facet_wrap(~Season, scales = "free_y") +
  geom_vline(xintercept = 2024.5, linetype = "dashed") +
  labs(y = "Logbook CPUE", x = "Fishing year") +
  theme(legend.position = "right")

#ggsave("obs&proj_CPUE.pdf", width = 10, height = 6)
cpue_input <- cpue_plot
save(cpue_input, file = "Input_cpue.rda")

## adjusted vulnerable biomass (by season) and B0 to compare outputs in the EM phase 
head(df_UB)
names(df_list)
head(df_list$B0)

VB_true <- df_UB %>%
  left_join(df_list$B0 %>% 
              select(Iteration, B0 = value),by = "Iteration") 
head(VB_true)
tail(VB_true)
save(VB_true, file = "Adj_VB_compare.rda")

