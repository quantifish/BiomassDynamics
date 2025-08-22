get_cpue_proj <- function(data, parameters, log_q_3, log_pow_2, n_proj = n_proj, numbers_cpue_ytrsl=numbers_cpue_ytrsl,
                          selectivity_ytrsfl = selectivity_ytrsfl) {
  log_q_3 <- log_q_3
  log_pow_2 <- log_pow_2
  getAll(data, parameters, warn = FALSE)
  log_q[3] <- log_q_3
  n_year <- data$n_year
  n_proj <- n_proj
  n_tot <- n_year + n_proj
  nq <- length(log_q)
  q_yi <- matrix(1, nrow = n_tot, ncol = nq)
  for (i in 1:nq) {
    for (y in unique(cpue_first_qcreep)[i] - 1) {
      q_yi[y, i] <- exp(log_q[i]) # remove?
    }
    for (y in unique(cpue_first_qcreep)[i]:n_tot) {
      q_yi[y, i] <- q_yi[y - 1, i] * (1 + qcreep[i])
    }
  }
  q_yi <- q_yi[c((n_year+1):n_tot),] #projected
  
  n_cpue <- n_proj*2*2 # I need it by season and fishery
  cpue_pred <- numeric(n_cpue) 
  cpue_year <- rep(1:n_proj, each = 2, times = 2)
  cpue_season <- rep(1:2, times = n_proj*2)
  cpue_region <- rep(1, times = length(cpue_year))
  #cpue_sex <- rep(1:2, times = n_proj*2)
  cpue_fishery <- rep(c(1, 2), each = n_proj*2)
  cpue_q <- rep(3, times = length(cpue_year)) # only for logbook cpue
  for (i in 1:n_cpue) {
      iy <- cpue_year[i]
      it <- cpue_season[i]
      ir <- cpue_region[i]
      ix <- cpue_sex[i]
      ij <- cpue_fishery[i]
      iq <- cpue_q[i]
      # units == 1 for logbook CPUE and sexes combined: 
        cpue_tmp <- numbers_cpue_ytrsl[iy, it, ir, ,] * selectivity_ytrsfl[iy, it, ir, , ij,]
        cpue_tmp <- cpue_tmp * legal_ytrsfl[n_year, it, ir, , ij,] # legal equal to the last year 
        cpue_pred[i] <- q_yi[iy, iq] * sum(cpue_tmp)^exp(log_pow_2) + 1e-6
  }
  cpue_pro <-  cpue_pro[length(cpue_pro)] # last value = 0
  cpue_wt <- cpue_wt[length(cpue_wt)]  # # last value = 1.31 for logbook cpue
  sd_l <- logbook_CRA2$SD   # observed SD for logbook CPUE
  # Sample cpue_sd values for each cpue_pred
  cpue_sd <- sample(sd_l, size = length(cpue_pred), replace = TRUE)
  cpue_sigma <- sqrt(cpue_sd^2 + cpue_pro^2) / cpue_wt
  cpue_obs <- rlnorm(n = length(cpue_pred),
                     meanlog = log(cpue_pred),
                     sdlog   = cpue_sigma)      
  results <- data.frame(Year = cpue_year, Season= cpue_season, Fishery = cpue_fishery, 
                        cpue_obs = cpue_obs, cpue_sd = cpue_sigma)
        return(results)
}
