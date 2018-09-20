library(tidyverse)
library(xtable)

load("output/output/simulation_intervals.Rdata") 
load("output/output/sim_data.Rdata")

sim_data <- data_list$sim_data

# Extract intervals
obtain_intervals <- function(interval_list, sim_data){
  n_sims <- length(interval_list)
  n_dsets <- length(interval_list[[1]])
  combined_iterations <- lapply(1:n_sims, function(i){
    merged <- lapply(1:n_dsets, function(j){
      df_sim <- sim_data[[i]][[j]]
      df_interval <- interval_list[[i]][[j]]
      df_interval <- cbind(df_sim, df_interval)
      df_interval <- df_interval %>%
        mutate(sim = i, dataset = j)
      df_interval$area_num <- 1:nrow(df_interval)
      df_interval
    })
    do.call(rbind, merged)
  })
  combined_iterations
}

results_list <- obtain_intervals(interval_list, sim_data)
save(results_list, file = "output/results_list.Rdata")

load("output/results_list.Rdata")


# Convert to tidy data
ci_data <- lapply(results_list, function(sim){
  sim  <- sim %>%
    dplyr::select(-X1, -dir_est, -var_est, -samp_size)
  results <- sim %>%
    gather(interval, value, direct_int.lower:eb_nonspatial_nocov.upper) %>%
    #      select(id, interval, value) %>%
    mutate(type = case_when(
      grepl("direct", interval) ~ "Direct",
      grepl("fab", interval) ~ "FAB",
      T ~ "EB"
    ), spatial = case_when(
      grepl("nonspatial", interval) ~ "Non-Spatial",
      grepl("direct", interval) ~ "N/A",
      T ~ "Spatial"
    ), covariate =  case_when(
      grepl("nocov", interval) ~ "No Covariate",
      grepl("direct", interval) ~ "N/A",
      T ~ "Covariate"
    ), endpoint = case_when(
      grepl("lower", interval) ~ "Lower",
      T ~ "Upper"
    )) %>%
    spread(endpoint, value) %>%
    group_by(area_num, dataset, type, spatial, covariate) %>%
    mutate(Lower = mean(Lower, na.rm = T), Upper = mean(Upper, na.rm = T)) %>%
    ungroup() %>%
    mutate(length = Upper - Lower) %>%
    mutate(coverage = (theta < Upper) & (theta > Lower))
  results
})

save(ci_data, file = "output/ci_tidy.Rdata")
load("output/ci_tidy.Rdata")

# Extract lengths of intervals
ci_length <- lapply(ci_data, function(sim){
  sim %>%
    select(-interval) %>%
    distinct() %>%
    group_by(dataset) %>%
    mutate(n_unstable = sum(length > 3.92 & type == "EB")) %>%
    filter(n_unstable == 0) %>%
    select(-n_unstable) %>%
    ungroup() %>%
    group_by(type, spatial, covariate) %>%
    summarise(mean_length = mean(length))
})

ci_length_better <- lapply(ci_data, function(sim){
  sim %>%
    select(-interval) %>%
    distinct() %>%
    group_by(dataset) %>%
    mutate(n_unstable = sum(length > 3.92 & type == "EB")) %>%
    filter(n_unstable == 0) %>%
    select(-n_unstable) %>%
    ungroup() %>%
    group_by(type, spatial, covariate) %>%
    summarise(mean_length = mean(length < 3.92))
})

save(ci_length, file = "output/ci_length.Rdata")

coverage_plot_list <- lapply(1:length(
  ci_data), function(i){
  df <- ci_data[[i]]
  beta <- data_list$params[i, "beta"]
  df <- df %>%
    filter(spatial == "Spatial" & covariate == "Covariate") %>%
    filter(type != "Direct") %>%
    select(-interval) %>%
    distinct() %>%
    mutate(tmxb = theta - X1 * beta) %>%
    filter(tmxb < quantile(tmxb, 0.98) & tmxb > quantile(tmxb, 0.02)) %>%
    na.omit()
  
  xlimitbds <- quantile(df$tmxb, c(0.05, 0.95))
ggplot(df) +
  stat_summary_bin(aes(x= tmxb, y = as.numeric(coverage), color = type), 
                   fun.data='mean_cl_normal',fun.args = list(mult=2), bins = 30,
                   geom = "errorbar") +
  labs(x = expression(theta[j] - X*beta), y = "Interval Coverage") +
  geom_hline(yintercept = 0.95, lty = 2) +
  xlim(xlimitbds) + 
  scale_colour_discrete(name = "Interval type",
                      labels = c("Spatial EB", "Spatial FAB")) +
  theme(legend.position = "top", axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))
})

save(coverage_plot_list, file = "output/coverage_plot_list.Rdata")

params <- data_list$params

# Length Tables
lengthtable_df <- lapply(1:length(ci_length), function(i){
  ci_length[[i]] %>%
    mutate(intercept = params[i, "intercept"], beta = params[i, "beta"],
           rho = params[i, "rho"], tau2 = params[i, "tau2"]) %>%
    select(beta, rho, tau2, type, spatial, covariate, mean_length)
})

lengthtable_df <- lapply(1:length(ci_length), function(i){
  namecol <- paste0("beta", params[i, "beta"], "rho", params[i, "rho"], "tau", params[i, "tau2"])
  df <- ci_length[[i]] %>%
    mutate(intercept = params[i, "intercept"], beta = params[i, "beta"],
           rho = params[i, "rho"], tau2 = params[i, "tau2"]) %>%
    unite(interval, type, spatial, covariate, sep = "_") %>%
    select(interval, mean_length)
  rownames(df) <- df$interval
  df %>%
    select(mean_length)
})

lengthtable_df <- do.call(cbind, lengthtable_df)

colnames(lengthtable_df) <- sapply(1:length(ci_length), function(i){paste0("beta", params[i, "beta"],
                                                                           "rho", params[i, "rho"], "tau", 
                                                                           params[i, "tau2"]) })

xtable(lengthtable_df, digits = 3)

table1 <- lengthtable_df %>%
  filter(spatial != "Spatial")

print(xtable(table1, digits = 3), include.rownames = F)

table2 <- lengthtable_df %>%
  filter(spatial != "Non-Spatial")

print(xtable(table2, digits = 3), include.rownames = F)

print(xtable(lengthtable_df, digits = 3), include.rownames = F)

# Coverage Figures
coverage_plot_list[[1]]
ggsave("figures/fig_coverage_lowtau.pdf", device = "pdf")
coverage_plot_list[[7]]
ggsave("figures/fig_coverage_hightau.pdf", device = "pdf")
