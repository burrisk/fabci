# Load Outside Scripts
source("code/generate_data.R")
source("code/area_ci.R")

# Create Datasets
set.seed(315)
rows <- 7
cols <- 7
d <- 1
p <- rows*cols
R <- createNeighborhoodMat(rows, cols)
W <- R
diag(W) <- 0
W <- sweep(W, 1, rowSums(W), "/")
n_dsets <- 5000
n_min <- 1
n_max <- 1


# Tuning Parameters
rho <- c(0, 0.9)
intercept <- c(0)
beta <- c(0, 10)
tau2 <- c(0.5, 5)
params <- as.matrix(expand.grid(rho, intercept, beta, tau2))
colnames(params) <- c("rho", "intercept", "beta", "tau2")

set.seed(1)
sim_data <- lapply(1:nrow(params), function(i){
  lapply(1:n_dsets, function(j){
    genDataSAR(p, W = W, beta = params[i, "beta"], 
               intercept = params[i, "intercept"], rho = params[i, "rho"], tau2 = params[i, "tau2"],
               var = 1, n_min = n_min, n_max = n_max)
  })
})

data_list <- list(W = W, params = params, sim_data = sim_data)
save(data_list, file = "output/sim_data.Rdata")

load("output/sim_data.Rdata")
W <- data_list$W
params <- data_list$params
sim_data <- data_list$sim_data

interval_list <- lapply(sim_data, function(sim){
  lapply(1:length(sim), function(k){
    print(paste(c("Dataset", k, "of", n_dsets), sep = " "))
    df_sim <- sim[[k]]
    tryCatch({
      direct_int <- areaCI(formula = dir_est ~ X1, data = df_sim, W = W, direct_var = df_sim$var_est, 
             alpha = 0.05, method = "Direct", type = "FH", varknown = T, maxiter = 20)
      fab_spatial <- areaCI(formula = dir_est ~ X1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                            alpha = 0.05, method = "FAB", type = "SAR", varknown = T, maxiter = 20)
      fab_nonspatial <- areaCI(formula = dir_est ~ X1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                            alpha = 0.05, method = "FAB", type = "FH", varknown = T, maxiter = 20)
      fab_spatial_nocov <- areaCI(formula = dir_est ~ 1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                            alpha = 0.05, method = "FAB", type = "SAR", varknown = T, maxiter = 20)
      fab_nonspatial_nocov <- areaCI(formula = dir_est ~ 1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                               alpha = 0.05, method = "FAB", type = "FH", varknown = T, maxiter = 20)
      eb_spatial <- areaCI(formula = dir_est ~ X1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                            alpha = 0.05, method = "EB", type = "SAR", varknown = T)
      eb_nonspatial <- areaCI(formula = dir_est ~ X1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                               alpha = 0.05, method = "EB", type = "FH", varknown = T)
      eb_spatial_nocov <- areaCI(formula = dir_est ~ 1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                                  alpha = 0.05, method = "EB", type = "SAR", varknown = T)
      eb_nonspatial_nocov <- areaCI(formula = dir_est ~ 1, data = df_sim, W = W, direct_var = df_sim$var_est, 
                                     alpha = 0.05, method = "EB", type = "FH", varknown = T)
      list(direct_int = direct_int, fab_spatial = fab_spatial, fab_nonspatial = fab_nonspatial,
           fab_spatial_nocov = fab_spatial_nocov, fab_nonspatial_nocov = fab_nonspatial_nocov,
           eb_spatial = eb_spatial, eb_nonspatial = eb_nonspatial, eb_spatial_nocov = eb_spatial_nocov,
           eb_nonspatial_nocov = eb_nonspatial_nocov)
    },
    error = function(e) matrix(0, nrow = p, ncol = 9)
    )
  })
})

save(interval_list, file = "output/simulation_intervals.Rdata")

