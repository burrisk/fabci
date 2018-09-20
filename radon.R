# Load Packages ----
library(tidyverse)
library(xtable)
library(mvtnorm)

# Load Scripts
source("code/fab_functions.R")
source("code/area_ci.R")

# Load Data ----
radon <- read_csv("Data/SSRS.csv")
radium <- read_csv("Data/radium.csv")

# Clean Data ----
radium <- radium %>%
  rename(cntyfips = ctfips, State = st) %>%
  group_by(cntyfips, State, cty, Uppm) %>%
  summarise(lon = median(lon), lat = median(lat))

nrow(radon)
radon %>%
  group_by(state2, cntyfips) %>%
  filter(row_number() == 1) %>%
  nrow()
radon %>%
  group_by(state2) %>%
  filter(row_number() == 1) %>%
  nrow()

radon <- radon %>%
  dplyr::select(state2, cntyfips, county, activity) %>%
  rename(State = state2) %>%
  left_join(radium) %>%
  rename(County = cty) %>%
  filter(State %in% c("MN", "WI", "MI", "IN")) %>%
  mutate(activity = activity/2 + sqrt(activity ^ 2 / 4 + 0.25 ^ 2)) %>%
  mutate(lactivity = log(activity))

nrow(radon)
radon %>%
  group_by(State, cntyfips) %>%
  filter(row_number() == 1) %>%
  nrow()


g_hist <- ggplot(radon, aes(x = lactivity)) +
  geom_histogram(binwidth = 0.2, color = "black", fill = "orange") +
  xlab("Log radon concentration (pCi/L)") +
  ylab("Number of Households")

ggsave("output/lactivity_hist.png", g_hist, device = "png", width = 7, height = 5)

radon <- radon %>%
  group_by(State, County, cntyfips, lon, lat, Uppm) %>%
  summarise(mean_lradon = mean(lactivity), sd_lradon = sd(lactivity),
            sample_size = n()) %>%
  ungroup() %>%
  mutate(var = sd_lradon ^ 2 / sample_size) %>%
  droplevels() %>%
  na.omit()
  
# Calculate confidence intervals
lonlatmat <- radon %>%
  dplyr::select(lon, lat)

distmat <- exp(-as.matrix(dist(lonlatmat, diag = T, upper = T)^2))
diag(distmat) <- 0
distmat <- sweep(distmat, 1, rowSums(distmat), "/")

alpha = 0.05

radon <- as.data.frame(radon)

# Direct Interval
direct_int <- areaCI(formula = mean_lradon ~ 1, data = radon, direct_var = radon$var,
                    sample_size = radon$sample_size, alpha = 0.05, method = "Direct",
                    varknown = F)

# FAB Intervals
fab_nonspatial_nocov <- areaCI(formula = mean_lradon ~ 1, data = radon, direct_var = radon$var,
                    sample_size = radon$sample_size, alpha = 0.05, method = "FAB", type = "FH",
                    varknown = F)

fab_spatial_nocov <- areaCI(formula = mean_lradon ~ 1, data = radon, W = distmat, direct_var = radon$var,
                                 sample_size = radon$sample_size, alpha = 0.05, method = "FAB",
                                 type = "SAR", varknown = F)

fab_nonspatial <- areaCI(formula = mean_lradon ~ Uppm, data = radon, direct_var = radon$var,
                                 sample_size = radon$sample_size, alpha = 0.05, method = "FAB", type = "FH",
                                 varknown = F)

fab_spatial <- areaCI(formula = mean_lradon ~ Uppm, data = radon, W = distmat, direct_var = radon$var,
                              sample_size = radon$sample_size, alpha = 0.05, method = "FAB",
                              type = "SAR", varknown = F)

# Empirical Bayes Intervals
eb_nonspatial_nocov <- areaCI(formula = mean_lradon ~ 1, data = radon, direct_var = radon$var,
                                 sample_size = radon$sample_size, alpha = 0.05, method = "EB", type = "FH",
                                 varknown = F)

eb_spatial_nocov <- areaCI(formula = mean_lradon ~ 1, data = radon, W = distmat, direct_var = radon$var,
                              sample_size = radon$sample_size, alpha = 0.05, method = "EB",
                              type = "SAR", varknown = F)

eb_nonspatial <- areaCI(formula = mean_lradon ~ Uppm, data = radon, direct_var = radon$var,
                           sample_size = radon$sample_size, alpha = 0.05, method = "EB", type = "FH",
                           varknown = F)

eb_spatial <- areaCI(formula = mean_lradon ~ Uppm, data = radon, W = distmat, direct_var = radon$var,
                        sample_size = radon$sample_size, alpha = 0.05, method = "EB",
                        type = "SAR", varknown = F)

# Aggregate Results
radon_intervals <- list(direct_int = direct_int, fab_spatial = fab_spatial, fab_nonspatial = fab_nonspatial,
                  fab_spatial_nocov = fab_spatial_nocov, fab_nonspatial_nocov = fab_nonspatial_nocov,
                  eb_spatial = eb_spatial, eb_nonspatial = eb_nonspatial, eb_spatial_nocov = eb_spatial_nocov,
                  eb_nonspatial_nocov = eb_nonspatial_nocov)

radonint <- cbind(radon, radon_intervals)

results <- radonint %>%
  mutate(id = 1:nrow(radonint)) %>%
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
  group_by(id, type, spatial, covariate) %>%
  mutate(Lower = mean(Lower, na.rm = T), Upper = mean(Upper, na.rm = T)) %>%
  ungroup() %>%
  mutate(length = Upper - Lower) %>%
  select(-interval) %>%
  distinct()

save(results, file = "output/radon_gl_intervals.Rdata")

load("output/radon_gl_intervals.Rdata")

direct <- results %>%
  filter(type == "Direct") %>%
  arrange(cntyfips)

fab_exchangeable <- results %>%
  filter(type == "FAB" & covariate == "No Covariate" & spatial == "Non-Spatial") %>%
  arrange(cntyfips)

fab_covariate <- results %>%
  filter(type == "FAB" & covariate == "Covariate" & spatial == "Non-Spatial") %>%
  arrange(cntyfips)

fab_spatial <- results %>%
  filter(type == "FAB" & covariate == "No Covariate" & spatial == "Spatial") %>%
  arrange(cntyfips)

fab_full <- results %>%
  filter(type == "FAB" & covariate == "Covariate" & spatial == "Spatial") %>%
  arrange(cntyfips)

eb_exchangeable <- results %>%
  filter(type == "EB" & covariate == "No Covariate" & spatial == "Non-Spatial") %>%
  arrange(cntyfips)

eb_covariate <- results %>%
  filter(type == "EB" & covariate == "Covariate" & spatial == "Non-Spatial") %>%
  arrange(cntyfips)

eb_spatial <- results %>%
  filter(type == "EB" & covariate == "No Covariate" & spatial == "Spatial") %>%
  arrange(cntyfips)

eb_full <- results %>%
  filter(type == "EB" & covariate == "Covariate" & spatial == "Spatial") %>%
  arrange(cntyfips)

# Improvements over direct interval
mean((fab_exchangeable$length - direct$length) < 0)
mean((fab_covariate$length - direct$length) < 0)
mean((fab_spatial$length - direct$length) < 0)
mean((fab_full$length - direct$length) < 0)

length_table <- results %>%
  group_by(type, spatial, covariate) %>%
  summarise(mean_length = mean(length))

print(xtable(length_table, digits = 3), include.rownames = F)

length_table$mean_length / 1.701

# Plot Direct and FAB Intervals
results_plot <- results %>%
  filter(type == "FAB" | type == "Direct",
         spatial != "Non-Spatial",
         covariate != "No Covariate") %>%
  rename(ybar = mean_lradon)

g1 <- ggplot(results_plot) +
    geom_segment(mapping = aes(x = ybar, y = Lower, xend = ybar, yend =Upper,
                             colour = type), size = 1, alpha = 0.7) +
  xlab(expression(bar(y))) +
  ylab("Confidence Interval") +
  scale_color_discrete(name = 'Interval Type', labels = c("Direct", "FAB (Spatial)"))
ggsave(file = "figures/radon_compare.pdf", g1, device = "pdf", width = 6, height = 5)
