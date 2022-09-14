# This code implements the simulations described in Picardi et al., 
# "Defining null expectations for animal site fidelity". Package fidelity 
# is available on GitHub (https://github.com/picardis/fidelity) and can be 
# installed as follows: remotes::install_github("picardis/fidelity")

# Load packages ####

library(fidelity)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(raster)

# Save files to disk? ####

save <- FALSE

# Output folder path ####

path_out <- "output/"

# Set seed ####

set.seed(7)

# Generate landscapes ####

lands <- sim_land(res = 200, 
                  ext = 100000, 
                  orig = c(x = 292332, y = 4642013),
                  kappa = seq(1, 5000, length.out = 3), 
                  chi = seq(0.1, 1.5, length.out = 3))

if(save) {
  raster::writeRaster(lands, 
                      filename = paste0(path_out, "Landscapes/land.tif"), 
                      bylayer = TRUE, 
                      suffix = names(lands))
}

land_names <- paste0(path_out, 
                     "Landscapes/land_",
                     names(lands), 
                     ".tif")

# Create simulation scenarios ####

# CRWs

scen_crw <- create_scenarios_crw(rho = seq(0, 0.95, length.out = 4))

# CCRWs

scen_ccrw <- create_scenarios_ccrw(rho = 0.2, 
                                   boundary_size = c(1, 5, 10, 20) * 1000)

# MCRWs

scen_mcrw <- create_scenarios_mcrw(rho = 0.2,
                                   habitat_effect = seq(0, 1, length.out = 5)[2:5],
                                   lands = land_names)

# BCRWs

scen_bcrw <- create_scenarios_bcrw(rho = 0.2,
                                   beta = c(0.05, 0.15, 0.25, 0.5),
                                   lands = land_names)

# Bind
scen <- rbind(scen_crw, scen_ccrw, scen_mcrw, scen_bcrw)

# Calculate cell neighborhoods ####

# Define perceptual range as 95th percentile of step length distribution
prange <- quantile(rweibull(n = 1000000, shape = 1, scale = 200), .95)

neigh <- get_neighbors(lands[[1]], prange = prange)

# Make cluster for parallelization ####

no_cores <- parallel::detectCores()-1
clust <- parallel::makeCluster(no_cores) 

parallel::clusterExport(clust, varlist = c("scen", 
                                           "neigh",
                                           "path_out"))

# Run simulations ####

system.time(sims <- do.call(
  "rbind", 
  parallel::parLapply(cl = clust, X = 1:nrow(scen), 
                      fun = function(x) {
                        res <- fidelity::simulate_tracks(
                          scenarios = scen[x, ],
                          # 2 years at 4h resolution 
                          n_steps = 6 * 365 * 2,
                          # 500 individuals
                          n_tracks = 500,
                          sl_par = c(1, 200), 
                          # this is the center of the raster
                          start_loc = data.frame(x = 292432, y = 4641913),
                          neighbors = neigh)
                        saveRDS(res, paste0(path_out, 
                                            "Simulation_results/scenario_" , 
                                            x, 
                                            ".rds"))
                        return(res)
                      })))

# Calculate return rate ####

# Load simulations back in
f <- list.files(paste0(path_out, "Simulation_results/"), 
           pattern = "scenario", full.names = TRUE)

sims <- data.frame()

for (i in 1:length(f)) {
  print(i)
  s <- readRDS(f[i])
  sims <- rbind(sims, s)
  }

# Distance at which returns count is 50% percentile of step length
d <- quantile(rweibull(n = 100000, 1, 200), 0.5)

# Split up simulations by ID and scenario
sim_list <- split(sims, f = paste(sims$id, sims$scenario_id))

# Clear up environment
rm(list = ls()[!ls() %in% c("clust", "sim_list", "d", "path_out")])

# Set up cluster
no_cores <- parallel::detectCores()-1
clust <- parallel::makeCluster(no_cores) 
parallel::clusterExport(clust, varlist = c("sim_list", "d", "path_out"))

# Run
system.time(rets <- do.call(
  "rbind",
  parallel::parLapply(cl = clust, X = 1:length(sim_list),
                      fun = function(x) {
                        res <- fidelity::calc_returns(
                          tracks = sim_list[[x]],
                          dist = d,
                          lag = c(6 * 1, # one day
                                  6 * 14, # two weeks
                                  6 * 30, # one month
                                  6 * 90, # three months
                                  6 * 365)) # one year
                        write.table(res, paste0(path_out, 
                                                "Returns_calculated/returns_" , 
                                                x, 
                                                ".txt"), sep = ",")
                        return(res)
                      })
))

write.table(rets, paste0(path_out, "Returns_calculated/all_returns.txt"), sep = ",")

# Reference simulation ####

# Compare and contrast scenarios with a reference rate of return arising from a 
# simple CRW with the same value of rho as the simulation but no boundary,
# no habitat effect, or no beta. 

scen_ref <- create_scenarios_crw(rho = 0.2)

ref <- fidelity::simulate_tracks(
                          scenarios = scen_ref,
                          n_steps = 6 * 365 * 2, # 2 years at 4h resolution 
                          n_tracks = 500,
                          sl_par = c(1, 200), 
                          start_loc = data.frame(x = 292432, y = 4641913))

ret_ref <- fidelity::calc_returns(
  tracks = ref,
  dist = d,
  lag = c(6 * 1, # one day
          6 * 14, # two weeks
          6 * 30, # one month
          6 * 90, # three months
          6 * 365)) # one year

write.table(ret_ref, paste0(path_out, "Returns_calculated/reference_returns.txt"), sep = ",")

# Visualization ####

rets <- read.table(paste0(path_out, "Returns_calculated/all_returns.txt"), sep = ",")
ret_ref <- read.table(paste0(path_out, "Returns_calculated/reference_returns.txt"), sep = ",")

# Left join to scenario data.frame to get back values of kappa and chi
rets <- left_join(rets, scen)

# Split by model
rets_ccrw <- rets %>% 
  filter(grepl("CCRW", scenario_id))
rets_bcrw <- rets %>% 
  filter(grepl("BCRW", scenario_id))
rets_mcrw <- rets %>% 
  filter(grepl("MCRW", scenario_id))
rets_crw <- rets %>% 
  filter(!grepl("CCRW", scenario_id) &
           !grepl("BCRW", scenario_id) &
           !grepl("MCRW", scenario_id))

## Rearrange in plottable format ####

### Reference scenario ####

fid_ref <- ret_ref %>%     
  group_by(step) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(boundary_size = NA_integer_, 
         habitat_effect = NA_real_,
         beta = NA_real_,
         kappa = NA_real_,
         chi = NA_real_,
         lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

fid_ref$lag_words <- factor(fid_ref$lag_words,
                            levels = c("1 Day Lag", 
                                       "2 Weeks Lag",
                                       "1 Month Lag",
                                       "3 Months Lag",
                                       "1 Year Lag"))

### CRW ####

# Get mean returns across IDs at each step for each scenario
rr_crw <- rets_crw %>%     
  group_by(step, rho) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE))

# Format for plotting
f_crw <- rr_crw %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(rho = round(rho, 2), 
         lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

# Spell out lags
f_crw$lag_words <- factor(f_crw$lag_words,
                          levels = c("1 Day Lag", 
                                     "2 Weeks Lag",
                                     "1 Month Lag",
                                     "3 Months Lag",
                                     "1 Year Lag"))

### CCRW ####

# Get mean returns across IDs at each step for each scenario
rr_ccrw <- rets_ccrw %>%     
  group_by(step, boundary_size) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_ccrw <- rr_ccrw %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6),
         boundary_size = boundary_size/1000) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_ccrw$lag_words <- factor(f_ccrw$lag_words,
                           levels = c("1 Day Lag", 
                                      "2 Weeks Lag",
                                      "1 Month Lag",
                                      "3 Months Lag",
                                      "1 Year Lag"))

f_ccrw$boundary_size <- factor(f_ccrw$boundary_size,
                               levels = c("1", "5", "10", "20"))

### MCRW ####

# Get mean returns across IDs at each step for each scenario
rr_mcrw <- rets_mcrw %>%     
  group_by(step, habitat_effect) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_mcrw <- rr_mcrw %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_mcrw$habitat_effect <- factor(f_mcrw$habitat_effect)

f_mcrw$lag_words <- factor(f_mcrw$lag_words,
                           levels = c("1 Day Lag", 
                                      "2 Weeks Lag",
                                      "1 Month Lag",
                                      "3 Months Lag",
                                      "1 Year Lag"))

### BCRW ####

# Get mean returns across IDs at each step for each scenario
rr_bcrw <- rets_bcrw %>%     
  group_by(step, beta) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_bcrw <- rr_bcrw %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_bcrw$beta <- factor(f_bcrw$beta)

f_bcrw$lag_words <- factor(f_bcrw$lag_words,
                           levels = c("1 Day Lag", 
                                      "2 Weeks Lag",
                                      "1 Month Lag",
                                      "3 Months Lag",
                                      "1 Year Lag"))

### MCRW + chi ####

# Get mean returns across IDs at each step for each scenario
rr_mcrw_chi <- rets_mcrw %>%     
  group_by(step, chi) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_mcrw_chi <- rr_mcrw_chi %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_mcrw_chi$lag_words <- factor(f_mcrw_chi$lag_words,
                              levels = c("1 Day Lag", 
                                         "2 Weeks Lag",
                                         "1 Month Lag",
                                         "3 Months Lag",
                                         "1 Year Lag"))

f_mcrw_chi$chi <- factor(f_mcrw_chi$chi)

### MCRW + kappa ####

# Get mean returns across IDs at each step for each scenario
rr_mcrw_kappa <- rets_mcrw %>%     
  group_by(step, kappa) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_mcrw_kappa <- rr_mcrw_kappa %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_mcrw_kappa$lag_words <- factor(f_mcrw_kappa$lag_words,
                              levels = c("1 Day Lag", 
                                         "2 Weeks Lag",
                                         "1 Month Lag",
                                         "3 Months Lag",
                                         "1 Year Lag"))

f_mcrw_kappa$kappa <- factor(f_mcrw_kappa$kappa)

### BCRW + chi ####

# Get mean returns across IDs at each step for each scenario
rr_bcrw_chi <- rets_bcrw %>%     
  group_by(step, chi) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_bcrw_chi <- rr_bcrw_chi %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_bcrw_chi$lag_words <- factor(f_bcrw_chi$lag_words,
                              levels = c("1 Day Lag", 
                                         "2 Weeks Lag",
                                         "1 Month Lag",
                                         "3 Months Lag",
                                         "1 Year Lag"))

f_bcrw_chi$chi <- factor(f_bcrw_chi$chi)

### BCRW + kappa ####

# Get mean returns across IDs at each step for each scenario
rr_bcrw_kappa <- rets_bcrw %>%     
  group_by(step, kappa) %>% 
  summarize(across(contains("lag_"), mean, na.rm = TRUE)) %>% 
  ungroup() 

# Pivot longer to get in plottable format
f_bcrw_kappa <- rr_bcrw_kappa %>% 
  pivot_longer(cols = contains("lag_"), 
               names_to = "lag", values_to = "return_rate") %>% 
  mutate(lag = as.numeric(word(lag, 2, 2, "_")),
         day = ceiling(step/6)) %>% 
  ungroup() %>% 
  mutate(lag_words = factor(case_when(
    lag == 6 ~ "1 Day Lag",
    lag == 84 ~ "2 Weeks Lag",
    lag == 180 ~ "1 Month Lag",
    lag == 540 ~ "3 Months Lag",
    lag == 2190 ~ "1 Year Lag"
  ))) %>% 
  filter(lag == 6 & step >= 6 |
           lag == 84 & step >= 84 |
           lag == 180 & step >= 180 |
           lag == 540 & step >= 540 |
           lag == 2190 & step >= 2190)

f_bcrw_kappa$lag_words <- factor(f_bcrw_kappa$lag_words,
                              levels = c("1 Day Lag", 
                                         "2 Weeks Lag",
                                         "1 Month Lag",
                                         "3 Months Lag",
                                         "1 Year Lag"))

f_bcrw_kappa$kappa <- factor(f_bcrw_kappa$kappa)

## Example tracks ####

ggcolors <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
ggcolors_landscape <- c("#F8766D", "#00BA38", "#619CFF")

### Example tracks CRW ####

lims_crw <- rets_crw %>% 
  filter(scenario_id %in% c("CRW_1", "CRW_2", "CRW_3", "CRW_4")
         & id == 1) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y))

et_crw1 <- rets_crw %>% 
  filter(scenario_id == "CRW_1" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[1]) +
  coord_cartesian(xlim = c(lims_crw$minx, lims_crw$maxx), 
                  ylim = c(lims_crw$miny, lims_crw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_crw$maxx, y = lims_crw$maxy, hjust = 1, vjust = 1, 
           label = "rho == 0", parse = TRUE, fill = "white", alpha = 0.8)

et_crw2 <- rets_crw %>% 
  filter(scenario_id == "CRW_2" & id == 17) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[2]) +
  coord_cartesian(xlim = c(lims_crw$minx, lims_crw$maxx), 
                  ylim = c(lims_crw$miny, lims_crw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_crw$maxx, y = lims_crw$maxy, hjust = 1, vjust = 1, 
           label = "rho == 0.32", parse = TRUE, fill = "white", alpha = 0.8)

et_crw3 <- rets_crw %>% 
  filter(scenario_id == "CRW_3" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[3]) +
  coord_cartesian(xlim = c(lims_crw$minx, lims_crw$maxx), 
                  ylim = c(lims_crw$miny, lims_crw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_crw$maxx, y = lims_crw$maxy, hjust = 1, vjust = 1, 
           label = "rho == 0.63", parse = TRUE, fill = "white", alpha = 0.8)

et_crw4 <- rets_crw %>% 
  filter(scenario_id == "CRW_4" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[4]) +
  coord_cartesian(xlim = c(lims_crw$minx, lims_crw$maxx), 
                  ylim = c(lims_crw$miny, lims_crw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_crw$maxx, y = lims_crw$maxy, hjust = 1, vjust = 1, 
           label = "rho == 0.95", parse = TRUE, fill = "white", alpha = 0.8)

### Example tracks CCRW ####

lims_ccrw <- rets_ccrw %>% 
 summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y))

et_ccrw1 <- rets_ccrw %>% 
  filter(scenario_id == "CCRW_1" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[1]) +
  coord_cartesian(xlim = c(lims_ccrw$minx, lims_ccrw$maxx), 
                  ylim = c(lims_ccrw$miny, lims_ccrw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_ccrw$maxx, y = lims_ccrw$maxy, hjust = 1, vjust = 1, 
           label = "r = 1", fill = "white", alpha = 0.8)

et_ccrw2 <- rets_ccrw %>% 
  filter(scenario_id == "CCRW_2" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[2]) +
  coord_cartesian(xlim = c(lims_ccrw$minx, lims_ccrw$maxx), 
                  ylim = c(lims_ccrw$miny, lims_ccrw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_ccrw$maxx, y = lims_ccrw$maxy, hjust = 1, vjust = 1, 
           label = "r = 5", fill = "white", alpha = 0.8)

et_ccrw3 <- rets_ccrw %>% 
  filter(scenario_id == "CCRW_3" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[3]) +
  coord_cartesian(xlim = c(lims_ccrw$minx, lims_ccrw$maxx), 
                  ylim = c(lims_ccrw$miny, lims_ccrw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_ccrw$maxx, y = lims_ccrw$maxy, hjust = 1, vjust = 1, 
           label = "r = 10", fill = "white", alpha = 0.8)

et_ccrw4 <- rets_ccrw %>% 
  filter(scenario_id == "CCRW_4" & id == 1) %>% 
  ggplot() +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[4]) +
  coord_cartesian(xlim = c(lims_ccrw$minx, lims_ccrw$maxx), 
                  ylim = c(lims_ccrw$miny, lims_ccrw$maxy)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_ccrw$maxx, y = lims_ccrw$maxy, hjust = 1, vjust = 1, 
           label = "r = 20", fill = "white", alpha = 0.8)

### Example tracks MCRW ####

lims_mcrw <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           ((id == 1 & habitat_effect == 0.25) |
              (id == 3 & habitat_effect == 0.50) |
              (id == 5 & habitat_effect == 0.75) |
              (id == 5 & habitat_effect == 1.00))) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_mcrw <- data.frame(diffx = max(c(lims_mcrw$diffx, lims_mcrw$diffy)),
                        diffy = max(c(lims_mcrw$diffx, lims_mcrw$diffy)))

land_m <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_2500.5_chi_0.8.tif")), 
  xy = TRUE)

lims_m1 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_mcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_mcrw$diffx/2 + 5000,
         limy_low = center_y - lims_mcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_mcrw$diffy/2 + 5000)

et_mcrw1 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m1$limx_low, lims_m1$limx_hig), 
                  ylim = c(lims_m1$limy_low, lims_m1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m1$limx_hig, y = lims_m1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "eta == 0.25", parse = TRUE, fill = "white", alpha = 0.8)

lims_m2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.50 &
           id == 3) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_mcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_mcrw$diffx/2 + 5000,
         limy_low = center_y - lims_mcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_mcrw$diffy/2 + 5000)

et_mcrw2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.50 &
           id == 3) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y),  color = ggcolors[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m2$limx_low, lims_m2$limx_hig), 
                  ylim = c(lims_m2$limy_low, lims_m2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_m2$limx_hig, y = lims_m2$limy_hig, hjust = 1, vjust = 1, 
           label = "eta == 0.5", parse = TRUE, fill = "white", alpha = 0.8)

lims_m3 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.75 &
           id == 5) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_mcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_mcrw$diffx/2 + 5000,
         limy_low = center_y - lims_mcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_mcrw$diffy/2 + 5000)

et_mcrw3 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.75 &
           id == 5) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m3$limx_low, lims_m3$limx_hig), 
                  ylim = c(lims_m3$limy_low, lims_m3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_m3$limx_hig, y = lims_m3$limy_hig, hjust = 1, vjust = 1, 
           label = "eta == 0.75", parse = TRUE, fill = "white", alpha = 0.8)

lims_m4 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 1.00 &
           id == 5) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_mcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_mcrw$diffx/2 + 5000,
         limy_low = center_y - lims_mcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_mcrw$diffy/2 + 5000)

et_mcrw4 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 1.00 &
           id == 5) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[4]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m4$limx_low, lims_m4$limx_hig), 
                  ylim = c(lims_m4$limy_low, lims_m4$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", x = lims_m4$limx_hig, y = lims_m4$limy_hig, hjust = 1, vjust = 1, 
           label = "eta == 1", parse = TRUE, fill = "white", alpha = 0.8)

### Example tracks BCRW ####

lims_bcrw <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           id == 1 & !is.na(beta)) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_bcrw <- data.frame(diffx = max(c(lims_bcrw$diffx, lims_bcrw$diffy)),
                        diffy = max(c(lims_bcrw$diffx, lims_bcrw$diffy)))

lims_b1 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.05 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_bcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_bcrw$diffx/2 + 5000,
         limy_low = center_y - lims_bcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_bcrw$diffy/2 + 5000)

et_bcrw1 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.05 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b1$limx_low, lims_b1$limx_hig), 
                  ylim = c(lims_b1$limy_low, lims_b1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b1$limx_hig, y = lims_b1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "beta == 0.05", parse = TRUE, fill = "white", alpha = 0.8)

lims_b2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_bcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_bcrw$diffx/2 + 5000,
         limy_low = center_y - lims_bcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_bcrw$diffy/2 + 5000)

et_bcrw2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b2$limx_low, lims_b2$limx_hig), 
                  ylim = c(lims_b2$limy_low, lims_b2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b2$limx_hig, y = lims_b2$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "beta == 0.15", parse = TRUE, fill = "white", alpha = 0.8)

lims_b3 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_bcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_bcrw$diffx/2 + 5000,
         limy_low = center_y - lims_bcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_bcrw$diffy/2 + 5000)

et_bcrw3 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b3$limx_low, lims_b3$limx_hig), 
                  ylim = c(lims_b3$limy_low, lims_b3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b3$limx_hig, y = lims_b3$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "beta == 0.25", parse = TRUE, fill = "white", alpha = 0.8)

lims_b4 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.5 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_bcrw$diffx/2 - 5000,
         limx_hig = center_x + lims_bcrw$diffx/2 + 5000,
         limy_low = center_y - lims_bcrw$diffy/2 - 5000,
         limy_hig = center_y + lims_bcrw$diffy/2 + 5000)

et_bcrw4 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.5 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m, mapping = aes(x = x, y = y, fill = land_m[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors[4]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b4$limx_low, lims_b4$limx_hig), 
                  ylim = c(lims_b4$limy_low, lims_b4$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b4$limx_hig, y = lims_b4$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "beta == 0.5", parse = TRUE, fill = "white", alpha = 0.8)

### Example tracks landscape effects (MCRW and BCRW) ####

land_m1 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_1_chi_0.8.tif")), 
  xy = TRUE)

land_m2 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_2500.5_chi_0.8.tif")), 
  xy = TRUE)

land_m3 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_5000_chi_0.8.tif")), 
  xy = TRUE)

lims_m_kappa <- rets_mcrw %>% 
  filter((landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" |
            grepl("kappa_2500.5_chi_0.8.tif", landscape) |
            landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif") 
         & (id == 1 & habitat_effect == 0.25)) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_m_kappa <- data.frame(diffx = max(c(lims_m_kappa$diffx, lims_m_kappa$diffy)),
                          diffy = max(c(lims_m_kappa$diffx, lims_m_kappa$diffy)))

lims_m_kappa1 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_m_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_m_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_m_kappa$diffy/2 + 5000)

et_m_kappa1 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m1, mapping = aes(x = x, y = y, fill = land_m1[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_kappa1$limx_low, lims_m_kappa1$limx_hig), 
                  ylim = c(lims_m_kappa1$limy_low, lims_m_kappa1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_kappa1$limx_hig, y = lims_m_kappa1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 1", parse = TRUE, fill = "white", alpha = 0.8)

lims_m_kappa2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_m_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_m_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_m_kappa$diffy/2 + 5000)

et_m_kappa2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m2, mapping = aes(x = x, y = y, fill = land_m2[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_kappa2$limx_low, lims_m_kappa2$limx_hig), 
                  ylim = c(lims_m_kappa2$limy_low, lims_m_kappa2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_kappa2$limx_hig, y = lims_m_kappa2$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 2500", parse = TRUE, fill = "white", alpha = 0.8)

lims_m_kappa3 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_m_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_m_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_m_kappa$diffy/2 + 5000)

et_m_kappa3 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m3, mapping = aes(x = x, y = y, fill = land_m3[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_kappa3$limx_low, lims_m_kappa3$limx_hig), 
                  ylim = c(lims_m_kappa3$limy_low, lims_m_kappa3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_kappa3$limx_hig, y = lims_m_kappa3$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 5000", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_kappa <- rets_bcrw %>% 
  filter((landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" |
            grepl("kappa_2500.5_chi_0.8.tif", landscape) |
            landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif") 
         & (id == 1 & beta == 0.15)) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_b_kappa <- data.frame(diffx = max(c(lims_b_kappa$diffx, lims_b_kappa$diffy)),
                          diffy = max(c(lims_b_kappa$diffx, lims_b_kappa$diffy)))

lims_b_kappa1 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_b_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_b_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_b_kappa$diffy/2 + 5000)

et_b_kappa1 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_1_chi_0.8.tif" &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m1, mapping = aes(x = x, y = y, fill = land_m1[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_kappa1$limx_low, lims_b_kappa1$limx_hig), 
                  ylim = c(lims_b_kappa1$limy_low, lims_b_kappa1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b_kappa1$limx_hig, y = lims_b_kappa1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 1", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_kappa2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_b_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_b_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_b_kappa$diffy/2 + 5000)

et_b_kappa2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m2, mapping = aes(x = x, y = y, fill = land_m2[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_kappa2$limx_low, lims_b_kappa2$limx_hig), 
                  ylim = c(lims_b_kappa2$limy_low, lims_b_kappa2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b_kappa2$limx_hig, y = lims_b_kappa2$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 2500", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_kappa3 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif" &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_kappa$diffx/2 - 5000,
         limx_hig = center_x + lims_b_kappa$diffx/2 + 5000,
         limy_low = center_y - lims_b_kappa$diffy/2 - 5000,
         limy_hig = center_y + lims_b_kappa$diffy/2 + 5000)

et_b_kappa3 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_5000_chi_0.8.tif" &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_m3, mapping = aes(x = x, y = y, fill = land_m3[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_kappa3$limx_low, lims_b_kappa3$limx_hig), 
                  ylim = c(lims_b_kappa3$limy_low, lims_b_kappa3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b_kappa3$limx_hig, y = lims_b_kappa3$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "kappa == 5000", parse = TRUE, fill = "white", alpha = 0.8)

land_c1 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_2500.5_chi_0.1.tif")), 
  xy = TRUE)

land_c2 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_2500.5_chi_0.8.tif")), 
  xy = TRUE)

land_c3 <- as.data.frame(
  raster(
    paste0(path_out, "Landscapes/land_kappa_2500.5_chi_1.5.tif")), 
  xy = TRUE)

lims_m_chi <- rets_mcrw %>% 
  filter(((grepl("kappa_2500.5_chi_0.8.tif", landscape) |
             landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif") 
          & (id == 1 & habitat_effect == 0.25)) |
           (landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" 
            & id == 2 & habitat_effect == 0.25)) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_m_chi <- data.frame(diffx = max(c(lims_m_chi$diffx, lims_m_chi$diffy)),
                          diffy = max(c(lims_m_chi$diffx, lims_m_chi$diffy)))

lims_m_chi1 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" &
           habitat_effect == 0.25 &
           id == 2) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_m_chi$diffx/2 + 5000,
         limy_low = center_y - lims_m_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_m_chi$diffy/2 + 5000)

et_m_chi1 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" &
           habitat_effect == 0.25 &
           id == 2) %>% 
  ggplot() +
  geom_raster(data = land_c1, mapping = aes(x = x, y = y, fill = land_c1[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_chi1$limx_low, lims_m_chi1$limx_hig), 
                  ylim = c(lims_m_chi1$limy_low, lims_m_chi1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_chi1$limx_hig, y = lims_m_chi1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 0.1", parse = TRUE, fill = "white", alpha = 0.8)

lims_m_chi2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_m_chi$diffx/2 + 5000,
         limy_low = center_y - lims_m_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_m_chi$diffy/2 + 5000)

et_m_chi2 <- rets_mcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_c2, mapping = aes(x = x, y = y, fill = land_c2[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_chi2$limx_low, lims_m_chi2$limx_hig), 
                  ylim = c(lims_m_chi2$limy_low, lims_m_chi2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_chi2$limx_hig, y = lims_m_chi2$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 0.8", parse = TRUE, fill = "white", alpha = 0.8)

lims_m_chi3 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_m_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_m_chi$diffx/2 + 5000,
         limy_low = center_y - lims_m_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_m_chi$diffy/2 + 5000)

et_m_chi3 <- rets_mcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif" &
           habitat_effect == 0.25 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_c3, mapping = aes(x = x, y = y, fill = land_c3[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_m_chi3$limx_low, lims_m_chi3$limx_hig), 
                  ylim = c(lims_m_chi3$limy_low, lims_m_chi3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_m_chi3$limx_hig, y = lims_m_chi3$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 1.5", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_chi <- rets_bcrw %>% 
  filter(((grepl("kappa_2500.5_chi_0.8.tif", landscape) |
             landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif") 
          & (id == 1 & beta == 0.15)) |
           (landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" 
            & id == 2 & beta == 0.15)) %>% 
  group_by(scenario_id) %>% 
  summarize(maxx = max(x), maxy = max(y), 
            minx = min(x), miny = min(y),
            diffx = maxx - minx,
            diffy = maxy - miny)

lims_b_chi <- data.frame(diffx = max(c(lims_b_chi$diffx, lims_b_chi$diffy)),
                          diffy = max(c(lims_b_chi$diffx, lims_b_chi$diffy)))

lims_b_chi1 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" &
           beta == 0.15 &
           id == 2) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_b_chi$diffx/2 + 5000,
         limy_low = center_y - lims_b_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_b_chi$diffy/2 + 5000)

et_b_chi1 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_0.1.tif" &
           beta == 0.15 &
           id == 2) %>% 
  ggplot() +
  geom_raster(data = land_c1, mapping = aes(x = x, y = y, fill = land_c1[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[1]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_chi1$limx_low, lims_b_chi1$limx_hig), 
                  ylim = c(lims_b_chi1$limy_low, lims_b_chi1$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label",
           x = lims_b_chi1$limx_hig, y = lims_b_chi1$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 0.1", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_chi2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_b_chi$diffx/2 + 5000,
         limy_low = center_y - lims_b_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_b_chi$diffy/2 + 5000)

et_b_chi2 <- rets_bcrw %>% 
  filter(grepl("kappa_2500.5_chi_0.8.tif", landscape) &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_c2, mapping = aes(x = x, y = y, fill = land_c2[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[2]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_chi2$limx_low, lims_b_chi2$limx_hig), 
                  ylim = c(lims_b_chi2$limy_low, lims_b_chi2$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b_chi2$limx_hig, y = lims_b_chi2$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 0.8", parse = TRUE, fill = "white", alpha = 0.8)

lims_b_chi3 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif" &
           beta == 0.15 &
           id == 1) %>% 
  summarize(center_x = mean(x), center_y = mean(y)) %>% 
  mutate(limx_low = center_x - lims_b_chi$diffx/2 - 5000,
         limx_hig = center_x + lims_b_chi$diffx/2 + 5000,
         limy_low = center_y - lims_b_chi$diffy/2 - 5000,
         limy_hig = center_y + lims_b_chi$diffy/2 + 5000)

et_b_chi3 <- rets_bcrw %>% 
  filter(landscape == "output/Landscapes/land_kappa_2500.5_chi_1.5.tif" &
           beta == 0.15 &
           id == 1) %>% 
  ggplot() +
  geom_raster(data = land_c3, mapping = aes(x = x, y = y, fill = land_c3[, 3])) +
  geom_path(mapping = aes(x = x, y = y), color = ggcolors_landscape[3]) +
  scale_fill_gradient2(low = "white", high = 'black') +  
  coord_cartesian(xlim = c(lims_b_chi3$limx_low, lims_b_chi3$limx_hig), 
                  ylim = c(lims_b_chi3$limy_low, lims_b_chi3$limy_hig)) +
  theme_void() + 
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
  annotate("label", 
           x = lims_b_chi3$limx_hig, y = lims_b_chi3$limy_hig, 
           hjust = 1, vjust = 1, 
           label = "chi == 1.5", parse = TRUE, fill = "white", alpha = 0.8)

## Plots ####

fid_ref_1mo <- fid_ref %>% 
  filter(lag_words == "1 Month Lag")

f_crw_1mo <- f_crw %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "CRW")
f_ccrw_1mo <- f_ccrw %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "CCRW")
f_mcrw_1mo <- f_mcrw %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "MCRW")
f_bcrw_1mo <- f_bcrw %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "BCRW")

### CRW ####

crw_1mo <- ggplot(f_crw_1mo, aes(x = day, y = return_rate, 
                                 color = factor(rho), group = factor(rho))) +
  geom_line(alpha = 0.8) +
  labs(x = "Days", y = "Probability of return", color = expression(rho)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("CRW") +
  theme(legend.position = "bottom", text = element_text(size = 16))

crw_plot <- (crw_1mo | (et_crw1 / et_crw2 / et_crw3 / et_crw4)) +  plot_layout(widths = c(3, 1))

### CCRW ####

ccrw_1mo <- ggplot(f_ccrw_1mo, aes(x = day, y = return_rate, 
                                   color = factor(boundary_size), group = factor(boundary_size))) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = "r") +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("CCRW") +
  theme(legend.position = "bottom", text = element_text(size = 16))

ccrw_plot <- (ccrw_1mo | (et_ccrw1 / et_ccrw2 / et_ccrw3 / et_ccrw4)) +  plot_layout(widths = c(3, 1))

### MCRW ####

mcrw_1mo <- ggplot(f_mcrw_1mo, aes(x = day, y = return_rate, 
                                   color = factor(habitat_effect), group = factor(habitat_effect))) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(eta)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("MCRW") +
  theme(legend.position = "bottom", text = element_text(size = 16))

mcrw_plot <- (mcrw_1mo | (et_mcrw1 / et_mcrw2 / et_mcrw3 / et_mcrw4)) +  plot_layout(widths = c(3, 1))

### BCRW ####
 
bcrw_1mo <- ggplot(f_bcrw_1mo, aes(x = day, y = return_rate, 
                                   color = factor(beta), group = factor(beta))) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(beta)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("BCRW") +
  theme(legend.position = "bottom", text = element_text(size = 16))

bcrw_plot <- (bcrw_1mo | (et_bcrw1 / et_bcrw2 / et_bcrw3 / et_bcrw4)) +  plot_layout(widths = c(3, 1))

### MCRW + chi ####

f_mcrw_chi_1mo <- f_mcrw_chi %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "MCRW")

mcrw_chi_plot <- ggplot(f_mcrw_chi_1mo, aes(x = day, y = return_rate, 
                                          color = chi, 
                                          group = chi)) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(kappa)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("MCRW") +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), text = element_text(size = 18))

(mcrw_chi_plot | (et_m_chi1 / et_m_chi2 / et_m_chi3)) + plot_layout(widths = c(2, 1))

### MCRW + kappa ####

f_mcrw_kappa_1mo <- f_mcrw_kappa %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "MCRW")

mcrw_kappa_plot <- ggplot(f_mcrw_kappa_1mo, aes(x = day, y = return_rate, 
                                          color = kappa, 
                                          group = kappa)) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(chi)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("MCRW") +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), text = element_text(size = 18))

(mcrw_kappa_plot | (et_m_kappa1 / et_m_kappa2 / et_m_kappa3)) + plot_layout(widths = c(2, 1))

### BCRW + chi ####

f_bcrw_chi_1mo <- f_bcrw_chi %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "BCRW")

bcrw_chi_plot <- ggplot(f_bcrw_chi_1mo, aes(x = day, y = return_rate, 
                                          color = chi, 
                                          group = chi)) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(kappa)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("BCRW") +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), text = element_text(size = 18))

(bcrw_chi_plot | (et_b_chi1 / et_b_chi2 / et_b_chi3)) + plot_layout(widths = c(2, 1))

### BCRW + kappa ####

f_bcrw_kappa_1mo <- f_bcrw_kappa %>% 
  filter(lag_words == "1 Month Lag") %>% 
  mutate(model = "BCRW")

bcrw_kappa_plot <- ggplot(f_bcrw_kappa_1mo, aes(x = day, y = return_rate, 
                                          color = kappa, 
                                          group = kappa)) +
  geom_line(alpha = 0.8) +
  geom_smooth(data = fid_ref_1mo, method = "glm", formula = "y ~ log(x)",
              color = "gray50") + 
  labs(x = "Days", y = "Probability of return", color = expression(chi)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle("BCRW") +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_blank(), text = element_text(size = 18))

(bcrw_kappa_plot | (et_b_kappa1 / et_b_kappa2 / et_b_kappa3)) + plot_layout(widths = c(2, 1))


## Combined plots (MS figures) ####

### Figure 2 ####

# CRW
(crw_1mo + (et_crw1 / et_crw2 / et_crw3 / et_crw4) +
    # CCRW
    ccrw_1mo + (et_ccrw1 / et_ccrw2 / et_ccrw3 / et_ccrw4) +
    # MCRW
    mcrw_1mo + (et_mcrw1 / et_mcrw2 / et_mcrw3 / et_mcrw4) +
    # BCRW
    bcrw_1mo + (et_bcrw1 / et_bcrw2 / et_bcrw3 / et_bcrw4)) +
  # Layout
  plot_layout(nrow = 2, ncol = 4, widths = c(2, 1, 2, 1)) 

### Figure 3 ####

### Figure 4 ####

(mcrw_kappa_plot + (et_m_kappa1 / et_m_kappa2 / et_m_kappa3) +
    mcrw_chi_plot + (et_m_chi1 / et_m_chi2 / et_m_chi3) +
    bcrw_kappa_plot + (et_b_kappa1 / et_b_kappa2 / et_b_kappa3) +
    bcrw_chi_plot + (et_b_chi1 / et_b_chi2 / et_b_chi3)) +
  # Layout
  plot_layout(nrow = 2, ncol = 4, widths = c(2, 1, 2, 1)) 

### Figure 5 ####

ggplot(fid_ref, aes(x = day, y = return_rate)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ lag_words, nrow = 1, scales = "free_x") +
  labs(x = "Days", y = "Probability of return", color = expression(rho)) +
  theme_bw() +
  ylim(0, 1) +
  ggtitle(" ") +
  theme(legend.position = "bottom", text = element_text(size = 16))

### Figure 6 ####
