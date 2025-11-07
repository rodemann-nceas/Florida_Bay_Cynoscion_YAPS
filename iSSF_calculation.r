##iSSF with HMM data - Let's try it!
library(amt)
library(tidyverse)
library(survival)
library(terra)
library(sf)
library(furrr)
library(purrr)
library(ggeffects)

#Start with Rankin
rdat <- read.csv('Results/Rankin_full_hmm.csv')
rmap <- rast('Data/Rankin_map.tif')
str(rdat)

?make_track

fish_trackr <- rdat %>%
  mutate(ID = as.factor(ID), 
         class = as.factor(class),
         ts = as.POSIXct(ts, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")) %>% 
  make_track(
    .x = x, # X coordinate column
    .y = y, # Y coordinate column
    .t = ts, # Timestamp column
    id = ID, # Individual ID column
    crs = st_crs(rmap),
    state = state,
    pland = pland,
    ed = ed,
    pd = pd,
    shape_mn = shape_mn,
    class = class
  )

# Convert the track to observed steps, grouped by track
# This calculates step length (sl_) and turn angle (ta_) for the observed steps.
steps_obsr <- fish_trackr %>%
  steps(keep_cols = 'end') %>% 
  dplyr::filter(!is.na(ta_))

# Generate random steps (10 random steps for every 1 observed step)
steps_randr <- steps_obsr %>%
  random_steps(n_control = 10)

#We have our random steps, now we need to extract the habitat variables at the end of each step
head(steps_randr)

sr <- steps_randr %>% drop_na() %>% dplyr::filter(case_ == FALSE) %>% dplyr::select(-c(pland, ed, pd, shape_mn)) %>% mutate(ID = 1:n())
st <- steps_randr %>% drop_na() %>% dplyr::filter(case_ == TRUE)
#create a spatialpointsdataframe and extract from map
sr2 <- st_as_sf(sr, coords = c("x2_", "y2_"), crs = st_crs(rmap)) %>% terra::extract(x = rmap, y = ., ID = F) %>% bind_cols(sr) %>% 
mutate(class = if_else(savclass_2 == 1,2,1)) %>% drop_na() %>% dplyr::select(-savclass_2)

#run landscape metrics on points
rcoord <- as.data.frame(cbind(sr2$x2_, sr2$y2_))
rcoord_sf <- st_as_sf(rcoord, coords = c("V1", "V2"), crs = st_crs(rmap))

plan(multisession)  # Use multisession to run on multiple cores


result_metr <- 1:nrow(rcoord_sf) %>% future_map(
  function(i) {
    rmap_r <- rast("Data/Rankin_map.tif")
    sample_lsm(
      rmap_r,
      y = st_geometry(rcoord_sf[i, ]),  # Ensure geometry is passed
      size = 5,
      what = c(
        "lsm_c_pd", "lsm_c_ed",
        "lsm_c_division",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)
plan(sequential)

rt <-lapply(result_metr, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})
metr <- bind_rows(rt, .id = 'source')

unique_sources <- unique(metr$source)
unique_classes <- unique(metr$class)
unique_metric <- unique(metr$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(metr, by = c("source", "class", 'metric')) %>% distinct()

metr1 <- mettry %>% pivot_wider(names_from = metric, values_from = value) %>% dplyr::filter(class == 1)
names(metr1)

rankrand <- metr1 %>% dplyr::select(c(division, ed, pd, pland, shape_mn)) %>% bind_cols(sr2) %>% dplyr::select(-c(division,ID))

strue <- rbind(st,rankrand)
names(st)
names(rankrand)
# Standardize continuous habitat variables (z-score transformation)
# This helps with model convergence and interpretation of interaction terms.
steps_finalr <- strue %>%
  mutate(
    # Standardize continuous habitat covariates
    pland_z = scale(pland)[, 1],
    ed_z = scale(ed)[, 1],
    pd_z = scale(pd)[, 1],
    shape_mn_z = scale(shape_mn)[, 1],
    
    # Convert state to a factor for proper modeling interaction
    state_f = factor(state)
  ) %>% drop_na()
str(steps_finalr)
# Fit the iSSF using conditional logistic regression (clogit).
issf_modelr <- steps_finalr %>%
  drop_na() %>% 
  fit_clogit(case_ ~
               # Movement terms 
               sl_ + I(sl_^2) +
               cos(ta_) +
               # Habitat selection with interaction with HMM state
               pland_z * state_f + 
               ed_z * state_f +
               pd_z * state_f +
               shape_mn_z * state_f +
               class * state_f +
               # The required stratification term for iSSF
               strata(step_id_), model = T
  )

summary(issf_modelr)

#Plotting!
# --- Define Constants (Assuming these are still in your R environment) ---
SL_REF <- 6.15 
SL_RANGE <- seq(0.1, 12.2, length.out = 100)
Z_RANGE <- seq(-2, 2, length.out = 100)
N_POINTS <- length(SL_RANGE) # 100
N_POINTS_Z <- length(Z_RANGE)
ALL_HABITAT_VARS <- c("pland_z", "ed_z", "pd_z", "shape_mn_z")

# -----------------------------------------------------------------
# --- 1. Step Length Selection Curve (sl_ and I(sl_^2)) ---
# --- (Needs to be run first to define rss_sl) ---
# -----------------------------------------------------------------

# --- x1: Test data (varying sl_) - 100 rows ---
p_sl_x1 <- data.frame(
  sl_ = SL_RANGE, 
  ta_ = rep(0, N_POINTS),             
  pland_z = rep(0, N_POINTS), ed_z = rep(0, N_POINTS), pd_z = rep(0, N_POINTS), shape_mn_z = rep(0, N_POINTS), 
  state_f = factor(rep("1", N_POINTS), levels = c("1", "2")), 
  class = factor(rep("1", N_POINTS), levels = c("1", "2"))  
)

# --- x2: Reference data (fixed at median sl_) - 1 row ---
p_sl_x2 <- data.frame(
  sl_ = SL_REF, 
  ta_ = 0, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0, 
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

# Calculate log-RSS for Step Length
rss_sl <- log_rss(issf_modelr, x1 = p_sl_x1, x2 = p_sl_x2, ci = 'se')

# Plot the Step Length log-RSS curve
plot_sl <- ggplot(rss_sl$df, aes(x = sl_x1, y = log_rss)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "Log-RSS for Step Length (zero turn angle)",
    x = "Step Length (meters)",
    y = "Log Relative Selection Strength"
  )
print(plot_sl)

#habitat selection
#' Calculates and plots the difference in log-RSS (State 2 - State 1) for a given habitat variable
plot_interaction_effect <- function(var_name, title_name) {
  
  N_POINTS_HAB <- 2 * N_POINTS_Z # 200 points total
  
  # --- 1. Create x1 (Test Data: State 1 and State 2) ---
  # This uses the stable 200-row method from the previous successful function
  p_hab_x1 <- data.frame(
    VARYING_VAR = rep(Z_RANGE, 2), 
    state_f = factor(rep(c("1", "2"), each = N_POINTS_Z), levels = c("1", "2"))
  )
  names(p_hab_x1)[1] <- var_name
  
  p_const <- data.frame(
    sl_ = rep(SL_REF, N_POINTS_HAB), 
    ta_ = rep(0, N_POINTS_HAB),
    class = factor(rep("1", N_POINTS_HAB), levels = c("1", "2"))
  )
  
  fixed_hab_vars <- setdiff(ALL_HABITAT_VARS, var_name)
  for (h_var in fixed_hab_vars) {
    p_const[[h_var]] <- rep(0, N_POINTS_HAB)
  }
  
  p_hab_x1 <- bind_cols(p_hab_x1, p_const)
  
  # --- 2. Create x2 (Reference Data: 1 row, fixed at mean/neutral value in State 1) ---
  p_hab_x2 <- data.frame(
    sl_ = SL_REF, ta_ = 0, 
    class = factor("1", levels = c("1", "2")),
    state_f = factor("1", levels = c("1", "2"))
  )
  for (h_var in ALL_HABITAT_VARS) {
    p_hab_x2[[h_var]] <- 0
  }
  
  # --- 3. Calculate log-RSS for both states relative to the baseline ---
  rss_hab <- log_rss(issf_modelr, x1 = p_hab_x1, x2 = p_hab_x2, ci = 'se')
  df <- rss_hab$df
  
  df_state1 <- df %>% filter(state_f_x1 == "1") %>% arrange(get(paste0(var_name, "_x1")))
  df_state2 <- df %>% filter(state_f_x1 == "2") %>% arrange(get(paste0(var_name, "_x1")))
  
  # Calculate SE of the difference (using a conservative approximation, ignoring covariance)
  # SE_diff = sqrt(SE_1^2 + SE_2^2)
  se_diff <- sqrt(((df_state1$upr-df_state1$lwr)/3.92)^2 + ((df_state2$upr-df_state2$lwr)/3.92)^2)
  
  # Create a new data frame for the difference plot
  df_diff <- data.frame(
    z_score = df_state1[[paste0(var_name, "_x1")]],
    log_rss_diff = df_state2$log_rss - df_state1$log_rss,
    # Calculate 95% CI bounds (Z-score 1.96 for 95%)
    lwr_diff = (df_state2$log_rss - df_state1$log_rss) - 1.96 * se_diff,
    upr_diff = (df_state2$log_rss - df_state1$log_rss) + 1.96 * se_diff
  )
  
  # --- 5. Plot the Difference with CI ---
  ggplot(df_diff, aes(x = z_score, y = log_rss_diff)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    # Add the confidence interval ribbon
    geom_ribbon(aes(ymin = lwr_diff, ymax = upr_diff), alpha = 0.2, fill = "darkgreen") +
    theme_bw() +
    labs(
      title = paste("State-Dependent Effect: Difference in", title_name, "Selection"),
      subtitle = "log-RSS(State 2) - log-RSS(State 1). CI calculated via SE propagation.",
      x = paste0(title_name, " (Standardized Z-Score)"),
      y = "log-RSS Difference (State 2 - State 1)"
    )
}

# -----------------------------------------------------------------
# --- Generate Interaction Effect Plots for All Habitat Variables ---
# -----------------------------------------------------------------

#PLAND
plot_pland_interaction <- plot_interaction_effect("pland_z", "PLAND")
print(plot_pland_interaction)

#PD
plot_pd_interaction <- plot_interaction_effect("pd_z", "Patch Density")
print(plot_pd_interaction)

#Shape_mn
plot_shape_interaction <- plot_interaction_effect("shape_mn_z", "Mean Shape")
print(plot_shape_interaction)

#amazing!!!! Let's save these plots and run for Johnson
ggsave(filename = 'Results/iSSF_sl_Rankin.png', plot_sl, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_pland_HMM_Rankin.png', plot_pland_interaction, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_pd_HMM_Rankin.png', plot_pd_interaction, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_shape_HMM_Rankin.png', plot_shape_interaction, dpi = 300, height = 5, width = 5)


####Johnson
jdat <- read.csv('Results/Johnson_full_hmm.csv')
jmap <- rast('Data/Johnson_map.tif')
str(jdat)

?make_track

fish_trackj <- jdat %>%
  mutate(ID = as.factor(ID), 
         class = as.factor(class),
         ts = as.POSIXct(ts, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")) %>% 
  make_track(
    .x = x, # X coordinate column
    .y = y, # Y coordinate column
    .t = ts, # Timestamp column
    id = ID, # Individual ID column
    crs = st_crs(jmap),
    state = state,
    pland = pland,
    ed = ed,
    pd = pd,
    shape_mn = shape_mn,
    class = class
  )

# Convert the track to observed steps, grouped by track
# This calculates step length (sl_) and turn angle (ta_) for the observed steps.
steps_obsj <- fish_trackj %>%
  steps(keep_cols = 'end') %>% 
  dplyr::filter(!is.na(ta_))

# Generate random steps (10 random steps for every 1 observed step)
steps_randj <- steps_obsj %>%
  random_steps(n_control = 10)

#We have our random steps, now we need to extract the habitat variables at the end of each step
head(steps_randr)

sj <- steps_randj %>% drop_na() %>% dplyr::filter(case_ == FALSE) %>% dplyr::select(-c(pland, ed, pd, shape_mn)) %>% mutate(ID = 1:n())
stj <- steps_randj %>% drop_na() %>% dplyr::filter(case_ == TRUE)
#create a spatialpointsdataframe and extract from map
sj2 <- st_as_sf(sj, coords = c("x2_", "y2_"), crs = st_crs(jmap)) %>% terra::extract(x = jmap, y = ., ID = F) %>% bind_cols(sj) %>% 
  mutate(class = if_else(johnclass2 == 1,2,1)) %>% drop_na() %>% dplyr::select(-johnclass2)

#run landscape metrics on points
jcoord <- as.data.frame(cbind(sj2$x2_, sj2$y2_))
jcoord_sf <- st_as_sf(jcoord, coords = c("V1", "V2"), crs = st_crs(jmap))

plan(multisession)  # Use multisession to run on multiple cores


result_metj <- 1:nrow(jcoord_sf) %>% future_map(
  function(i) {
    jmap_r <- rast("Data/Johnson_map.tif")
    sample_lsm(
      jmap_r,
      y = st_geometry(jcoord_sf[i, ]),  # Ensure geometry is passed
      size = 15,
      what = c(
        "lsm_c_pd", "lsm_c_ed",
        "lsm_c_division",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)
plan(sequential)

rt <-lapply(result_metj, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})
metj <- bind_rows(rt, .id = 'source')

unique_sources <- unique(metj$source)
unique_classes <- unique(metj$class)
unique_metric <- unique(metj$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(metj, by = c("source", "class", 'metric')) %>% distinct()

metj1 <- mettry %>% pivot_wider(names_from = metric, values_from = value) %>% dplyr::filter(class == 1)
names(metr1)

johnrand <- metj1 %>% dplyr::select(c(division, ed, pd, pland, shape_mn)) %>% bind_cols(sj2) %>% dplyr::select(-c(division,ID))

sjtrue <- rbind(stj,johnrand)

# Standardize continuous habitat variables (z-score transformation)
# This helps with model convergence and interpretation of interaction terms.
steps_finalj <- sjtrue %>%
  mutate(
    # Standardize continuous habitat covariates
    pland_z = scale(pland)[, 1],
    ed_z = scale(ed)[, 1],
    pd_z = scale(pd)[, 1],
    shape_mn_z = scale(shape_mn)[, 1],
    
    # Convert state to a factor for proper modeling interaction
    state_f = factor(state)
  ) %>% drop_na()

write.csv(steps_finalj, file = 'Data/iSSF_Johnson.csv')
# Fit the iSSF using conditional logistic regression (clogit).
issf_modelj <- steps_finalj %>%
  drop_na() %>% 
  fit_clogit(case_ ~
               # Movement terms 
               sl_ + I(sl_^2) +
               cos(ta_) +
               # Habitat selection with interaction with HMM state
               pland_z * state_f + 
               ed_z * state_f +
               pd_z * state_f +
               shape_mn_z * state_f +
               class * state_f +
               # The required stratification term for iSSF
               strata(step_id_), model = T
  )

summary(issf_modelj)
max(steps_finalj$sl_)
#Plotting!
SL_REF <- 25.1 # Example: Replace with the actual mean step length from your data
SL_RANGE <- seq(0, 380, length.out = 100) # Example: Replace with the actual min/max step length range
Z_RANGE <- seq(-2, 2, length.out = 100) # Standard z-score range
ALL_HABITAT_VARS <- c("pland_z", "ed_z", "pd_z", "shape_mn_z")
CLASS_REF <- factor("1", levels = c("1", "2")) 
STATE_REF <- factor("1", levels = c("1", "2"))

#plot for step length
p_sl_x1 <- data.frame(
  sl_ = SL_RANGE, 
  ta_ = rep(0, N_POINTS),             
  pland_z = rep(0, N_POINTS), ed_z = rep(0, N_POINTS), pd_z = rep(0, N_POINTS), shape_mn_z = rep(0, N_POINTS), 
  state_f = factor(rep("1", N_POINTS), levels = c("1", "2")), 
  class = factor(rep("1", N_POINTS), levels = c("1", "2"))  
)

# --- x2: Reference data (fixed at median sl_) - 1 row ---
p_sl_x2 <- data.frame(
  sl_ = SL_REF, 
  ta_ = 0, 
  pland_z = 0, ed_z = 0, pd_z = 0, shape_mn_z = 0, 
  state_f = factor("1", levels = c("1", "2")),
  class = factor("1", levels = c("1", "2"))
)

# Calculate log-RSS for Step Length
rss_sl <- log_rss(issf_modelj, x1 = p_sl_x1, x2 = p_sl_x2, ci = 'se')

# Plot the Step Length log-RSS curve
plot_slj <- ggplot(rss_sl$df, aes(x = sl_x1, y = log_rss)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
  theme_bw() +
  labs(
    title = "Log-RSS for Step Length (zero turn angle)",
    x = "Step Length (meters)",
    y = "Log Relative Selection Strength"
  )
print(plot_slj)

#plot for turn angle
plot_ta_log_rss <- function() {
  # Range of turning angle in radians (-pi to pi)
  TA_RANGE <- seq(-pi, pi, length.out = 100)
  N_POINTS <- length(TA_RANGE)
  
  # --- x1: Test data (varying ta_, fixed sl_=SL_REF, state_f=1) ---
  p_ta_x1 <- data.frame(
    sl_ = rep(SL_REF, N_POINTS), 
    ta_ = TA_RANGE, 
    state_f = rep(STATE_REF, N_POINTS),
    class = rep(CLASS_REF, N_POINTS)
  )
  
  # Add ALL habitat variables fixed at Z=0
  for (h_var in ALL_HABITAT_VARS) {
    p_ta_x1[[h_var]] <- rep(0, N_POINTS)
  }
  
  # --- x2: Reference data (1 row, fixed at mean/neutral value) ---
  # Reference: Mean step length (SL_REF), maximum turning angle (ta_=pi), State 1, Class 1, Z=0 for habitat
  p_ta_x2 <- data.frame(
    sl_ = SL_REF, ta_ = pi, # Reference is max turn (cos(ta) = -1)
    state_f = STATE_REF,
    class = CLASS_REF
  )
  for (h_var in ALL_HABITAT_VARS) {
    p_ta_x2[[h_var]] <- 0
  }
  
  rss_ta <- log_rss(issf_modelj, x1 = p_ta_x1, x2 = p_ta_x2, ci = 'se')
  
  # Plot
  ggplot(rss_ta$df, aes(x = ta_x1, y = log_rss)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2, color = "darkblue") +
    # Add the CI ribbon (assuming lwr/upr are present)
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkblue") +
    theme_bw() +
    scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                       labels = c("-π (max turn)", "-π/2", "0 (straight)", "π/2", "π (max turn)")) +
    labs(
      title = "Log-RSS for Turning Angle (mean step lenth)",
      x = "Turning Angle (Radians)",
      y = "Log Relative Selection Strength"
    )
}

#plot for significant singular habitat variables
plot_habitat_log_rss <- function(var_name, title_name) {
  
  N_POINTS <- length(Z_RANGE)
  
  # --- x1: Test data (varying habitat var, fixed sl_=SL_REF, state_f=1) ---
  p_hab_x1 <- data.frame(
    VARYING_VAR = Z_RANGE, 
    sl_ = rep(SL_REF, N_POINTS), 
    ta_ = rep(0, N_POINTS), 
    state_f = rep(STATE_REF, N_POINTS),
    class = rep(CLASS_REF, N_POINTS)
  )
  names(p_hab_x1)[1] <- var_name
  
  # Add OTHER habitat variables fixed at Z=0
  fixed_hab_vars <- setdiff(ALL_HABITAT_VARS, var_name)
  for (h_var in fixed_hab_vars) {
    p_hab_x1[[h_var]] <- rep(0, N_POINTS)
  }
  
  # --- x2: Reference data (1 row, fixed at Z=0) ---
  # Reference: Mean sl, straight ta, State 1, Class 1, Z=0 for ALL habitat vars
  p_hab_x2 <- data.frame(
    sl_ = SL_REF, ta_ = 0, 
    state_f = STATE_REF,
    class = CLASS_REF
  )
  for (h_var in ALL_HABITAT_VARS) {
    p_hab_x2[[h_var]] <- 0
  }
  
  rss_hab <- log_rss(issf_modelj, x1 = p_hab_x1, x2 = p_hab_x2, ci = 'se')
  
  # Plot
  ggplot(rss_hab$df, aes(x = get(paste0(var_name, "_x1")), y = log_rss)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    # Add the CI ribbon (assuming lwr/upr are present)
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkgreen") +
    theme_bw() +
    labs(
      title = paste("Log-RSS for", title_name),
      x = paste0(title_name, " (Standardized Z-Score)"),
      y = "Log Relative Selection Strength"
    )
}

#turn angle plot
plot_ta <- plot_ta_log_rss()
plot_ta

#PLAND
plot_plandj <- plot_habitat_log_rss("pland_z", "PLAND")
print(plot_plandj)

# 4. Edge Density (ed_z)
plot_edj <- plot_habitat_log_rss("ed_z", "Edge Density")
print(plot_edj)

# 5. Patch Density (pd_z)
plot_pdj <- plot_habitat_log_rss("pd_z", 'Patch Density')
print(plot_pdj)

# 6. Mean Patch Shape (shape_mn_z)
plot_shape_mnj <- plot_habitat_log_rss("shape_mn_z", "Mean Patch Shape")
print(plot_shape_mnj)

#interactions between states
plot_interaction_effect <- function(var_name, title_name) {
  
  N_POINTS_HAB <- 2 * N_POINTS_Z # 200 points total
  
  # --- 1. Create x1 (Test Data: State 1 and State 2) ---
  # This uses the stable 200-row method from the previous successful function
  p_hab_x1 <- data.frame(
    VARYING_VAR = rep(Z_RANGE, 2), 
    state_f = factor(rep(c("1", "2"), each = N_POINTS_Z), levels = c("1", "2"))
  )
  names(p_hab_x1)[1] <- var_name
  
  p_const <- data.frame(
    sl_ = rep(SL_REF, N_POINTS_HAB), 
    ta_ = rep(0, N_POINTS_HAB),
    class = factor(rep("1", N_POINTS_HAB), levels = c("1", "2"))
  )
  
  fixed_hab_vars <- setdiff(ALL_HABITAT_VARS, var_name)
  for (h_var in fixed_hab_vars) {
    p_const[[h_var]] <- rep(0, N_POINTS_HAB)
  }
  
  p_hab_x1 <- bind_cols(p_hab_x1, p_const)
  
  # --- 2. Create x2 (Reference Data: 1 row, fixed at mean/neutral value in State 1) ---
  p_hab_x2 <- data.frame(
    sl_ = SL_REF, ta_ = 0, 
    class = factor("1", levels = c("1", "2")),
    state_f = factor("1", levels = c("1", "2"))
  )
  for (h_var in ALL_HABITAT_VARS) {
    p_hab_x2[[h_var]] <- 0
  }
  
  # --- 3. Calculate log-RSS for both states relative to the baseline ---
  rss_hab <- log_rss(issf_modelj, x1 = p_hab_x1, x2 = p_hab_x2, ci = 'se')
  df <- rss_hab$df
  
  df_state1 <- df %>% filter(state_f_x1 == "1") %>% arrange(get(paste0(var_name, "_x1")))
  df_state2 <- df %>% filter(state_f_x1 == "2") %>% arrange(get(paste0(var_name, "_x1")))
  
  # Calculate SE of the difference (using a conservative approximation, ignoring covariance)
  # SE_diff = sqrt(SE_1^2 + SE_2^2)
  se_diff <- sqrt(((df_state1$upr-df_state1$lwr)/3.92)^2 + ((df_state2$upr-df_state2$lwr)/3.92)^2)
  
  # Create a new data frame for the difference plot
  df_diff <- data.frame(
    z_score = df_state1[[paste0(var_name, "_x1")]],
    log_rss_diff = df_state2$log_rss - df_state1$log_rss,
    # Calculate 95% CI bounds (Z-score 1.96 for 95%)
    lwr_diff = (df_state2$log_rss - df_state1$log_rss) - 1.96 * se_diff,
    upr_diff = (df_state2$log_rss - df_state1$log_rss) + 1.96 * se_diff
  )
  
  # --- 5. Plot the Difference with CI ---
  ggplot(df_diff, aes(x = z_score, y = log_rss_diff)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_line(linewidth = 1.2, color = "darkgreen") +
    # Add the confidence interval ribbon
    geom_ribbon(aes(ymin = lwr_diff, ymax = upr_diff), alpha = 0.2, fill = "darkgreen") +
    theme_bw() +
    labs(
      title = paste("State-Dependent Effect: Difference in", title_name, "Selection"),
      x = paste0(title_name, " (Standardized Z-Score)"),
      y = "log-RSS Difference (State 2 - State 1)"
    )
}

#PD
plot_pd_interactionj <- plot_interaction_effect("pd_z", "Patch Density")
print(plot_pd_interactionj)

#ED
plot_ed_interactionj <- plot_interaction_effect("ed_z", "Edge Density")
print(plot_ed_interactionj)

ggsave(filename = 'Results/iSSF_PD_HMM_Johnson.png', plot_pd_interactionj, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_ed_HMM_Johnson.png', plot_ed_interactionj, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_shape_Johnson.png', plot_shape_mnj, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_pland_Johnson.png', plot_plandj, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_pd_Johnson.png', plot_pdj, dpi = 300, height = 5, width = 5)
ggsave(filename = 'Results/iSSF_ed_Johnson.png', plot_edj, dpi = 300, height = 5, width = 5)
