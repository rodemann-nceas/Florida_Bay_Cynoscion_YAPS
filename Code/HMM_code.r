###Hidden Markov Models####
###Author: Jon Rodemann
#Load in Libraries

library(tidyverse)
library(lubridate)
library(moveHMM)
library(purrr)
library(readr)
library(parallel)
library(circular)

##Rankin YAPS####
#1) load in data

dat <- read_csv("Data/Rankin_track_dat.csv", show_col_types = FALSE)

glimpse(dat)

# The script expects at least these columns:
#   tag          : individual identifier
#   track        : track number/ID for that individual's track segment
#   ts           : timestamp string (ISO-like)
#   x, y         : coordinates (assumed UTM/metres). If they're lat/long, convert to projected meters first.

# Parse timestamp
dat <- dat %>%
  mutate(ts = as.POSIXct(.data[[ts_col]],
                         format = "%Y-%m-%d %H:%M:%OS",
                         tz = "UTC"))

#2) Interpolate to regular 3-minute intervals by individual and track
interpolate_track <- function(df_group, dt_seconds = 180) {
  # df_group must contain columns: ts, x, y, plus id/track identifiers
  if (nrow(df_group) < 2) return(NULL)  # can't interpolate single-point tracks
  
  tmin <- min(df_group$ts, na.rm = TRUE)
  tmax <- max(df_group$ts, na.rm = TRUE)
  
  new_times <- seq(from = tmin, to = tmax, by = dt_seconds)
  if (length(new_times) < 2) return(NULL)
  
  # numeric times for approx
  orig_t_num <- as.numeric(df_group$ts)
  new_t_num <- as.numeric(new_times)
  
  # ensure sorted by time
  ord <- order(orig_t_num)
  orig_t_num <- orig_t_num[ord]
  x_orig <- df_group[[x_col]][ord]
  y_orig <- df_group[[y_col]][ord]
  
  # linear interpolation (NA handling: approx requires finite values)
  ok_x <- is.finite(x_orig) & is.finite(orig_t_num)
  ok_y <- is.finite(y_orig) & is.finite(orig_t_num)
  if (!any(ok_x) || !any(ok_y)) return(NULL)
  
  x_interp <- approx(x = orig_t_num[ok_x], y = x_orig[ok_x], xout = new_t_num, rule = 2)$y
  y_interp <- approx(x = orig_t_num[ok_y], y = y_orig[ok_y], xout = new_t_num, rule = 2)$y
  
  tibble(
    ts = new_times,
    !!id_col := df_group[[id_col]][1],
    !!track_col := df_group[[track_col]][1],
    x = x_interp,
    y = y_interp
  )
}

head(dat)
unique(dat$tag)
#group by 20-40 second tags and 60-90 second tags 
dat_20 <- dat %>% dplyr::filter(tag %in% c(34497,34499,34501,34502,34505,34507,34508,34509,34510,34513,34514,34518))
dat_60 <- dat %>% dplyr::filter(tag %in% c(854,871,875,880,881,904,912,47421))

##First look at YAPS tags (20-40 second random delay)
# Apply interpolation grouped by ID and track per time delay of tag
interp_list_20 <- dat_20 %>%
  group_by(.data[[id_col]], .data[[track_col]]) %>%
  group_map(~ interpolate_track(.x, dt_seconds = 30), .keep = TRUE)

# Remove NULLs (tracks that couldn't be interpolated)
interp_list_20 <- keep(interp_list_20, ~ !is.null(.x) && nrow(.x) > 0)

interp_df_20 <- bind_rows(interp_list_20)
interp_df_20 <- interp_df_20 %>% mutate(ID = paste(tag, track, sep = "_"))


#3) Step length and turn angles for HMM
# moveHMM expects columns: ID, x, y (named exactly "ID" by default for prepData),
# and it will create 'step' and 'angle' (in radians), and 'x'/'y' must be in meters if you want step in meters.


# Run prepData to compute steps/angles on the new ID column
# The function prepData can accept a "type" argument; set type = "UTM" if coords are planar meters.
# moveHMM's prepData expects the data sorted by ID and time
interp_df_20 <- interp_df_20 %>% arrange(ID, ts)
interp_df_20 <- interp_df_20 %>% group_by(ID) %>% filter(n()>= 3) %>% ungroup
# prepData: coordinate names are "x" and "y" in our frame
hmm_data_20 <- prepData(interp_df_20, type = "UTM", coordNames = c("x", "y"))

# moveHMM uses a column 'step' and 'angle' (radians) created by prepData
# Remove rows with zero or NA step (stay/missing movement) because gamma requires positive steps
hmm_data_20 <- hmm_data_20 %>% filter(!is.na(step) & step > 0) %>% filter(!is.na(angle))

hmm_data_20 <- hmm_data_20 %>%
  filter(is.finite(step), is.finite(angle), step > 0) %>%
  mutate(
    step = pmax(step, quantile(step, 0.01, na.rm = TRUE)),
    step = pmin(step, quantile(step, 0.99, na.rm = TRUE))
  )


summary(hmm_data_20$step)
summary(hmm_data_20$angle)
hist(hmm_data_20$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x,0.8,rate=1), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,1.75,rate=1), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_20$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.5), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,7.5), n=200, add=TRUE, col="cadetblue", lwd=4)


#4) fit HMMs!
#look at histograms and get numbers
sl_init_mean <- c(1.75/1, 0.8/1)
sl_init_sd <- c(sqrt(1.75)/1, sqrt(0.8)/1)
ta_init_mean <- c(0, pi)
ta_init_con <- c(7.5, 0.5)

mod <- fitHMM(data = hmm_data_20,
              nbStates = 2,
              stepPar0 = c(sl_init_mean, sl_init_sd),
              anglePar0 = c(ta_init_mean, ta_init_con),
              formula = ~1,
              stepDist = "gamma",
              angleDist = "vm")

plot(mod)
CI(mod)
plot(mbest)
m$mod
states <- viterbi(mod)
?viterbi
states
data.frame(val = rle(states)$values, n = rle(states)$lengths) %>%
  ggplot(aes(val %>% factor, n)) + geom_violin()
mod %>% plotPR()



#we are going to try to fit a bunch of models, both with 2 and 3 states. We will do this with different starting values in parallel to maximize model
ncores <- detectCores()-1
cl <- makeCluster(getOption('cl.cores', ncores))
clusterExport(cl, list('hmm_data', 'fitHMM'))

#number of tries with different starting values
niter <- 25
set.seed(1919)

# Create list of starting values
allPar0 <- lapply(as.list(1:niter), function(x) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.5, 3),
                     max = c(1.5, 5))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.5, 3),
                   max = c(1.5, 5))
  # Turning angle mean
  angleMean0 <- c(pi, 0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.5, 5),
                     max = c(1.5, 10))
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0)
  anglePar0 <- c(angleMean0, angleCon0)
  return(list(step = stepPar0, angle = anglePar0))
})

# Fit the niter models in parallel
allm_parallel <- parLapply(cl = cl, X = allPar0, fun = function(par0) {
  m <- fitHMM(data = hmm_data, nbStates = 2, stepPar0 = par0$step,
              anglePar0 = par0$angle)
  return(m)
})

allnllk <- unlist(lapply(allm_parallel, function(m) m$mod$minimum))
allnllk

whichbest <- which.min(allnllk)
mbest <- allm_parallel[[whichbest]]
mbest

#5) Decode states and summarize metrics

# Viterbi sequence
viterbi_states <- viterbi(mbest)
# viterbi returns vector aligned to rows of the data in best_fit$data
hmm_data_with_state <- mbest$data %>%
  mutate(state = viterbi_states)

# Convert angle to degrees for interpretability (angles are in radians)
hmm_data_with_state <- hmm_data_with_state %>%
  mutate(angle_deg = angle * 180 / pi,
         abs_angle_deg = abs(angle_deg))

#5a) average metrics per state
state_summary <- hmm_data_with_state %>%
  group_by(state) %>%
  summarize(
    mean_step_m = mean(step, na.rm = TRUE),
    sd_step_m = sd(step, na.rm = TRUE),
    mean_abs_turn_deg = mean(abs_angle_deg, na.rm = TRUE),
    n_positions = n()
  ) %>%
  arrange(state)

#5b) average number of positions per track assigned to each state
# here 'burst' is the original track id; we compute, per burst, counts per state, then average across bursts
positions_per_track <- hmm_data_with_state %>%
  group_by(burst, state) %>%
  summarize(n_pos = n(), .groups = "drop")

avg_positions_per_track <- positions_per_track %>%
  group_by(state) %>%
  summarize(mean_positions_per_burst = mean(n_pos, na.rm = TRUE),
            median_positions_per_burst = median(n_pos, na.rm = TRUE),
            n_bursts = n())

# Combine summaries
state_metrics <- left_join(state_summary, avg_positions_per_track, by = "state")

#6) save results, print summaries

# add ID, burst, ts, x, y, step, angle, state to output CSV
output_df <- hmm_data_with_state %>%
  dplyr::select(ID, burst, ts, x, y, step, angle, angle_deg, state)

write_csv(output_df, "Results/interpolated_hmm_states.csv")
