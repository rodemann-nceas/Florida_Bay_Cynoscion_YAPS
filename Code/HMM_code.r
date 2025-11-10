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
#1) load in data for 20 seconds

dat <- read_csv('Data/Rankin_30s_hmmdat.csv', show_col_types = FALSE)

glimpse(dat)

ts_col <- "ts"
id_col <- "tag"
track_col <- "track"
x_col <- "x"
y_col <- "y"


# The script expects at least these columns:
#   tag          : individual identifier
#   track        : track number/ID for that individual's track segment
#   ts           : timestamp string (ISO-like)
#   x, y         : coordinates (assumed UTM/metres). If they're lat/long, convert to projected meters first.

# Parse timestamp
dat <- dat %>%
  mutate(ts = as.POSIXct(ts,
                         format = "%Y-%m-%d %H:%M:%OS",
                         tz = "UTC"))
str(dat)
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
dat %>% mutate(ID = paste(tag, track, sep = "_")) %>% dplyr::filter(ID == '34510_56') %>%
  ggplot()+
  geom_point(aes(x = x, y = y))+
  theme_classic()
d <- dat %>% mutate(ID = paste(tag, track, sep = "_")) %>% dplyr::filter(ID == '34510_56') %>% 
  mutate(tot = as.numeric(max(ts))-as.numeric(min(ts)))
##First look at YAPS tags (20-40 second random delay)
# Apply interpolation grouped by ID and track per time delay of tag
interp_list_20 <- dat %>%
  group_by(.data[[id_col]], .data[[track_col]]) %>%
  group_map(~ interpolate_track(.x, dt_seconds = 60), .keep = TRUE)

# Remove NULLs (tracks that couldn't be interpolated)
interp_list_20 <- keep(interp_list_20, ~ !is.null(.x) && nrow(.x) > 0)

interp_df_20 <- bind_rows(interp_list_20)
interp_df_20 <- interp_df_20 %>% mutate(ID = paste(tag, track, sep = "_"))
head(interp_df_20)

interp_df_20 %>% dplyr::filter(ID == '34510_56') %>% 
  ggplot()+
  geom_point(aes(x = x, y = y))+
  theme_classic()

id20 <- interp_df_20 %>% group_by(ID) %>% summarize(tott = max(ts)-min(ts))
id20$tottime <- as.numeric(id20$tott)
hist(id20$tottime/60, freq = F, breaks = 30)
#3) Step length and turn angles for HMM
# moveHMM expects columns: ID, x, y (named exactly "ID" by default for prepData),
# and it will create 'step' and 'angle' (in radians), and 'x'/'y' must be in meters if you want step in meters.


# Run prepData to compute steps/angles on the new ID column
# The function prepData can accept a "type" argument; set type = "UTM" if coords are planar meters.
# moveHMM's prepData expects the data sorted by ID and time
interp_df_20 <- interp_df_20 %>% arrange(ID, ts)
interp_df_20 <- interp_df_20 %>% group_by(ID) %>% filter(n()>= 20) %>% ungroup
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
curve(dgamma(x,0.75,rate=0.9), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,2,rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_20$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.01), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,7), n=200, add=TRUE, col="cadetblue", lwd=4)


#4) fit HMMs!
#look at histograms and get numbers
sl_init_mean <- c(2/0.6, 0.75/0.9)
sl_init_sd <- c(sqrt(2)/0.6, sqrt(0.75)/0.9)
ta_init_mean <- c(0, pi)
ta_init_con <- c(7, 0.01)

mod <- fitHMM(data = hmm_data_20,
              nbStates = 2,
              stepPar0 = c(sl_init_mean, sl_init_sd),
              anglePar0 = c(ta_init_mean, ta_init_con),
              formula = ~1)
unique(hmm_data_20$ID)
mod
plot(mod)

#remove tracks that do not look good
hmm_data_20_2 <- hmm_data_20 %>% dplyr::filter(!ID %in% c('34501_39','34501_40','34501_44','34501_51','34501_7','34507_1','34507_14','34507_15',
                                                       '34507_27','34507_29','34507_3','34507_33','34507_40','34507_43','34507_57','34507_59',
                                                       '34507_61','34507_64','34507_66','34507_72','34507_76','34509_2','34509_20','34509_30',
                                                       '34510_103','34510_108','34510_117','34510_121','34510_122','34510_127','34510_131','34510_140',
                                                       '34510_381','34510_155','34510_166','34510_167','34510_171','34510_176','34510_177','34510_180',
                                                       '34510_19','34510_194','34510_197','34510_199','34510_207','34510_219','34510_222','34510_228',
                                                       '34510_233','34510_253','34510_263','34510_273','34510_274','34510_276','34510_28','34510_283',
                                                       '34510_287','34510_289','34510_299','34510_30','34510_302','34510_303','34510_313','34510_315',
                                                       '34510_318','34510_32','34510_383','34510_387','34510_388','34510_392','34510_395','34510_400',
                                                       '34510_407','34510_412','34510_416','34510_418','34510_424','34510_428','34510_435','34510_441',
                                                       '34510_445','34510_453','34510_455','34510_457','34510_458','34510_461','34510_495','34510_525',
                                                       '34510_526','34510_527','34510_528','34510_548','34510_55','34510_558','34510_61,','34510_61',
                                                       '34510_88','34510_91','34510_93','34510_323','34510_324','34510_326','34510_329','34510_33',
                                                       '34510_333','34510_337','34510_346','34510_348','34510_352','34510_353','34510_36','34510_361',
                                                       '34510_368','34510_37','34510_370','34510_371','34510_373','34510_38','34510_465','34510_467',
                                                       '34510_478','34510_49','34510_5','34510_501','34510_507','34510_519','34510_522','34510_580',
                                                       '34510_77','34510_80','34510_82','34518_14','34518_19','34518_25','34518_27','34518_28','34518_30',
                                                       '34518_36','34518_37','34518_4','34518_40','34518_42','34518_43','34518_46','34518_8'))


unique(hmm_data_20_2$ID)

summary(hmm_data_20_2$step)
summary(hmm_data_20_2$angle)
hist(hmm_data_20_2$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x,0.75,rate=0.9), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,2,rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_20_2$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.2), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,6), n=200, add=TRUE, col="cadetblue", lwd=4)


#4) fit HMMs!
#look at histograms and get numbers
sl_init_mean <- c(2/0.6, 0.75/0.9)
sl_init_sd <- c(sqrt(2)/0.6, sqrt(0.75)/0.9)
ta_init_mean <- c(0, pi)
ta_init_con <- c(7, 0.2)

mod2 <- fitHMM(data = hmm_data_20_2,
              nbStates = 2,
              stepPar0 = c(sl_init_mean, sl_init_sd),
              anglePar0 = c(ta_init_mean, ta_init_con),
              formula = ~1)

mod2
plot(mod2)

#another round of removing tracks
hmm_data_20_3 <- hmm_data_20_2 %>% dplyr::filter(!ID %in% c('34501_1','34501_13','34501_2','34507_13','34507_26','34507_4',
                                                            '34507_48','34507_53','34507_65','34509_17','34509_19','34509_4',
                                                            '34509_8','34518_33','34510_102','34510_107','34510_181','34510_183','34510_208',
                                                            '34510_211','34510_223','34510_255','34510_259','34510_265','34510_295',
                                                            '34510_31','34510_316','34510_35','34510_389','34510_393','34510_422','34510_425',
                                                            '34510_446','34510_511','34510_538','34510_7'))

summary(hmm_data_20_3$step)
summary(hmm_data_20_3$angle)
hist(hmm_data_20_3$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x,0.75,rate=0.9), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,2,rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_20_3$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.2), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,5), n=200, add=TRUE, col="cadetblue", lwd=4)


#4) fit HMMs!
#look at histograms and get numbers
sl_init_mean <- c(2/0.6, 0.75/0.9)
sl_init_sd <- c(sqrt(2)/0.6, sqrt(0.75)/0.9)
ta_init_mean <- c(0, pi)
ta_init_con <- c(7, 0.2)

mod3 <- fitHMM(data = hmm_data_20_3,
               nbStates = 2,
               stepPar0 = c(sl_init_mean, sl_init_sd),
               anglePar0 = c(ta_init_mean, ta_init_con),
               formula = ~1)

mod3
statecol <- c('lightblue', 'purple1')
plot(mod3, col = statecol)
plot(rmap)





states <- viterbi(mod3)

data.frame(val = rle(states)$values, n = rle(states)$lengths) %>%
  ggplot(aes(val %>% factor, n)) + geom_violin()
mod %>% plotPR()



# viterbi returns vector aligned to rows of the data in best_fit$data
state_20 <- mod3$data %>%
  mutate(state = states)

# Convert angle to degrees for interpretability (angles are in radians)
state_20 <- state_20 %>%
  mutate(angle_deg = angle * 180 / pi,
         abs_angle_deg = abs(angle_deg))

#5a) average metrics per state
state_summary_20 <- state_20 %>%
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
positions_per_track_20 <- state_20 %>%
  group_by(ID, state) %>%
  summarize(n_pos = n(), .groups = "drop")

avg_positions_per_track_20 <- positions_per_track_20 %>%
  group_by(state) %>%
  summarize(mean_positions_per_track = mean(n_pos, na.rm = TRUE),
            median_positions_per_track = median(n_pos, na.rm = TRUE),
            n_track = n())

# Combine summaries
state_metrics_20 <- left_join(state_summary_20, avg_positions_per_track_20, by = "state")

#6) save results, print summaries

# add ID, burst, ts, x, y, step, angle, state to output CSV
output_df <- state_20 %>%
  dplyr::select(ID, track, tag, ts, x, y, step, angle, angle_deg, state)

write_csv(output_df, "Results/interpolated_hmm_states_Rankin_20.csv")

#That was Rankin 20-40 second. Let's try Johnson now w/ 20-40 second

jdat <- data <- read.table('Data/2024-03-18_Rodemann_yaps_johnson.txt', sep = ',', header = T)
jdat <- jdat %>% filter(E < 10)
unique(jdat$tag)
str(jdat)
#anything over 10474 is 20-40 second! Let's subset and redo the tracks

jdat_20 <- jdat %>% dplyr::filter(tag >= 10474)
head(jdat_20)

jdat_20$ts <- as.POSIXct(jdat_20$ts, format = '%Y-%m-%dT%H:%M:%OS')

jdat20 <- jdat_20  %>% dplyr::filter(nobs >= 2)

jdat20t <- jdat20 %>% group_by(tag) %>% mutate(trackstart = ifelse(as.numeric(unclass(ts)-unclass(lag(ts))) > 300, 'Y', 'N')) %>% ungroup()
jdat20t[is.na(jdat20t)] <- 'Y'

str(jdat20t)

jdat20tt <- jdat20t %>% group_by(tag) %>% mutate(track = ifelse(row_number() == 1, 1, NA)) %>% ungroup()

tracknum1 <- function(x){
  for (i in 2:nrow(x)){
    if (x$trackstart[i] == 'Y'){
      x$track[i] <- x$track[i-1]+1
    } else {
      x$track[i] <- x$track[i-1]
    }
  }
  x
}

jdat20ttt <- jdat20tt %>% group_by(tag) %>% group_modify(~tracknum1(.x)) %>% ungroup()

head(jdat20ttt)

j20 <- jdat20ttt %>% mutate(ID = paste(tag, track, sep = "_"))
head(j20)

interp_list_j20 <- j20 %>%
  group_by(.data[[id_col]], .data[[track_col]]) %>%
  group_map(~ interpolate_track(.x, dt_seconds = 40), .keep = TRUE)

# Remove NULLs (tracks that couldn't be interpolated)
interp_list_j20 <- keep(interp_list_j20, ~ !is.null(.x) && nrow(.x) > 0)

interp_df_j20 <- bind_rows(interp_list_j20)
interp_df_j20 <- interp_df_j20 %>% mutate(ID = paste(tag, track, sep = "_"))
head(interp_df_j20)

interp_df_j20 <- interp_df_j20 %>% arrange(ID, ts)
interp_df_j20 <- interp_df_j20 %>% group_by(ID) %>% filter(n()>= 30) %>% ungroup() #track for at least 20 minutes
# prepData: coordinate names are "x" and "y" in our frame
hmm_data_j20 <- prepData(interp_df_j20, type = "UTM", coordNames = c("x", "y"))

# moveHMM uses a column 'step' and 'angle' (radians) created by prepData
# Remove rows with zero or NA step (stay/missing movement) because gamma requires positive steps
hmm_data_j20 <- hmm_data_j20 %>% filter(!is.na(step) & step > 0) %>% filter(!is.na(angle))

hmm_data_j20 <- hmm_data_j20 %>%
  filter(is.finite(step), is.finite(angle), step > 0) %>%
  mutate(
    step = pmax(step, quantile(step, 0.01, na.rm = TRUE)),
    step = pmin(step, quantile(step, 0.99, na.rm = TRUE))
  )

head(hmm_data_j20)

summary(hmm_data_j20$step)
summary(hmm_data_j20$angle)
hist(hmm_data_j20$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x,0.6,rate=1), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,1.75,rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_j20$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.01), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,6), n=200, add=TRUE, col="cadetblue", lwd=4)

#4) fit HMMs!
  #look at histograms and get numbers
  sl_init_mean <- c(1.75/0.6, 0.6/1)
  sl_init_sd <- c(sqrt(1.75)/0.6, sqrt(0.6)/1)
  ta_init_mean <- c(0, pi)
  ta_init_con <- c(6, 0.01)
  
  jmod <- fitHMM(data = hmm_data_j20,
                nbStates = 2,
                stepPar0 = c(sl_init_mean, sl_init_sd),
                anglePar0 = c(ta_init_mean, ta_init_con),
                formula = ~1)
  unique(hmm_data_j20$ID)
  jmod
  plot(jmod)
  
#Take out bad tracks
  
hmm_data_j20_2 <- hmm_data_j20 %>% dplyr::filter(!ID %in% c('10483_100','10483_105','10483_20','10483_21','10483_26','10483_27','10483_41',
                                                            '10483_5','10483_53','10483_62','10483_79','10483_8','10483_80','10483_83','10483_88',
                                                            '10483_90','10486_17','10486_3','10486_8','10487_12','10487_15','10487_4','10487_40','10487_5'))

summary(hmm_data_j20_2$step)
summary(hmm_data_j20_2$angle)
hist(hmm_data_j20_2$step, breaks = 20, freq = F, main = "Step lengths")
curve(dgamma(x,0.6,rate=1), n=200, add=TRUE, col="royalblue", lwd=4)
curve(dgamma(x,1.75,rate=0.6), n=200, add=TRUE, col="cadetblue", lwd=4)
#curve(dgamma(x,1.2,rate=0.02), n=200, add=TRUE, col="navyblue", lwd=4)
hist(hmm_data_j20_2$angle, breaks = 30, freq = F, main = "Turn angles")
curve(circular::dvonmises(x,pi,0.01), n=200, add=TRUE, col="royalblue", lwd=4)
curve(circular::dvonmises(x,0,6), n=200, add=TRUE, col="cadetblue", lwd=4)

sl_init_mean <- c(1.75/0.6, 0.6/1)
sl_init_sd <- c(sqrt(1.75)/0.6, sqrt(0.6)/1)
ta_init_mean <- c(0, pi)
ta_init_con <- c(6, 0.01)

jmod2 <- fitHMM(data = hmm_data_j20_2,
               nbStates = 2,
               stepPar0 = c(sl_init_mean, sl_init_sd),
               anglePar0 = c(ta_init_mean, ta_init_con),
               formula = ~1)

jmod2
plot(jmod2)

jstates <- viterbi(jmod2)

data.frame(val = rle(states)$values, n = rle(states)$lengths) %>%
  ggplot(aes(val %>% factor, n)) + geom_violin()
mod %>% plotPR()



# viterbi returns vector aligned to rows of the data in best_fit$data
jstate_20 <- jmod2$data %>%
  mutate(state = jstates)

# Convert angle to degrees for interpretability (angles are in radians)
jstate_20 <- jstate_20 %>%
  mutate(angle_deg = angle * 180 / pi,
         abs_angle_deg = abs(angle_deg))

#5a) average metrics per state
jstate_summary_20 <- jstate_20 %>%
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
jpositions_per_track_20 <- jstate_20 %>%
  group_by(ID, state) %>%
  summarize(n_pos = n(), .groups = "drop")

javg_positions_per_track_20 <- jpositions_per_track_20 %>%
  group_by(state) %>%
  summarize(mean_positions_per_track = mean(n_pos, na.rm = TRUE),
            median_positions_per_track = median(n_pos, na.rm = TRUE),
            n_track = n())

# Combine summaries
jstate_metrics_20 <- left_join(jstate_summary_20, javg_positions_per_track_20, by = "state")

#6) save results, print summaries

# add ID, burst, ts, x, y, step, angle, state to output CSV
output_jdf <- jstate_20 %>%
  dplyr::select(ID, track, tag, ts, x, y, step, angle, angle_deg, state)

write_csv(output_jdf, "Results/interpolated_hmm_states_Johnson_20.csv")
