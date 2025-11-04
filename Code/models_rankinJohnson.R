###YAPS Rankin II run summary stats####
###Author: Jon Rodemann
#Load in Libraries

library(tidyverse)
library(mlr)
library(tidyverse)
library(randomForest)
library(raster)
library(landscapemetrics)
library(future)
library(furrr)
library(sf)


#load in data
data <- read.table('Data/Rankin_pos.txt', sep = ',', header = T)

head(data)
str(data)

data <- data %>%
  mutate(E = sqrt(x_sd^2 + y_sd^2))

head(data)


#Error less than 10
dat_10 <- data %>% filter(E < 10)
write.csv(dat_10, file = 'Data/Seatrout_positions_sd10_Rankin.csv')


#summary stats
#number of positions per fish
datnfish <- dat_10 %>% group_by(tag) %>% summarize(n = n())
write.csv(datnfish, file = 'Seatrout_positions_per_fish_Rankin.csv')

#abacus plot for positions
dat_10$ts <- as.POSIXct(dat_10$ts, format = '%Y-%m-%dT%H:%M:%OS')


library(splitstackshape) ###package to use cSplit
library(chron)
library(ggplot2)
library(scales)

# Break up the new datetime column into  a date, and a time
#Jon's Data
data <- dat_10

data$Date_Time2 <- data$ts

data <- cSplit(data, "Date_Time2", sep = " ", type.convert = F)

head(data)


# Rename columns so they make sense, it'll just be your last 2 column numbers, in this case the 15th and 16th column s
colnames(data)[11:13]<-c("E", "Date", "Time")
head(data)


# Then I repeat this and parse the date into year, month, and day (then hour, minute, second), so I can easily subset by year
# Copy date to parse out, maintain original date
data$Date2 <- data$Date


data<-cSplit(data, "Date2", sep = "-", type.convert = FALSE)
head(data)
colnames(data)<-c("tag", "ts", "x", "y", "z", "x_sd", "y_sd", "z_sd", "nobs", "k", "E", "Date", "Time", "Year", "Month", "Day")
head(data)

data$Time2 <- data$Time
head(data)

data<-cSplit(data, "Time2", sep = ":", type.convert = FALSE)
head(data)
colnames(data)[17:19]<-c("Hour", "Minute", "Second")
head(data)

# Assign classes so data plays nice for plotting and summaries
data$Date <- as.Date(data$Date)
data$Time <- as.times(data$Time)
data$Year <- as.numeric(data$Year)
data$Month <- as.numeric(data$Month)
data$Day <- as.numeric(data$Day)
data$Hour <- as.numeric(data$Hour)
data$Minute <- as.numeric(data$Minute)
data$Second <- as.numeric(data$Second)
data$tag <- as.factor(data$tag)


data$ts <- as.POSIXct(data$ts, format='%Y-%m-%d %H:%M:%S', tz='America/New_York')

head(data)

trout_det <- ggplot(data,aes(ts, tag)) + 
  #geom_vline(aes(xintercept = (as.Date("2020-03-01"))), linetype=4, size=0.75, color="red")+    
  #geom_vline(aes(xintercept = (as.Date("2022-02-23"))), linetype=4, size=0.75, color="red")+    
  geom_point() +
  #scale_color_gradient(low="darkblue", high="lightgreen") +
  #scale_x_date(date_breaks = "4 months", limits = as.Date(c("2022-01-15", "2022-03-01")))+
  scale_x_datetime(labels = date_format("%Y-%m-%d"),
                   date_breaks = "2 weeks", limits = as.POSIXct(c('2021-11-13 00:00:00', '2022-05-13 00:00:00')))+
  labs(x = "Date",  y = "Transmitter") +
  # Change the angle and text size of x-axis tick labels to make more readable in final plots
  theme_bw()+
  theme(axis.text.x=element_text(angle= 50, size=10, vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y=element_text(size=8))+
  labs(title = "Trout Detections YAPS")
trout_det

ggsave(trout_det, file = 'YAPS_abacus_Rankin.png')

#tracks for each tag

#create columns for track investigation
dat_20 <- data
dat_20$track <- NA

str(dat_20)

dat_10$ts1 <- as.POSIXct(dat_10$ts, format='%Y-%m-%dT%H:%M:%OS', tz='America/New_York', origin = '1970-01-01 00:00:00')

dat_10t <- dat_10 %>% group_by(tag) %>% mutate(trackstart = ifelse(as.numeric(unclass(ts1)-unclass(lag(ts1))) > 3600, 'Y', 'N')) %>% ungroup()
dat_10t[is.na(dat_10t)] <- 'Y'

str(dat_10t)

dat_10tt <- dat_10t %>% group_by(tag) %>% mutate(track = ifelse(row_number() == 1, 1, NA)) %>% ungroup()

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

dat_10ttt <- dat_10tt %>% group_by(tag) %>% group_modify(~tracknum1(.x)) %>% ungroup()

#Add step length and turn angle
dat_10ang <- dat_10ttt %>% group_by(tag, track) %>% mutate(step_length = sqrt((x - lag(x, default = NA))^2 + (y - lag(y, default = NA))^2), bear = atan2(
  y - lag(y, default = NA),
  x - lag(x, default = NA)
) * (180 / pi), bear_utm = ifelse(bear< 0, bear + 360, bear), angle = lag(bear_utm, default = NA)-bear_utm) %>% ungroup()

write.csv(dat_10ang, file = 'Data/Rankin_track_dat.csv')
dat10 <- dat_10ttt %>% group_by(tag, track) %>% summarize(n = n(), time = max(as.numeric(ts)) - min(as.numeric(ts))) 

plot(map_dense)
plot(mapag1)
head(dat_10ang)

d10 <- dat_10ang %>% st_as_sf(coords = c("x", "y"), crs = st_crs(map_dense))

d10_1 <- d10 %>% mutate(class = extract(mapag1, d10))

d10_1_angle <- as.data.frame(d10_1 %>% mutate(angle_pos = abs(angle)) %>% dplyr::select(c(angle_pos, tag, class, step_length)) %>% mutate(angle = if_else(angle_pos == 0, 0.000001, angle_pos)))

d11 <- d10_1_angle %>% dplyr::select(-geometry)
max(d10$ts)

met_all <- read.csv('E:/FIU/Project/YAPS/Rankin_YAPS/met_all.csv')

data <- cbind(d11, met_all)

library(glmmTMB)
library(performance)

setwd('C:/Users/jonro/OneDrive/Desktop/Tom_SAV/data/')
source("HighstatLibV13.R")
options(na.action = "na.fail")

head(data)

dv = c('div','pland', 'area_mn', 'pd', 'ed', 'para_mn', 'gyrate_mn', 'shape_mn', 'ai')
Mypairs(data[, dv])
data1 <- data %>% drop_na()

head(data1)

data1 <- data1 %>% mutate(angle = if_else(angle_pos == 0, 0.000001, angle_pos))

model1 <- glmmTMB(step_length ~ as.factor(class) +pland + pd + ed + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model1)
check_model(model1)
performance(model1)

library(MuMIn)

sl = dredge(model1, rank = 'AIC') %>% 
  filter(delta < 2)
sl

model2 <- glmmTMB(step_length ~ as.factor(class) +pland + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
performance(model2)
summary(model2)


plot(model2)
library(ggeffects)
max(data1$pland)
min(data1$pland)
values_seq <- seq(0, 100, by = 5)

# Convert the sequence to a character string with values separated by commas
values_str <- paste(values_seq, collapse = ",")

# Use the sequence in ggpredict
pland_best <- ggpredict(model2, terms = paste("pland [", values_str, "]", sep=""))

ggplot(pland_best, aes(x, predicted))+
  geom_line(color = '#9E1B32', size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2)+
  labs(x = "PLAND",
       y = 'Step Length (m)')+
  theme_classic()

ggsave(filename = 'pland_SL.png')

max(data1$shape_mn)
min(data1$shape_mn)
values_seq <- seq(1, 3.5, by = 0.1)

# Convert the sequence to a character string with values separated by commas
values_str <- paste(values_seq, collapse = ",")

# Use the sequence in ggpredict
pland_best <- ggpredict(model2, terms = paste("shape_mn [", values_str, "]", sep=""))

ggplot(pland_best, aes(x, predicted))+
  geom_line(color = '#9E1B32', size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Step Length (m)')+
  theme_classic()

ggsave(filename = 'Shape_MN_SL.png')

class_best <- ggpredict(model2, terms = "class")
str(class_best)
class_best$tp <- as.numeric(class_best$group)

col = c('tan', 'lightgreen', 'darkgreen')

ggplot(class_best, aes(x, predicted, colour = x))+
  scale_colour_manual(values = col)+
  geom_point(size = 6) +
  geom_segment(aes(y = conf.low, yend = conf.high, x = x), size = 2)+
  labs(x = "Class",
       y = 'Step Length (m)')+
  theme_classic()+
  theme(legend.position="none")

ggsave(filename = 'E:/FIU/PostDoc/CESI/Final_report/YAPS/class_steplength.png')

data1$angle_pos <- as.numeric(data1$angle_pos)
model3 <- glmmTMB(angle ~ as.factor(class) +pland + pd + ed + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model3)
check_model(model1)
performance(model3)

an = dredge(model2, rank = 'AIC') %>% 
  filter(delta < 2)
an

model3 <- glmmTMB(angle ~ pland + ed + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model3)
performance(model3)



abs(d10_1$angle)

###GAMS. RF won't work for this data. Try RSFs later
library(mgcv)

head(data1)

model <- gam(angle_pos ~ s(pland, k = 4) + s(pd, k = 4)  + s(ed, k = 4) + class + s(shape_mn, k = 4) +
               s(tag, bs = "re"),
             family = nb, data = data1, method = "REML")
summary(model)
check_model(model)
r.squaredGLMM(model)
performance(model)
r2(model)

plot(model, pages = 1, all.terms = TRUE)

model1 <- gam(step_length ~ s(pland, k = 4) + s(pd, k = 4)  + s(ed, k = 4) + class + s(shape_mn, k = 4) +
               s(tag, bs = "re"),
             family = tw, data = data1, method = "REML")
summary(model1)
check_model(model)
r.squaredGLMM(model)
performance(model1)
r2(model)


compare_performance(model1, model2)
plot(model1, pages = 1, all.terms = TRUE)


##IDK about GAMs... Let's do RSF framework. This will require multiple cores for the landscape metrics
library(terra, exclude = 'resample') #work with spatial data - rasters
library(raster) #work with spatial data - rasters
library(sf) #Work with spatial data - shapefiles
library(sfheaders) #work with spatial data
library(chron) #visualization of acoustic data
library(splitstackshape) #break up data into smaller columns
library(scales)#visualization of acoustic data
library(mlr3verse) # Machine learning
library(mlr3spatial) #spatial machine learning
library(randomForest) #Machine learning
library(iml) #result interpretation
library(ranger) #Machine learning
library(ggmap) #plotting onto map
library(beepr) #beeps when code is done running
library(tidyverse)

#make random points

extent <- resf

randLocs <- sf::st_sample(extent, size = 45540, type = 'random') %>% st_transform(2958) %>% sfc_to_df()

#have random points. Need to extract classes, then run landscape metrics
rand <- randLocs %>% st_as_sf(coords = c("x", "y"), crs = st_crs(map_dense))

rand1 <- rand %>% mutate(class = extract(mapag1, rand))

rand1 <- rand1 %>% st_as_sf(coords = c("x", "y"), crs = "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs")

coord <- as.data.frame(cbind(randLocs$x, randLocs$y))
coord_sf <- st_as_sf(coord, coords = c("V1", "V2"), crs = "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs")
coord1 <- coord_sf[st_within(coord_sf, resf, sparse = FALSE), ]
rand2 <- rand1[st_within(rand1, resf, sparse = FALSE),]

library(furrr)
future::plan(multisession, workers = 9)
result_metrand <- 1:nrow(coord1) %>% future_map(
  function(i) {
    sample_lsm(
      map_dense,
      y = st_geometry(coord1[i, ]),  # Ensure geometry is passed
      size = 5,
      what = c(
        "lsm_c_pd", "lsm_c_ed",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)
future::plan(sequential)

rt <-lapply(result_metrand, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})
met <- bind_rows(rt, .id = 'source')
unique_sources <- unique(met$source)
unique_classes <- unique(met$class)
unique_metric <- unique(met$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(met, by = c("source", "class", 'metric')) %>% distinct()
met_1 <- mettry %>% dplyr::filter(class == 1) %>% mutate(source = as.numeric(source)) %>% arrange(source)
rep4<-function(d){
  d$newcol =  rep(1:ceiling(nrow(d)/4), each = 4)[1:nrow(d)]
  d
}

met_1 <- rep4(met_1)
#met_2 <- rep1(met_3)
met_1 <- met_1 %>% dplyr::select(-plot_id) %>% dplyr::rename(plot_id = newcol)

met_pland <- met_1 %>% dplyr::filter(metric == 'pland') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pland' = value)
met_pd <- met_1 %>% dplyr::filter(metric == 'pd') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pd' = value)
met_ed <- met_1 %>% dplyr::filter(metric == 'ed') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('ed' = value)
met_shape_mn <- met_1 %>% dplyr::filter(metric == 'shape_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('shape_mn' = value)

met_all <- list(met_pland, met_pd, met_ed, met_shape_mn)
#met_all1 <- list(met_pd, met_ed, met_np, met_te, met_ca, met_area_mn, met_shape_mn, met_div)
met_all <- met_all %>% purrr::reduce(full_join, by='plot_id')

rand_dat <- cbind(rand2, met_all)
rand_dat <- rand_dat %>% dplyr::select(-c(sfg_id, point_id))

head(data1)

data1$pa <- 1
rand_dat$pa <- 0

head(rand_dat)

data2 <- data1 %>% dplyr::select(c(class, pland, plot_id, pd, ed, shape_mn, pa))
rand_dat2 <- rand_dat %>% as.data.frame() %>% dplyr::select(-geometry)

head(data2)
head(rand_dat2)

rsf_dat <- rbind(data2, rand_dat2)

head(rsf_dat)

set.seed(19)

write.csv(file = 'rsf_rankin.csv', rsf_dat)

rsf_dat$pa <- as.factor(rsf_dat$pa)

# Randomly select 70% of the data frame for the training dataset
RSF_ar.train <- rsf_dat[sample(1:nrow(rsf_dat), nrow(rsf_dat) * 0.7, replace = FALSE), ]
# Remainder (i.e., 30%) becomes the test dataset.
RSF_ar.test <- rsf_dat[!(rsf_dat$plot_id %in% RSF_ar.train$plot_id), ] 

head(RSF_ar.test)

#take out ID column
RSF_ar.test <- RSF_ar.test %>% dplyr::select(-plot_id) %>% drop_na()
RSF_ar.train <- RSF_ar.train %>% dplyr::select(-plot_id) %>% drop_na()


# Set tasks for training and test datasets.
task_trout.train <- as_task_classif(
  RSF_ar.train, target = "pa", positive = '1',
)
str(task_trout.train)
task_trout.train

task_trout.test <- as_task_classif(
  x = RSF_ar.test, target = "pa",
  positive = "1"
)

# Make learner.
learner <-lrn(
  "classif.ranger",
  predict_type = "prob",
  mtry  = to_tune(1, ncol(RSF_ar.train) - 4),
  sample.fraction = to_tune(0.2, 0.9),
  min.node.size = to_tune(1,10),
  importance = 'impurity'
)

#tune hyperparameters
instance = ti(
  task = task_trout.train,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measures = msr("classif.ce"),
  terminator = trm("none")
)

tuner = mlr3tuning::tnr("grid_search", resolution = 5, batch_size = 5)
tuner$optimize(instance) # Takes ~ 4 minutes on my relatively fast computer
beep(1)

#store tuned hyperparameters in learner
learner$param_set$values <- instance$result_learner_param_vals

#finally! We can train our model with the train function
learner$train(task_trout.train)

#let's quickly look at the model
learner$model

#Accuracy of model - first on training data, then testing data
measures <- msrs(c('classif.acc'))

pred_train <- learner$predict(task_trout.train)
pred_train$confusion
pred_train$score(measures)


pred_test <- learner$predict(task_trout.test)
pred_test$confusion
pred_test$score(measures)


#importance with iml package - this is looking at the most influencial predictors in the model
x_trout <- RSF_ar.train %>% dplyr::select(-pa) 
# Create "Predictor" object to interpret findings via the iml package.
predictor_trout <- Predictor$new(learner, data = x_trout, y = RSF_ar.train$pa) 
options(future.globals.maxSize = 1024 * 1024 * 1024) 
imp_trout <- FeatureImp$new(predictor_trout, loss = "ce") # Calculate importance.

imp_trout$plot()+
  geom_bar(stat = 'identity')+
  scale_color_continuous()+
  theme_classic()

imp$feature <- factor(imp$feature, levels = c('class', 'pd', 'shape_mn', 'ed', 'pland'))
x = c('darkgreen', 'forestgreen', 'green3', 'lightgreen', 'lightblue')
ggplot(imp, aes(x = feature, y = importance))+
  geom_bar(stat = 'identity', fill = x)+
  coord_flip()+
  theme_classic()


#pland most important
effect_pland <- FeatureEffect$new(predictor_trout, feature = c('pland'), method = 'pdp')
beep(1)
effect_pland$plot()
effect_pland
p_res <- effect_pland$results
p_res <- p_res %>% filter(.class == 1)

ggplot(p_res, aes(x = pland, y = .value))+
  geom_smooth(color = 'darkgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "PLAND",
       y = 'Probability of detection')+
  theme_classic()

#
effect_ed <- FeatureEffect$new(predictor_trout, feature = c('ed'), method = 'pdp')
beep(1)
effect_ed$plot()

ed_res <- effect_ed$results
ed_res <- ed_res %>% filter(.class == 1)

ggplot(ed_res, aes(x = ed, y = .value))+
  geom_smooth(color = 'forestgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Edge Density",
       y = 'Probability of detection')+
  theme_classic()

effect_pd <- FeatureEffect$new(predictor_trout, feature = c('pd'), method = 'pdp')
effect_pd$plot()

pd_res <- effect_pd$results
pd_res <- pd_res %>% filter(.class == 1)

ggplot(pd_res, aes(x = pd, y = .value))+
  geom_smooth(color = 'lightgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Patch Density",
       y = 'Probability of detection')+
  theme_classic()

effect_shape <- FeatureEffect$new(predictor_trout, feature = c('shape_mn'), method = 'pdp')
effect_shape$plot()


sh_res <- effect_shape$results
sh_res <- sh_res %>% filter(.class == 1)

ggplot(sh_res, aes(x = shape_mn, y = .value))+
  geom_smooth(color = 'green3', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Probability of detection')+
  theme_classic()
#random vs our tracks

##### Visualization########
dat20 <- merge(dat20, timing, by = 'tag')
dat20 <- dat20 %>% dplyr::select(-X)
dat20$hours <- dat20$time/3600


d <- dat20 %>% group_by(tag) %>% filter(time > 3600) %>% summarize(n = n())
d2 <- dat20 %>% group_by(tag) %>% filter(time > 7200) %>% summarize(n = n())
d3 <- dat20 %>% group_by(tag) %>% summarize(max = max(time)/3600)
dshort <- dat20 %>% filter(timing == '20-40')
dlong <- dat20 %>% filter(timing == '60-90')

ggplot(data = dat20, aes(hours, colour = timing))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_all_timing.png')

ggplot(data = dat20, aes(hours, colour = timing))+
  geom_density()+
  theme_classic()+
  xlim(0,20)
ggsave(filename = 'density_plot_20_timing.png')

ggplot(data = dat20, aes(hours, colour = timing))+
  geom_density()+
  theme_classic()+
  xlim(0,5)
ggsave(filename = 'density_plot_5_timing.png')

dat8 <- dat20 %>% filter(tag %in% c('875', '880', '881', '912', '34501', '34505', '34507', '34509', '34510', '34518', '47421'))
head(dat8)
dat8$tag <- as.factor(dat8$tag)

ggplot(data = dat8, aes(hours, colour = tag))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_all_tags.png')

ggplot(data = dat8, aes(hours, colour = tag))+
  geom_density()+
  theme_classic()+
  xlim(0,20)
ggsave(filename = 'density_plot_20_tags.png')

d <- dat8 %>% filter(tag == '875')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_875.png')

d <- dat8 %>% filter(tag == '880')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_880.png')

d <- dat8 %>% filter(tag == '881')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_881.png')

d <- dat8 %>% filter(tag == '912')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_912.png')

d <- dat8 %>% filter(tag == '34501')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34501.png')

d <- dat8 %>% filter(tag == '34505')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34505.png')

d <- dat8 %>% filter(tag == '34507')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34507.png')

d <- dat8 %>% filter(tag == '34509')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34509.png')

d <- dat8 %>% filter(tag == '34510')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34510.png')

d <- dat8 %>% filter(tag == '34518')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_34518.png')

d <- dat8 %>% filter(tag == '47421')
ggplot(data = d, aes(hours))+
  geom_density()+
  theme_classic()
ggsave(filename = 'density_plot_47421.png')

datnt <- dat20 %>% group_by(tag) %>% summarize(n =n())

write.csv(datnt, file = 'num_tracks_fish.csv')

datav <- dat20 %>% group_by(tag) %>% summarize(time = mean(time), pos = mean(n))

write.csv(datav, file = 'av_tracks_fish.csv')

###Johnson models####
setwd('E:/FIU/PostDoc/CESI/Final_report/YAPS/')
data <- read.csv('E:/FIU/PostDoc/CESI/Final_report/YAPS/jdat_final.csv')
library(glmmTMB)
library(performance)
library(tidyverse)

setwd('C:/Users/jonro/OneDrive/Desktop/Tom_SAV/data/')
source("HighstatLibV13.R")
options(na.action = "na.fail")

head(data)

data1 <- data %>% drop_na()

head(data1)

data1 <- data1 %>% mutate(angle = abs(if_else(angle == 0, 0.000001, angle)))
data1 <- data1 %>% mutate(step_length = abs(if_else(step_length == 0, 0.000001, step_length)))
names(data1)

model1 <- glmmTMB(step_length ~ as.factor(johnclass2) +pland + pd + ed + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model1)
check_model(model1, simulate = FALSE, check = c("normality", "linearity", "homogeneity"), residual_type = 'normal')
performance(model1)

library(MuMIn)

sl = dredge(model1, rank = 'AIC') %>% 
  filter(delta < 2)
sl

model2 <- glmmTMB(step_length ~ as.factor(johnclass2) + pd + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
check_model(model2)
performance(model2)
summary(model2)


plot(model2)
library(ggeffects)
max(data1$pd)
min(data1$pd)
values_seq <- seq(1205, 13290, by = 5)

# Convert the sequence to a character string with values separated by commas
values_str <- paste(values_seq, collapse = ",")

# Use the sequence in ggpredict
pd_best <- ggpredict(model2, terms = paste("pd [", values_str, "]", sep=""))

ggplot(pd_best, aes(x, predicted))+
  geom_line(color = '#9E1B32', size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2)+
  labs(x = "Patch Density",
       y = 'Step Length (m)')+
  theme_classic()

ggsave(filename = 'Figs_john/pd_SL.png')

max(data1$shape_mn)
min(data1$shape_mn)
values_seq <- seq(1, 3.7, by = 0.1)

# Convert the sequence to a character string with values separated by commas
values_str <- paste(values_seq, collapse = ",")

# Use the sequence in ggpredict
shape_best <- ggpredict(model2, terms = paste("shape_mn [", values_str, "]", sep=""))

ggplot(shape_best, aes(x, predicted))+
  geom_line(color = '#9E1B32', size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Step Length (m)')+
  theme_classic()

ggsave(filename = 'Figs_john/Shape_MN_SL.png')

class_best <- ggpredict(model2, terms = "johnclass2")
str(class_best)
class_best$tp <- as.numeric(class_best$group)

col = c('lightgreen', 'darkgreen')

ggplot(class_best, aes(x, predicted, colour = x))+
  scale_colour_manual(values = col)+
  geom_point(size = 6) +
  geom_segment(aes(y = conf.low, yend = conf.high, x = x), size = 2)+
  labs(x = "Class",
       y = 'Step Length (m)')+
  scale_x_discrete(labels = c("0" = "Sparse (<45%)", "1" = "Dense (>45%)")) +
  theme_classic()+
  theme(legend.position="none")

ggsave(filename = 'Figs_john/class_steplength.png')

data1$angle <- as.numeric(data1$angle)
model3 <- glmmTMB(angle ~ as.factor(johnclass2) +pland + pd + ed + shape_mn + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model3)
check_model(model3)
performance(model3)

an = dredge(model3, rank = 'AIC') %>% 
  filter(delta < 2)
an

model3 <- glmmTMB(angle ~ shape_mn + ed + (1|tag), data = data1, family = Gamma(link = 'log'))
summary(model3)
performance(model3)




##Let's do RSF framework. This will require multiple cores for the landscape metrics
library(terra, exclude = 'resample') #work with spatial data - rasters
library(raster) #work with spatial data - rasters
library(sf) #Work with spatial data - shapefiles
library(sfheaders) #work with spatial data
library(chron) #visualization of acoustic data
library(splitstackshape) #break up data into smaller columns
library(scales)#visualization of acoustic data
library(mlr3verse) # Machine learning
library(mlr3spatial) #spatial machine learning
library(randomForest) #Machine learning
library(iml) #result interpretation
library(ranger) #Machine learning
library(ggmap) #plotting onto map
library(beepr) #beeps when code is done running
library(tidyverse)
library(landscapemetrics)
#make random points

jmap1 <- rast('jmap_dense.tif')
plot(jmap1)
raster_bbox <- as.polygons(jmap1) |> st_as_sf()
crs(raster_bbox)
raster_bbox$merged_class <- 1  # assign same class
merged <- st_union(raster_bbox)
extent <- merged
plot(merged)
randLocs <- sf::st_sample(extent, size = 24015, type = 'random') %>% st_transform(2958) %>% sfc_to_df()

#have random points. Need to extract classes, then run landscape metrics
rand <- randLocs %>% dplyr::select(c(x,y)) %>% st_as_sf(coords = c("x", "y"), crs = st_crs(jmap1))

randclass <- rand %>% mutate(class = extract(jmap1, rand))

randclass1 <- as.data.frame(randclass$class)
str(randclass1)

library(furrr)
future::plan(multisession)
result_metrand <- 1:nrow(rand) %>% future_map(
  function(i) {
    jmap_r <- rast("jmap_dense.tif")
    sample_lsm(
      jmap_r,
      y = st_geometry(rand[i, ]),  # Ensure geometry is passed
      size = 15,
      what = c(
        "lsm_c_pd", "lsm_c_ed",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)
future::plan(sequential)

rt <-lapply(result_metrand, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})
met <- bind_rows(rt, .id = 'source')
unique_sources <- unique(met$source)
unique_classes <- unique(met$class)
unique_metric <- unique(met$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(met, by = c("source", "class", 'metric')) %>% distinct()
mettry1 <- mettry %>% dplyr::filter(class == 1) %>% pivot_wider(names_from = metric, values_from = value)
rand_dat <- mettry1
head(data1)

data1$pa <- 1
rand_dat$pa <- 0
rand_dat <- cbind(rand_dat, randclass1)
head(rand_dat)
str(rand_dat)
head(data1)
data2 <- data1 %>% rename(class = johnclass2) %>% dplyr::select(c(class, pland, pd, ed, shape_mn, pa))
names(rand_dat)
rand_dat2 <- rand_dat %>%  dplyr::select(c(johnclass2, pland, pd, ed, shape_mn, pa))

head(data2)
names(rand_dat2)
rand_dat2 <- rand_dat2 %>% rename(class = johnclass2)
str(rand_dat2)
head(rand_dat2)

rsf_dat <- rbind(data2, rand_dat2)


set.seed(19)

write.csv(file = 'rsf_John.csv', rsf_dat)

set.seed(19)

rsf_dat$pa <- as.factor(rsf_dat$pa)

rsf_dat$plot_id <- 1:nrow(rsf_dat)

# Randomly select 70% of the data frame for the training dataset
RSF_ar.train <- rsf_dat[sample(1:nrow(rsf_dat), nrow(rsf_dat) * 0.7, replace = FALSE), ]
# Remainder (i.e., 30%) becomes the test dataset.
RSF_ar.test <- rsf_dat[!(rsf_dat$plot_id %in% RSF_ar.train$plot_id), ] 

head(RSF_ar.test)

#take out ID column
RSF_ar.test <- RSF_ar.test %>% dplyr::select(-plot_id) %>% drop_na()
RSF_ar.train <- RSF_ar.train %>% dplyr::select(-plot_id) %>% drop_na()


# Set tasks for training and test datasets.
task_trout.train <- as_task_classif(
  RSF_ar.train, target = "pa", positive = '1',
)
str(task_trout.train)
task_trout.train

task_trout.test <- as_task_classif(
  x = RSF_ar.test, target = "pa",
  positive = "1"
)

# Make learner.
learner <-lrn(
  "classif.ranger",
  predict_type = "prob",
  mtry  = to_tune(1, ncol(RSF_ar.train) - 4),
  sample.fraction = to_tune(0.2, 0.9),
  min.node.size = to_tune(1,10),
  importance = 'impurity'
)

#tune hyperparameters
instance = ti(
  task = task_trout.train,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measures = msr("classif.ce"),
  terminator = trm("none")
)

tuner = mlr3tuning::tnr("grid_search", resolution = 5, batch_size = 5)
tuner$optimize(instance) # Takes ~ 4 minutes on my relatively fast computer
beep(1)

#store tuned hyperparameters in learner
learner$param_set$values <- instance$result_learner_param_vals

#finally! We can train our model with the train function
learner$train(task_trout.train)

#let's quickly look at the model
learner$model

#Accuracy of model - first on training data, then testing data
measures <- msrs(c('classif.acc'))

pred_train <- learner$predict(task_trout.train)
pred_train$confusion
pred_train$score(measures)


pred_test <- learner$predict(task_trout.test)
pred_test$confusion
pred_test$score(measures)


#importance with iml package - this is looking at the most influencial predictors in the model
x_trout <- RSF_ar.train %>% dplyr::select(-pa) 
# Create "Predictor" object to interpret findings via the iml package.
predictor_trout <- Predictor$new(learner, data = x_trout, y = RSF_ar.train$pa) 
options(future.globals.maxSize = 1024 * 1024 * 1024) 
imp_trout <- FeatureImp$new(predictor_trout, loss = "ce") # Calculate importance.
x = c('darkgreen', 'forestgreen', 'green3', 'lightgreen', 'lightblue')
p <- imp_trout$plot()

# Remove all point layers
p$layers <- Filter(function(layer) {
  !inherits(layer$geom, c("GeomPoint", "GeomErrorbar", "GeomErrorbarh", 'GeomSegment'))
}, p$layers)

# Add your clean bars
p +
  geom_bar(stat = "identity", fill = x) +
  theme_classic()

ggsave(filename = 'Figs_john/rsf_importance_john.png')

#pd only importance
effect_pd <- FeatureEffect$new(predictor_trout, feature = c('pd'), method = 'pdp')
beep(1)
effect_pd$plot()
effect_pd
p_res <- effect_pd$results
p_res <- p_res %>% filter(.class == 1)

ggplot(p_res, aes(x = pd, y = .value))+
  geom_smooth(color = 'darkgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Patch Density",
       y = 'Probability of detection')+
  theme_classic()

ggsave(filename = 'Figs_john/rsf_pd_john.png')
#
effect_ed <- FeatureEffect$new(predictor_trout, feature = c('ed'), method = 'pdp')
beep(1)
effect_ed$plot()

ed_res <- effect_ed$results
ed_res <- ed_res %>% filter(.class == 1)

ggplot(ed_res, aes(x = ed, y = .value))+
  geom_smooth(color = 'forestgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Edge Density",
       y = 'Probability of detection')+
  theme_classic()

effect_pd <- FeatureEffect$new(predictor_trout, feature = c('pd'), method = 'pdp')
effect_pd$plot()

pd_res <- effect_pd$results
pd_res <- pd_res %>% filter(.class == 1)

ggplot(pd_res, aes(x = pd, y = .value))+
  geom_smooth(color = 'lightgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Patch Density",
       y = 'Probability of detection')+
  theme_classic()

effect_shape <- FeatureEffect$new(predictor_trout, feature = c('shape_mn'), method = 'pdp')
effect_shape$plot()


sh_res <- effect_shape$results
sh_res <- sh_res %>% filter(.class == 1)

ggplot(sh_res, aes(x = shape_mn, y = .value))+
  geom_smooth(color = 'green3', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Probability of detection')+
  theme_classic()
