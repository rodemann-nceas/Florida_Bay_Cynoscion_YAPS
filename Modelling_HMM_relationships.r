####Relating HMMs with environmental data
library(tidyverse)
library(randomForest)
library(raster)
library(landscapemetrics)
library(future)
library(furrr)
library(sf)
library(glmmTMB)
library(performance)
library(MuMIn)
library(ggeffects)
library(terra)
library(mlr3verse) # Machine learning
library(mlr3spatial) #spatial machine learning
library(iml) #result interpretation
library(ranger) #Machine learning
library(ggmap) #plotting onto map
library(beepr) #beeps when code is done running
library(mgcv)

#load in data
jmap <- rast('Data/Maps/jmap_dense.tif')
plot(jmap)
rmap <- rast('Data/Maps/savclass_2.tif')
plot(rmap)
jdat <- read_csv('Results/interpolated_hmm_states_Johnson_20.csv')
rdat <- read_csv('Results/interpolated_hmm_states_Rankin_20.csv')

#Aggregate and resample Rankin map
#aggregate to 1 m
str(rmap)
rmap <- terra::aggregate(rmap, fact = 3, fun='modal')
plot(rmap)

reclass_df <- c(0, 1, 0,
                1, 2, 0,
                2, 3, 1)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
rmap <- classify(rmap,
                        reclass_m)
plot(rmap)
writeRaster(jmap, filename = 'Data/Johnson_map.tif')
writeRaster(rmap, filename = 'Data/Rankin_map.tif')
#extract values to points from HMMs
#make hmm datasets spatial objects
jdats <- st_as_sf(jdat, coords = c("x", "y"), crs = st_crs(jmap))
plot(jmap)
plot(jdats, add = T, pch = 19, col = jdats$state)

rdats <- st_as_sf(rdat, coords = c("x", "y"), crs = st_crs(rmap))
plot(rmap)
plot(rdats, add = T, pch = 19, col = rdats$state)

#extract
rdate <- terra::extract(rmap, rdats, xy = T, ID = F) %>% bind_cols(rdats) %>% mutate(class = if_else(savclass_2 == 1,2,1)) %>% drop_na() %>% dplyr::select(-savclass_2)
jdate <- terra::extract(jmap, jdats, xy = T, ID = F) %>% bind_cols(jdats) %>% mutate(class = if_else(johnclass2 == 1,2,1)) %>% drop_na() %>% dplyr::select(-johnclass2)

#sucess! Let's look at first relationships between state and class
rdate$class <- as.factor(rdate$class)
rdate$state <- as.factor(rdate$state)
rdate %>% group_by(state) %>% mutate(n = n()) %>% ungroup() %>% group_by(class, state) %>% summarize(n = n()/n) %>% distinct() 
rchi <- chisq.test(rdate$class, rdate$state)
rchi$observed
rchi$expected
rchi$statistic
rchi$p.value


jdate$class <- as.factor(jdate$class)
jdate$state <- as.factor(jdate$state)
jdate %>% group_by(state) %>% mutate(n = n()) %>% ungroup() %>% group_by(class, state) %>% summarize(n = n()/n) %>% distinct() 
jchi <- chisq.test(jdate$class, jdate$state)
jchi$observed
jchi$expected
jchi$statistic
jchi$p.value

#Looks like more effects for Johnson. Ok, let's run some landscape metrics
#Rankin first
rcoord <- as.data.frame(cbind(rdate$x, rdate$y))
rcoord_sf <- st_as_sf(rcoord, coords = c("V1", "V2"), crs = st_crs(rmap))
?st_crs
library(furrr)
plan(multisession)  # Use multisession to run on multiple cores
library(purrr)

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

rankd <- metr1 %>% dplyr::select(c(division, ed, pd, pland, shape_mn)) %>% bind_cols(rdate)
rankdat <- rankd %>% as.data.frame %>% dplyr::select(-geometry)
write.csv(rankdat, file = 'Results/Rankin_full_hmm.csv')
#rankd is the Rankin dataset to run models on

#now Johnson
jcoord <- as.data.frame(cbind(jdate$x, jdate$y))
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

jt <-lapply(result_metj, as.data.frame)

jt<- lapply(jt, function(x) if (is.null(x)) data.frame() else x)

all_columnsj <- unique(unlist(lapply(jt, colnames)))
jt <- lapply(jt, function(df) {
  df[setdiff(all_columnsj, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columnsj]  # Reorder columns
  df
})
metj <- bind_rows(jt, .id = 'source')

unique_sources <- unique(metj$source)
unique_classes <- unique(metj$class)
unique_metric <- unique(metj$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettryj <- full_combinations %>%
  left_join(metj, by = c("source", "class", 'metric')) %>% distinct()

metj1 <- mettryj %>% pivot_wider(names_from = metric, values_from = value) %>% dplyr::filter(class == 1)
names(metj1)

johnd <- metj1 %>% dplyr::select(c(division, ed, pd, pland, shape_mn)) %>% bind_cols(jdate)
johndat <- johnd %>% as.data.frame() %>% dplyr::select(-geometry)
write.csv(johndat, file = 'Results/Johnson_full_hmm.csv')

#Let's try some modelling! - none are significant before RF

#Rankin
rankd2 <- rankd %>% mutate(state_int = if_else(state == '2',1,0)) %>% drop_na()
rankmod <- glmmTMB(state_int ~ class +pland + pd + ed + shape_mn + (1|tag), data = rankd2, family = binomial)
summary(rankmod)
performance(rankmod)
check_model(rankmod)

options(na.action = "na.fail")
sl = dredge(rankmod, rank = 'AIC') %>% 
  filter(delta < 2)
sl

johnd2 <- johnd %>% mutate(state_int = if_else(state == '2',1,0)) %>% drop_na()
johnmod <- glmmTMB(state_int ~ class +pland + pd + ed + shape_mn + (1|tag), data = johnd2, family = binomial)
summary(johnmod)
performance(johnmod)
check_model(johnmod)

johnmod <- glmmTMB(state_int ~ class +pland + ed + (1|tag), data = johnd2, family = binomial)
summary(johnmod)
performance(johnmod)
check_model(johnmod)
performance::check_overdispersion(johnmod)
johnmodgam <- gam(state_int ~ class + s(pland) + s(pd) + s(ed) + s(shape_mn) + s(tag, bs="re"), 
    family = binomial, data = johnd2)
summary(johnmodgam)
check_model(johnmodgam)
performance(johnmodgam)

#GAMMs do not have a great outcome. GLMMs are not working. I am going to try RF and see. It could be that there is no relationship
#do it for both, Rankin is first
rankd <- read_csv('Results/Rankin_full_hmm.csv')
johnd <- read_csv('Results/Johnson_full_hmm.csv')
rankd$state <- as.factor(rankd$state)
rankd$class <- as.factor(rankd$class)
johnd$state <- as.factor(johnd$state)
johnd$class <- as.factor(johnd$class)

#randomly select 70% of state 1 and state 2 for training
head(rankd)
rankd$ID <- 1:nrow(rankd)
set.seed(1919)
rsfr1 <- rankd %>% dplyr::filter(state == '1') %>% slice_sample(prop = 0.7)
rsfr2 <- rankd %>% dplyr::filter(state == '2') %>% slice_sample(prop = 0.7)

rsf.train <- rbind(rsfr1,rsfr2)
rsf.test <- rankd[!(rankd$ID %in%rsf.train$ID), ] 

#take out ID column
rsf.testr <- rsf.test %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class)) %>% drop_na() 
rsf.trainr <- rsf.train %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class)) %>% drop_na()


# Set tasks for training and test datasets.
task_trout.train <- as_task_classif(
  rsf.trainr, target = "state", positive = '2',
)
str(task_trout.train)


task_trout.test <- as_task_classif(
  x = rsf.testr, target = "state",
  positive = "2"
)

# Make learner.
learner <-lrn(
  "classif.ranger",
  predict_type = "prob",
  mtry  = to_tune(1, ncol(rsf.trainr) - 1),
  sample.fraction = to_tune(0.2, 0.9),
  min.node.size = to_tune(p_int(1, 10)),
  importance = 'impurity'
)

#tune hyperparameters
rm(instance)
instance = mlr3tuning::ti(
  task = task_trout.train,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measures = msr("classif.ce"),
  terminator = trm("none")
)

tuner = mlr3tuning::tnr("grid_search", resolution = 5, batch_size = 5)
tuner$optimize(instance) # Takes ~ 4 minutes on my relatively fast computer

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
x_trout <- rsf.trainr %>% dplyr::select(-state) 
# Create "Predictor" object to interpret findings via the iml package.
predictor_trout <- Predictor$new(learner, data = x_trout, y = rsf.trainr$state) 
options(future.globals.maxSize = 1024 * 1024 * 1024) 
imp_trout <- FeatureImp$new(predictor_trout, loss = "ce") # Calculate importance.
x = c('darkgreen', 'forestgreen', 'green3', 'lightgreen', 'lightblue', 'blue3')
p <- imp_trout$plot()

# Remove all point layers
p$layers <- Filter(function(layer) {
  !inherits(layer$geom, c("GeomPoint", "GeomErrorbar", "GeomErrorbarh", 'GeomSegment'))
}, p$layers)

# Add your clean bars
p +
  geom_bar(stat = "identity", fill = x) +
  theme_classic()

#ed is the best predictor, so let's look individually
effect_ed <- FeatureEffect$new(predictor_trout, feature = c('ed'), method = 'pdp')
beep(1)
effect_ed$plot()
effect_ed
p_res <- effect_ed$results
p_res <- p_res %>% filter(.class == 2)

ggplot(p_res, aes(x = ed, y = .value))+
  geom_smooth(color = 'darkgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Edge Density",
       y = 'Probability of detection')+
  theme_classic()

#not super ecologically significant. Let's look at shape_mn
effect_shape <- FeatureEffect$new(predictor_trout, feature = c('shape_mn'), method = 'pdp')
effect_shape$plot()


sh_res <- effect_shape$results
sh_res <- sh_res %>% filter(.class == 2)

ggplot(sh_res, aes(x = shape_mn, y = .value))+
  geom_smooth(color = 'green3', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Probability of detection')+
  theme_classic()

#again, not great. Lastly, pland
effect_pland <- FeatureEffect$new(predictor_trout, feature = c('pland'), method = 'pdp')
effect_pland$plot()


pl_res <- effect_pland$results
pl_res <- pl_res %>% filter(.class == 2)

ggplot(pl_res, aes(x = pland, y = .value))+
  geom_smooth(color = 'green3', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Percentage of Landscape",
       y = 'Probability of detection')+
  theme_classic()

#not great. So Rankin HMMs don't tell us anything about behavior vs. background habitat
#johnson
#randomly select 70% of state 1 and state 2 for training
head(johnd)
johnd$ID <- 1:nrow(johnd)
set.seed(1919)
rsfj1 <- johnd %>% dplyr::filter(state == '1') %>% slice_sample(prop = 0.7)
rsfj2 <- johnd %>% dplyr::filter(state == '2') %>% slice_sample(prop = 0.7)

rsf.train <- rbind(rsfj1,rsfj2)
rsf.test <- johnd[!(johnd$ID %in%rsf.train$ID), ] 

#take out ID column
rsf.testj <- rsf.test %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class)) %>% drop_na() 
rsf.trainj <- rsf.train %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class)) %>% drop_na()


# Set tasks for training and test datasets.
task_trout.train <- as_task_classif(
  rsf.trainj, target = "state", positive = '2',
)
str(task_trout.train)


task_trout.test <- as_task_classif(
  x = rsf.testj, target = "state",
  positive = "2"
)

# Make learner.
learner <-lrn(
  "classif.ranger",
  predict_type = "prob",
  mtry  = to_tune(1, ncol(rsf.trainr) - 1),
  sample.fraction = to_tune(0.2, 0.9),
  min.node.size = to_tune(p_int(1, 10)),
  importance = 'impurity'
)

#tune hyperparameters
rm(instance)
instance = mlr3tuning::ti(
  task = task_trout.train,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measures = msr("classif.ce"),
  terminator = trm("none")
)

tuner = mlr3tuning::tnr("grid_search", resolution = 5, batch_size = 5)
tuner$optimize(instance) # Takes ~ 4 minutes on my relatively fast computer

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
x_trout <- rsf.trainj %>% dplyr::select(-state) 
# Create "Predictor" object to interpret findings via the iml package.
predictor_trout <- Predictor$new(learner, data = x_trout, y = rsf.trainj$state) 
options(future.globals.maxSize = 1024 * 1024 * 1024) 
imp_trout <- FeatureImp$new(predictor_trout, loss = "ce") # Calculate importance.
x = c('darkgreen', 'forestgreen', 'green3', 'lightgreen', 'lightblue', 'blue3')
p <- imp_trout$plot()

# Remove all point layers
p$layers <- Filter(function(layer) {
  !inherits(layer$geom, c("GeomPoint", "GeomErrorbar", "GeomErrorbarh", 'GeomSegment'))
}, p$layers)

# Add your clean bars
p +
  geom_bar(stat = "identity", fill = x) +
  theme_classic()

#shape most important
effect_shape <- FeatureEffect$new(predictor_trout, feature = c('shape_mn'), method = 'pdp')
effect_shape$plot()


sh_res <- effect_shape$results
sh_res <- sh_res %>% filter(.class == 2)

ggplot(sh_res, aes(x = shape_mn, y = .value))+
  geom_smooth(color = 'green3', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Mean Shape",
       y = 'Probability of detection')+
  theme_classic()

#still not great. Let's try second important factor: ed
effect_ed <- FeatureEffect$new(predictor_trout, feature = c('ed'), method = 'pdp')
beep(1)
effect_ed$plot()
effect_ed
p_res <- effect_ed$results
p_res <- p_res %>% filter(.class == 2)

ggplot(p_res, aes(x = ed, y = .value))+
  geom_smooth(color = 'darkgreen', size = 2) +
  #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = .2)+
  labs(x = "Edge Density",
       y = 'Probability of detection')+
  theme_classic()
 #looks decent. Let's try a grouping to make sure this is real
johndj <- johnd %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class,tag)) %>% drop_na() 
rsf.testjg <- rsf.test %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class,tag)) %>% drop_na() 
rsf.trainjg <- rsf.train %>% as.data.frame %>% dplyr::select(c(division, ed, pd, pland, shape_mn, state, class,tag)) %>% drop_na()
rsf.trainjg$tag <- as.factor(rsf.trainjg$tag)
rsf.testjg$tag <- as.factor(rsf.testjg$tag)
task <- TaskClassif$new(
  id = "behavior_rf",
  backend =johndj,
  target = "state",
  positive = "2" # if your binary response is coded 0/1
)

set.seed(1919)  # for reproducibility

# Get unique tags
johndj$tag <- as.factor(johndj$tag)
unique_tags <- unique(johndj$tag)

# Assign each unique tag to one of 5 folds
folds <- sample(rep(1:5, length.out = length(unique_tags)))

# Map fold assignment back to all rows in your data
fold_ids <- setNames(folds, unique_tags)
johndj$fold <- fold_ids[johndj$tag]

resampling <- rsmp("custom")

resampling$instantiate(
  task,
  train_sets = lapply(1:5, function(i) which(johndj$fold != i)),
  test_sets  = lapply(1:5, function(i) which(johndj$fold == i))
)

learner <-lrn(
  "classif.ranger",
  predict_type = "prob",
  mtry  = to_tune(1, ncol(johndj) - 2),
  sample.fraction = to_tune(0.2, 0.9),
  min.node.size = to_tune(p_int(1, 10)),
  importance = 'impurity'
)

#tune hyperparameters
rm(instance)
instance = mlr3tuning::ti(
  task = task,
  learner = learner,
  resampling = rsmp("cv", folds = 5),
  measures = msr("classif.ce"),
  terminator = trm("none")
)

tuner = mlr3tuning::tnr("grid_search", resolution = 5, batch_size = 5)
tuner$optimize(instance) # Takes ~ 4 minutes on my relatively fast computer

#store tuned hyperparameters in learner
learner$param_set$values <- instance$result_learner_param_vals

rr <- resample(task, learner, resampling, store_models = TRUE)


acc_blocked <- rr$aggregate(msr("classif.acc"))
acc_blocked
#56.6% accuracy. Not good, shows there is autocorrelation. I don't think we should present these results