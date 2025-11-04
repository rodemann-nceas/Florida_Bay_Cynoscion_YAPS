#landscape classification using machine learning - random forests
#author: Jon Rodemann
#date last edited: 7/1/2022
#Install and load in libraries
library(mlr)
library(tidyverse)
library(randomForest)
library(raster)
library(landscapemetrics)
library(future)
library(furrr)
library(sf)
library(terra)

#load in datasets
setwd('E:/FIU/Project/YAPS/Rankin_YAPS/')

map <- stack('aerialmap.tif')
#plot(map)

training_points_bare <- shapefile('Mapping/Training_points_bare.shp')
training_points_sparse <- shapefile('Mapping/Training_points_sparse.shp')
training_points_dense <- shapefile('Mapping/Training_points_dense.shp')

#extract values from map
tpb_ex <- raster::extract(map, training_points_bare, df = T)
#tpb_ex
tps_ex <-raster::extract(map, training_points_sparse, df = T) 
#tps_ex
tpd_ex <- raster::extract(map, training_points_dense, df = T)
#tpd_ex

#classify and combine training points
tpb_ex$train <- 'bare'
tps_ex$train <- 'sparse'
tpd_ex$train <- 'dense'

training <- rbind(tpb_ex, tps_ex, tpd_ex)
training

training <- training %>% dplyr::select(-ID) %>% dplyr::rename('band1' = aerialmap_1, 'band2' = aerialmap_2, 'band3' = aerialmap_3, 'band4' = aerialmap_4)
#str(training)

#set up machine learning
rforestLearner <- makeLearner('classif.randomForest')
SAVTask <- makeClassifTask(data = training, target = 'train')
SAVtrained <- train(rforestLearner, SAVTask)
SAVtrained

p <- predict(SAVtrained, newdata = training) #overfitting easily can happen, evaluating on training set
performance(p, measures = list(acc, mmce))

#CV
kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 10, reps = 10)
rforest_cv <- mlr::resample(rforestLearner, task = SAVTask, resampling = kfold_cv, measures = list(acc, mmce))
rforest_cv
#much better than decision tree!
#optimize hyperparameters
getParamSet(rforestLearner)

rforestparam <- makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = ncol(training) - 1),
  makeIntegerParam("nodesize", lower = 1, upper = 10),
  makeIntegerParam('ntree', lower = 250, upper = 750)
)

randSearch <- makeTuneControlRandom(maxit = 50)
rforest_tuned <- tuneParams('classif.randomForest', 
                            SAVTask,
                            resampling = kfold_cv,
                            par.set = rforestparam,
                            control = randSearch)

rforest2 <- setHyperPars(rforestLearner, par.vals = rforest_tuned$x)
rforest2_trained <- train(rforest2, SAVTask)
rforest2_cv <- mlr::resample(rforest2, task = SAVTask, resampling = kfold_cv, measures = list(acc, mmce))

rforest2_trained
getFeatureImportance(rforest2_trained)
rforest2_model <- getLearnerModel(rforest2_trained)

#the tuned model is better! apply it across the image
#rename bands
names(map) <- c('band1', 'band2', 'band3', 'band4')
SAVmap <- raster::predict(map, rforest2_model)
plot(SAVmap)

writeRaster(SAVmap, filename = 'savclass_2.tif')

#create classification matrix for new map
reclass_df <- c(0, 1, 1,
                1, 2, 3,
                2, 3,2)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
SAVmap1 <- reclassify(SAVmap,
                             reclass_m)
plot(SAVmap1)
#resample to 1 meter resolution
mapag <- raster::aggregate(SAVmap1, fact = 3)
plot(mapag)

reclass_df <- c(0, 1, 1,
                1, 2, 2,
                2, 3,3)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
mapag1 <- reclassify(mapag,
                      reclass_m)
plot(mapag1)
check_landscape(mapag1)

#take out only sparse, then dense classes
reclass_df <- c(0, 1, 0,
                1, 2, 0,
                2, 3, 1)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
map_dense <- reclassify(mapag1,
                     reclass_m)
plot(map_dense)
check_landscape(map_dense)

reclass_df <- c(0, 1, NA,
                1, 2, 1,
                2, 3, NA)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
map_sparse <- reclassify(mapag1,
                        reclass_m)
plot(map_sparse)


#now we can run landscape metrics on them!
#since our positions are based on our current GPS, landscape metrics can't be more fine than 5 meters
#start with dense - 5 meter matrix
# moving_window <- matrix(1, nrow = 5, ncol = 5)
dim(map_dense)
#let's run it and hope...
future::plan(multisession, workers = 9)
coord <- data.frame(coordinates(map_dense))
coord <- as.data.frame(cbind(dat_10ang$x, dat_10ang$y))

re <- as(extent(map), "SpatialPolygons")
crs(re) <- crs(map)  # Assign CRS
resf <- st_as_sf(re)  # Convert to sf object


coord1 <- coord_sf[st_within(coord_sf, resf, sparse = FALSE), ]
d10 <- d10[st_within(d10, resf, sparse = FALSE), ]

coord_sf <- st_as_sf(coord, coords = c("V1", "V2"), crs = "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs")

st_bbox(map_dense)
st_bbox(coord_sf)

print(class(coords_sf))       # Should include "sf"
print(st_geometry(coords_sf))

if (st_crs(coords_sf) != st_crs(map_dense)) {
  coords_sf <- st_transform(coords_sf, crs = st_crs(map_dense))
}


result_test <- sample_lsm(
  map_dense,
  y = coord_sf[1, ],
  size = 5,
  what = c(
    "lsm_c_pd", "lsm_c_ed", "lsm_c_area_mn", "lsm_c_ai",
    "lsm_c_para_mn", "lsm_c_gyrate_mn", "lsm_c_division",
    "lsm_c_pland", "lsm_c_shape_mn"
  ),
  shape = "circle"
)
1:nrow(coord_sf)
# Apply the function

library(furrr)
plan(multisession)  # Use multisession to run on multiple cores
library(purrr)
coord[i,]
result_met <- 1:nrow(coord1) %>% future_map(
  function(i) {
    sample_lsm(
      map_dense,
      y = st_geometry(coord1[i, ]),  # Ensure geometry is passed
      size = 5,
      what = c(
        "lsm_c_pd", "lsm_c_ed", "lsm_c_area_mn", "lsm_c_ai",
        "lsm_c_para_mn", "lsm_c_gyrate_mn", "lsm_c_division",
        "lsm_c_pland", "lsm_c_shape_mn"
      ),
      shape = "circle"
    )
  }
)

# coord_sf[i, ] %>% class()
# ?future_map
# 
# future::plan(multisession, workers = 9)
# result_met <- furrr::future_map(1:nrow(coor), function(i)(sample_lsm(map_dense, y = coor[i, ], size = 5, what = c("lsm_c_pd", "lsm_c_ed", 'lsm_c_area_mn', 'lsm_c_ai', 'lsm_c_para_mn', 'lsm_c_gyrate_mn', 'lsm_c_division', 'lsm_c_pland', 'lsm_c_shape_mn'), shape = 'circle')))
# future::plan(sequential)

rt <-lapply(result_met, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})

problematic <- which(sapply(rt, function(df) !inherits(try(bind_rows(df), silent = TRUE), "data.frame")))
## or
## is.na(Result) <- !lengths(Result)


met <- bind_rows(rt, .id = 'source')

met_slice <- met %>% slice(1:1000)

unique_sources <- unique(met$source)
unique_classes <- unique(met$class)
unique_metric <- unique(met$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(met, by = c("source", "class", 'metric')) %>% distinct()

write.csv(file = 'metslice.csv', met_slice)
nrow(as.data.frame(unique(met$source)))
met_1 <- mettry %>% dplyr::filter(class == 1) %>% mutate(source = as.numeric(source)) %>% arrange(source)
met_2 <- dplyr::bind_rows(result_met1)
met_3 <- met_2 %>% dplyr::filter(class == 1)
rep8<-function(d){
  d$newcol =  rep(1:ceiling(nrow(d)/8), each = 8)[1:nrow(d)]
  d
}

rep9<-function(d){
  d$newcol =  rep(1:ceiling(nrow(d)/9), each = 9)[1:nrow(d)]
  d
}
rep1<-function(d){
  d$newcol =  rep(1:ceiling(nrow(d)/1), each = 1)[1:nrow(d)]
  d
}


met_1 <- rep9(met_1)
#met_2 <- rep1(met_3)
met_1 <- met_1 %>% dplyr::select(-plot_id) %>% dplyr::rename(plot_id = newcol)

met_div <- met_1 %>% dplyr::filter(metric == 'division') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('div' = value)
met_pland <- met_1 %>% dplyr::filter(metric == 'pland') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pland' = value)
met_ai <- met_1 %>% dplyr::filter(metric == 'ai') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('ai' = value)
met_pd <- met_1 %>% dplyr::filter(metric == 'pd') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pd' = value)
met_ed <- met_1 %>% dplyr::filter(metric == 'ed') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('ed' = value)
met_para_mn <- met_1 %>% dplyr::filter(metric == 'para_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('para_mn' = value)
met_gyrate <- met_1 %>% dplyr::filter(metric == 'gyrate_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('gyrate_mn' = value)
#met_cohesion <- met_1 %>% dplyr::filter(metric == 'cohesion') %>% dplyr::select(c(value, plot_id)) %>% rename('cohesion' = value)
#met_np <- met_1 %>% dplyr::filter(metric == 'np') %>% dplyr::select(c(value, plot_id)) %>% rename('np' = value)
#met_te <- met_1 %>% dplyr::filter(metric == 'te') %>% dplyr::select(c(value, plot_id)) %>% rename('te' = value)
#met_ca <- met_1 %>% dplyr::filter(metric == 'ca') %>% dplyr::select(c(value, plot_id)) %>% rename('ca' = value)
met_area_mn <- met_1 %>% dplyr::filter(metric == 'area_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('area_mn' = value)
met_shape_mn <- met_1 %>% dplyr::filter(metric == 'shape_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('shape_mn' = value)
#met_enn_mn <- met_1 %>% dplyr::filter(metric == 'enn_mn') %>% dplyr::select(c(value, plot_id)) %>% rename('enn_mn' = value)
# met_clump <- met_2 %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('clumpy' = value)
# met_pladj <- met_2 %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pladj' = value)

met_all <- list(met_div, met_pland, met_area_mn, met_pd, met_ed, met_para_mn, met_gyrate, met_shape_mn, met_ai)
#met_all1 <- list(met_pd, met_ed, met_np, met_te, met_ca, met_area_mn, met_shape_mn, met_div)
met_all <- met_all %>% purrr::reduce(full_join, by='plot_id')
#met_all1 <- met_all1 %>% purrr::reduce(full_join, by='plot_id')
#met_all2 <- list(met_pland, met_pd, met_ed, met_shape_mn, met_para_mn)
#met_all2 <- met_all2 %>% purrr::reduce(full_join, by='plot_id')

head(met_all)

write.csv(file = 'met_all.csv', met_all)

########RSF data frame#############
#rename id column in met_all
met_all <- met_all %>% dplyr::rename(ID = plot_id) %>% tidyr::drop_na()

##Now with Johnson data!
setwd('E:/FIU/PostDoc/CESI/Final_report/YAPS/')

image <- raster::stack('20230730_150134_07_2460_3B_AnalyticMS_SR_8b_harmonized_clip.tif')
mask <- shapefile('John_mask.shp')

#crop image
map1 <- raster::crop(image, mask)

#add in training data
td <- read.csv('Data_Johnson_YAPS_GCP_5-2023.csv')

#add in dense, sparse, bare
for (i in 1:nrow(td)) {
  if (td$Tot[i] <= 20) {
    td$class[i] <- 0
  } else if (td$Tot[i] > 20 & td$Tot[i] <= 40) {
    td$class[i] <- 1
  } else {
    td$class[i] <- 2
  }
}

#make into spatial data frame
sftd <- td %>% 
  st_as_sf(coords = c('UTM.Easting', 'UTM.Northing'), crs = 2958)

tpd <- extract(map1, sftd)

tpd <- cbind(tpd, td)

tpd$class <- as.factor(tpd$class)

head(tpd)

training <- tpd %>% dplyr::select(-c(Point, Date, Tot, Tt, Hw, Sf, Drift, GC, CH, UTM.Easting, UTM.Northing)) 
training2 <- tpd %>% dplyr::select(-c(Point, Date, Tt, Hw, Sf, Drift, GC, CH, UTM.Easting, UTM.Northing, class))
#str(training)

#set up machine learning
rforestLearner <- makeLearner('regr.randomForest')
SAVTask <- makeRegrTask(data = training2, target = 'Tot')
SAVtrained <- train(rforestLearner, SAVTask)
SAVtrained

p <- predict(SAVtrained, newdata = training2) #overfitting easily can happen, evaluating on training set
performance(p, measures = list(acc, mmce))

#CV
kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 10, reps = 10)
rforest_cv <- mlr::resample(rforestLearner, task = SAVTask, resampling = kfold_cv, measures = list(acc, mmce))
rforest_cv
#much better than decision tree!
#optimize hyperparameters
getParamSet(rforestLearner)

rforestparam <- makeParamSet(
  makeIntegerParam("mtry", lower = 1, upper = ncol(training) - 1),
  makeIntegerParam("nodesize", lower = 1, upper = 10),
  makeIntegerParam('ntree', lower = 250, upper = 750)
)

randSearch <- makeTuneControlRandom(maxit = 50)
rforest_tuned <- tuneParams('regr.randomForest', 
                            SAVTask,
                            resampling = kfold_cv,
                            par.set = rforestparam,
                            control = randSearch)

rforest2 <- setHyperPars(rforestLearner, par.vals = rforest_tuned$x)
rforest2_trained <- train(rforest2, SAVTask)
rforest2_cv <- mlr::resample(rforest2, task = SAVTask, resampling = kfold_cv, measures = list(acc, mmce))

rforest2_trained
getFeatureImportance(rforest2_trained)
rforest2_model <- getLearnerModel(rforest2_trained)
rforest2_modelsav <- getLearnerModel(SAVtrained)

SAVmap1 <- raster::predict(map1, rforest2_model)
plot(SAVmap1)

writeRaster(SAVmap1, filename = 'johnclass2.tif', overwrite = T)

setwd('E:/FIU/PostDoc/CESI/Final_report/YAPS/')

jmap <- rast('johnclass2.tif')
plot(jmap)

reclass_df <- c(0, 44.99999999, 0,
                45, 100, 1)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
jmap1 <- classify(jmap,
                      reclass_m)
plot(jmap1)
writeRaster(jmap1, 'jmap_dense.tif', overwrite = T)
dim(jmap1)
#let's run it and hope...
dat_john <- read.csv('Seatrout_positions_sd10.csv')

coord <- as.data.frame(cbind(dat_john$x, dat_john$y))
coord_sf <- st_as_sf(coord, coords = c("V1", "V2"), crs = "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs")

jmap1 <- rast('jmap_dense.tif')
crs(jmap1)
crs(coord_sf)
if (st_crs(coord_sf) != st_crs(jmap1)) {
  coord_sf <- st_transform(coord_sf, crs = st_crs(jmap1))
}
# Get extent of the raster as an sf polygon
raster_bbox <- as.polygons(jmap1) |> st_as_sf()
crs(raster_bbox)
raster_bbox$merged_class <- 1  # assign same class
merged <- st_union(raster_bbox)
plot(merged)
crs(merged)
# Subset points inside the raster extent
coord_sf_in <- coord_sf[st_within(coord_sf, merged, sparse = FALSE), ]

plot(jmap1)
points(coord_sf_in)

# Apply the function

library(furrr)
plan(multisession)  # Use multisession to run on multiple cores
library(purrr)
coord[i,]
result_met <- 1:nrow(coord_sf_in) %>% future_map(
  function(i) {
    jmap_r <- rast("jmap_dense.tif")
    sample_lsm(
      jmap_r,
      y = st_geometry(coord_sf_in[i, ]),  # Ensure geometry is passed
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
# coord_sf[i, ] %>% class()
# ?future_map
# 
# future::plan(multisession, workers = 9)
# result_met <- furrr::future_map(1:nrow(coor), function(i)(sample_lsm(map_dense, y = coor[i, ], size = 5, what = c("lsm_c_pd", "lsm_c_ed", 'lsm_c_area_mn', 'lsm_c_ai', 'lsm_c_para_mn', 'lsm_c_gyrate_mn', 'lsm_c_division', 'lsm_c_pland', 'lsm_c_shape_mn'), shape = 'circle')))
# future::plan(sequential)
# 
rt <-lapply(result_met, as.data.frame)

rt<- lapply(rt, function(x) if (is.null(x)) data.frame() else x)

all_columns <- unique(unlist(lapply(rt, colnames)))
rt <- lapply(rt, function(df) {
  df[setdiff(all_columns, colnames(df))] <- NA  # Add missing columns with NA
  df <- df[, all_columns]  # Reorder columns
  df
})

problematic <- which(sapply(rt, function(df) !inherits(try(bind_rows(df), silent = TRUE), "data.frame")))
## or
## is.na(Result) <- !lengths(Result)


met <- bind_rows(rt, .id = 'source')

met_slice <- met %>% slice(1:1000)

unique_sources <- unique(met$source)
unique_classes <- unique(met$class)
unique_metric <- unique(met$metric)



# Create a complete grid of all combinations of `source` and `class`
full_combinations <- expand.grid(source = unique_sources, class = unique_classes, metric = unique_metric)

# Perform a left join with the original dataset to add NAs for missing combinations
mettry <- full_combinations %>%
  left_join(met, by = c("source", "class", 'metric')) %>% distinct()

mettry1 <- mettry %>% dplyr::filter(class == 1) %>% pivot_wider(names_from = metric, values_from = value)


dat_john_s <- st_as_sf(dat_john, coords = c("x", "y"), crs = "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs")

if (st_crs(dat_john_s) != st_crs(jmap1)) {
  dat_john_s <- st_transform(dat_john_s, crs = st_crs(jmap1))
}

dat_john_s1 <- dat_john_s[st_within(dat_john_s, merged, sparse = FALSE), ]
dat_all <- as.data.frame(dat_john_s1)
names(mettry1)

metac <- mettry1 %>% dplyr::select(c(source, division, ed, pd, pland, shape_mn))

dat_all1 <- cbind(dat_all, metac)

classif <- extract(jmap1, coord_sf_in)
classif <- classif %>% dplyr::select(-ID)
dat_all2 <- cbind(dat_all1, classif)

names(dat_all1)
write.csv(dat_all, file = 'E:/FIU/PostDoc/CESI/Final_report/YAPS/johndat_yaps.csv')


# Bind coordinates back to the original data

dat_10 <- dat_all2
dat_10$track <- NA

str(dat_10)

dat_10$ts <- as.POSIXct(dat_10$ts, format='%Y-%m-%dT%H:%M:%OS', tz='America/New_York', origin = '1970-01-01 00:00:00')

dat_10t <- dat_10 %>% group_by(tag) %>% mutate(trackstart = ifelse(as.numeric(unclass(ts)-unclass(lag(ts))) > 3600, 'Y', 'N')) %>% ungroup()
dat_10t[is.na(dat_10t)] <- 'Y'

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

dat_10ttt <- dat_10tt %>% group_by(tag) %>% group_modify(~tracknum1(.x)) %>% ungroup() %>% dplyr::select(-X)
coords_dat <- st_coordinates(dat_10ttt$geometry)

# Combine coordinates with attributes
dat_10cor <- cbind(dat_10ttt, coords_dat) 
#Add step length and turn angle
dat_10ang2 <- dat_10cor %>% group_by(tag, track) %>% mutate(step_length = sqrt((X - lag(X, default = NA))^2 + (Y - lag(Y, default = NA))^2), bear = atan2(
  Y - lag(Y, default = NA),
  X - lag(X, default = NA)
) * (180 / pi), bear_utm = ifelse(bear< 0, bear + 360, bear), angle = lag(bear_utm, default = NA)-bear_utm) %>% ungroup()

dat_e <- extract(jmap, data.frame(dat_10ang2$X, dat_10ang2$Y))

dat_allfin <- cbind(dat_10ang2, dat_e) %>% dplyr::select(-ID)

write.csv(dat_10ang2, file = 'E:/FIU/PostDoc/CESI/Final_report/YAPS/jdat_final.csv')

# write.csv(file = 'metslice.csv', met_slice)
# nrow(as.data.frame(unique(met$source)))
# met_1 <- mettry %>% dplyr::filter(class == 1) %>% mutate(source = as.numeric(source)) %>% arrange(source)
# met_2 <- dplyr::bind_rows(result_met1)
# met_3 <- met_2 %>% dplyr::filter(class == 1)
# rep8<-function(d){
#   d$newcol =  rep(1:ceiling(nrow(d)/8), each = 8)[1:nrow(d)]
#   d
# }
# 
# rep9<-function(d){
#   d$newcol =  rep(1:ceiling(nrow(d)/9), each = 9)[1:nrow(d)]
#   d
# }
# rep1<-function(d){
#   d$newcol =  rep(1:ceiling(nrow(d)/1), each = 1)[1:nrow(d)]
#   d
# }
# 
# 
# met_1 <- rep9(met_1)
# #met_2 <- rep1(met_3)
# met_1 <- met_1 %>% dplyr::select(-plot_id) %>% dplyr::rename(plot_id = newcol)
# 
# met_div <- met_1 %>% dplyr::filter(metric == 'division') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('div' = value)
# met_pland <- met_1 %>% dplyr::filter(metric == 'pland') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pland' = value)
# met_ai <- met_1 %>% dplyr::filter(metric == 'ai') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('ai' = value)
# met_pd <- met_1 %>% dplyr::filter(metric == 'pd') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pd' = value)
# met_ed <- met_1 %>% dplyr::filter(metric == 'ed') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('ed' = value)
# met_para_mn <- met_1 %>% dplyr::filter(metric == 'para_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('para_mn' = value)
# met_gyrate <- met_1 %>% dplyr::filter(metric == 'gyrate_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('gyrate_mn' = value)
# #met_cohesion <- met_1 %>% dplyr::filter(metric == 'cohesion') %>% dplyr::select(c(value, plot_id)) %>% rename('cohesion' = value)
# #met_np <- met_1 %>% dplyr::filter(metric == 'np') %>% dplyr::select(c(value, plot_id)) %>% rename('np' = value)
# #met_te <- met_1 %>% dplyr::filter(metric == 'te') %>% dplyr::select(c(value, plot_id)) %>% rename('te' = value)
# #met_ca <- met_1 %>% dplyr::filter(metric == 'ca') %>% dplyr::select(c(value, plot_id)) %>% rename('ca' = value)
# met_area_mn <- met_1 %>% dplyr::filter(metric == 'area_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('area_mn' = value)
# met_shape_mn <- met_1 %>% dplyr::filter(metric == 'shape_mn') %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('shape_mn' = value)
# #met_enn_mn <- met_1 %>% dplyr::filter(metric == 'enn_mn') %>% dplyr::select(c(value, plot_id)) %>% rename('enn_mn' = value)
# # met_clump <- met_2 %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('clumpy' = value)
# # met_pladj <- met_2 %>% dplyr::select(c(value, plot_id)) %>% dplyr::rename('pladj' = value)
# 
# met_all <- list(met_div, met_pland, met_area_mn, met_pd, met_ed, met_para_mn, met_gyrate, met_shape_mn, met_ai)
# #met_all1 <- list(met_pd, met_ed, met_np, met_te, met_ca, met_area_mn, met_shape_mn, met_div)
# met_all <- met_all %>% purrr::reduce(full_join, by='plot_id')
# #met_all1 <- met_all1 %>% purrr::reduce(full_join, by='plot_id')
# #met_all2 <- list(met_pland, met_pd, met_ed, met_shape_mn, met_para_mn)
# #met_all2 <- met_all2 %>% purrr::reduce(full_join, by='plot_id')
# 
# head(met_all)
