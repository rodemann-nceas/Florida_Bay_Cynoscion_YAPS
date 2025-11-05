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
write.csv(rankd, file = 'Results/Rankin_full_hmm.csv')
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
write.csv(johnd, file = 'Results/Johnson_full_hmm.csv')

#Let's try some modelling!

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
