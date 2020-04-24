#############################################################
##Date: 4/2/2019
##Purpose: Vet ST-GPR input data in response to plots from run 45750
##Notes: 
##Updates: Check out issues with variance and outliering
# Do I drop as outlier or just drop the data point? possibly drop data point because I've heard weird stuff with outliers...
###########################################################
library(data.table)
library(readr)
source(get_outputs.R)
source('get_envelope.R')
source('get_draws.R')
source('get_location_metadata.R')
source('get_cause_metadata.R')
source('get_demographics.R')
source('utility.r')
locs <- get_location_metadata(35)

source('get_covariate_estimates.R')
test <- get_covariate_estimates(1099, gbd_round_id = 7, decomp_step = 'iterative')

## Set up in-facility delivery
# ifd <- get_covariate_estimates(51, location_id = unique(dt$location_id), gbd_round_id = 6, decomp_step = 'step1', sex_id = c(3), year_id = c(1990:2019))
# ## Update in-facility delivery
# # Do at draw level
# ifd_draw_files <- Sys.glob(FILEPATH)
# ifd_draws <- rbindlist(lapply(ifd_draw_files, function(x){
#   ifd_draw1 <- fread(x) %>%
#     melt(id.vars = c('run_id', 'covariate_id', 'location_id', 'year_id', 'age_group_id', 'sex_id'), value.vars = patterns('draw_'))
#   # SE calculation from this formula: sd(draws)/sqrt(1000)
#   ifd_draw1[, se := sd(value), by = c('location_id', 'year_id')]
#   ifd_draw1[, se := se/sqrt(1000)]
#   ifd_draw1[, mean := mean(value), by = c('location_id', 'year_id')]
#   ifd_draw1 <- ifd_draw1[year_id >= 1990, c('location_id', 'year_id', 'mean', 'se')] %>% unique()
#   return(ifd_draw1)
# }))
# setnames(ifd_draws, c('mean', 'se'), c('data', 'standard_error'))
# ifd_draws[, variance := standard_error^2]
# 
# # Double by sexes
# ifd_draws[, index := 1]
# sexes <- data.table(sex_id = c(1,2), index = 1)
# ifd_draws <- merge(ifd_draws, sexes, by = 'index', allow.cartesian = TRUE)
# ifd_draws[, me_name := 'ip_envelope'][, measure := 'continuous'][, is_outlier := 0][, age_group_id := 2]
# 
# ifd_nid <- 310156
# dt <- dt[nid != ifd_nid]
# ifd_draws[, nid := ifd_nid]
# dt <- rbind(dt, ifd_draws, fill = TRUE) # age cols and sample size are NA. shouldn't matter
# # 

# # Poland shows I think I used population as sample size for some. Going to check
# all_pops <- get_population(location_id = unique(dt$location_id), age_group_id = unique(dt$age_group_id), year_id = c(1990:2019), sex_id = c(1,2), gbd_round_id = 6, decomp_step = 'step1')
# dt <- merge(dt, all_pops, by = c('age_group_id', 'location_id', 'sex_id', 'year_id'))
# dt[sample_size > population, sample_size := population]
# # Change age group ID 1 because I think this is an error on my part for Austria data only
# dt[age_group_id == 1, age_group_id := 5]
# 
# # Outlier Armenia and Fujian and Gansu and Guangdong and Guangxi and Guizhou and Hainan and Hebei data point (Chinese all look like one source)
# dt <- dt[nid != 307778]
# dt <- dt[nid != 138595]
# 
# ##### Some outliers form Spencer ####
# # Kenya 2003 and VNM 2002
# dt[location_id ==180 & year_id == 2003 & nid != 310156, is_outlier := 1]
# dt[location_id == 20 & year_id == 2002 & nid != 310156, is_outlier := 1]

date <- Sys.Date()
write_csv(dt, paste0(FILEPATH))

# Previous model
#write_csv(dt, FILEPATH)

## Add all-cause mortality as custom covariate. Use get_outputs square to see if you can do it that way because the get_draws didn't really work
#dt <- fread(FILEPATH)


ages <- c(2,3,28, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,235,238,388,389)
ages <- sort(ages)
locs <- get_location_metadata(22, gbd_round_id = 6)
#locs <- locs[level >= 3]

library(mortdb)
#mort <- mortdb::get_mort_outputs(model_name = "with shock death number", model_type = "estimate", gbd_year = 2020, run_id = "best")
#mort <- mort[sex_id %in% c(1,2) & year_id %in% c(1990:2019) & location_id %in% locs$location_id & 
#               age_group_id %in% ages]

mort <- get_envelope(year_id=c(1990:2019), sex_id=c(1,2), gbd_round_id=7, with_hiv=1, rates=1,
                     decomp_step = 'iterative', location_id = unique(locs$location_id),
                     age_group_id= ages)
setnames(mort, 'mean', 'cv_mort')
write_csv(mort, FILEPATH)

dt <- merge(dt, mort[, c('location_id', 'age_group_id', 'sex_id', 'year_id', 'cv_mort')], 
            by= c('location_id', 'age_group_id', 'sex_id', 'year_id'), all.x = TRUE)



