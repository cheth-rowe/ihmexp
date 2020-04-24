#############################################################
##Author: Nick Roberts
##Date: 3/12/2019
##Purpose: Make the hospital data predictions dataset
##Notes: 
###########################################################

user <- Sys.info()[7][[1]]
library(data.table)
library(parallel)
library(readstata13)
library(magrittr)
source(paste0('/share/code/hospital/', user, '/clinical_info/clinical_info/Functions/db_utilities.R'))
source('/ihme/cc_resources/libraries/current/r/get_age_metadata.R')
source('/ihme/cc_resources/libraries/current/r/get_location_metadata.R')
source('/ihme/cc_resources/libraries/current/r/get_covariate_estimates.R')
source('/ihme/cc_resources/libraries/current/r/get_population.R')
source('/ihme/cc_resources/libraries/current/r/get_model_results.R')

ages <- get_age_metadata(12, gbd_round_id = 5) %>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
age_1 <- data.table(age_group_id = 28, age_group_years_start = 0, age_group_years_end = 1)%>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
ages_over1 <- ages[age_group_id > 4] %>% rbind(age_1, fill = TRUE)

locs <- get_hosp_countries() %>% as.data.table()
locs[, index := 1]
ages_over1[, index := 1]

preddf <- merge(locs, ages_over1, by = 'index', allow.cartesian = TRUE)

###### Useful functions ####################################
# Not a db utility
interp_under1_env <- function(env_df){
  ages <- get_age_metadata(12, gbd_round_id = 5) %>% 
    setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
  env_df <- env_df[, c('location_id', 'year_id', 'age_group_id', 'sex_id', 'envelope_mean')]
  env_df <- merge(env_df, ages, by = 'age_group_id', all.x = TRUE)
  colnames(env_df)[5] <- 'envelope_mean'
  env_df <- env_df[age_group_id == 164, age_start := 0]
  env_df <- env_df[age_group_id == 164, age_end := 0]
  
  under_one_env <- env_df[(env_df$age_group_id == 2) |
                            (env_df$age_group_id == 3) |
                            (env_df$age_group_id == 4)]
  # get years that are present in under one envelope
  years <- unique(under_one_env$year_id)
  # get location_ids that are present in the envelope
  loc_list <- unique(under_one_env$location_id)
  # get population data for under one years old age groups
  under_one_pop <- get_population(age_group_id=c(2:4), location_id=loc_list,
                                  sex_id=c(1,2), year_id=years, gbd_round_id= 5)
  # under_one_pop = get_population(age_group_id=[2, 3, 4], location_id=-1,
  #                                sex_id=[1, 2], year_id=years)
  # Merge population data onto under one envelope
  under_one_env <- merge(under_one_env, under_one_pop, by = c('location_id', 'sex_id', 'age_group_id', 'year_id'))
  under_one_env$envelope_mean <- under_one_env$envelope_mean*under_one_env$population
  
  ## aggregate to just under_one age group
  under_one_env <- under_one_env[, lapply(.SD, sum), by = c('location_id', 'sex_id', 'year_id'),
                                 .SD = c('envelope_mean', 'population')]
  ## divide by pop to get actual value
  under_one_env$envelope_mean <- under_one_env$envelope_mean/under_one_env$population
  
  ## ad age group id
  under_one_env$age_group_id <- 28
  under_one_env$age_start <- 0
  under_one_env$age_end <- 1
  under_one_env$population <- NULL
  
  ## add back to envelope and replace age_group_ids 2, 3, 4
  env_df <- env_df[age_group_id != 2 & age_group_id != 3 & age_group_id != 4]
  env_df[, age_group_weight_value := NULL]
  env_df <- rbind(env_df, under_one_env)
}

preddf <- fread('/home/j/temp/nrober75/hospital/cf_models/cf_1_2/hospital_countries_inputs.csv')

### Get outpatient envelope ####
# Lots of deprecated code, keepign it in case we need to re-run again with new sources. It's just v slow.
# These are point estimate results, need to be interpolated unfortunately
# I interpolated in another script becaues it takes ages, saved it to a csv
# op_env <- get_model_results(gbd_team = 'epi', 
#                             gbd_round_id = 5, location_set_id = 35,
#                             model_version_id = 246617,
#                             age_group_id = c(2:20, 28, 30, 31, 32, 33, 164, 235))
# 
# op_env <- op_env[, c('location_id', 'year_id', 'age_group_id', 'sex_id', 'mean')] %>% setnames('mean', 'envelope_mean')
# op_env <- interp_under1_env(op_env)
# setnames(op_env, 'year_id', 'year_start')
# op_env <- op_env[location_id %in% unique(preddf$location_id)]

# Interpolate it to get new years, call child srcipt
# JUST KIDDING: THis was breaking. Ran it a single time in an interactive session with a single interpolate call
# takes a long time, will save and we're good to go
# locs <- get_location_metadata(35)
# loc_ids <- locs[level == 3]$location_id
# mclapply(loc_ids[1:5], function(i){
#   
#   system(paste0('qsub -j y -o /share/temp/sgeoutput/nrober75/errors -cwd -l archive -l m_mem_free=4G -l h_rt=10:00:00 -l fthread=4 -q all.q -P proj_hospital -N job_',loc_ids[i],
#                 ' /share/singularity-images/health_fin/forecasting/shells/health_fin_forecasting_shell_singularity.sh /share/code/hospital/nrober75/clinical_info/clinical_info/Corrections/mr_brt/envelope_interp/interp_env_child.R ',
#                 locs[i]))
#   
# }, mc.cores = 5)
# hosp_locs <- fread('/home/j/temp/nrober75/hosp_locations_years.csv')
# hosp_locs <- hosp_locs[location_id != 'NULL']
# setnames(hosp_locs, 'year_start_id', 'year_start')
# hosp_locs[, location_id := as.integer(location_id)]
# # df <- interpolate('modelable_entity_id', gbd_id = 19797, version_id = 246617, source = 'epi', measure_id = 19, sex_id = c(1,2), 
# #                   gbd_round_id= 5, reporting_year_start = 1990, location_id = hosp_locs$location_id,
# #                   reporting_year_end = 2017)
# op_env <- fread('/home/j/temp/nrober75/op_env_interp.csv')
# setnames(op_env, 'year_id', 'year_start')
# op_env <- merge(op_env, hosp_locs, by = c('location_id', 'year_start')) # Merge on hospital and year to get it to a more reasonable size
# op_env <- melt(op_env, id.vars = c('age_group_id', 'location_id', 'year_start', 'sex_id'), 
#            measure.vars = patterns('draw_')) 
# op_env[, envelope_mean := mean(value), by = c('age_group_id', 'location_id', 'year_start', 'sex_id')] # Get mean of the draws
# op_env[, value := NULL][, variable := NULL]
# setnames(op_env, 'year_start', 'year_id')
# op_env <- interp_under1_env(op_env)
# op_env <- unique(op_env) %>% setnames('year_id', 'year_start')
# op_env <- op_env[age_group_id != 164]

#write_csv(op_env, '/home/j/temp/hospital/2016/data/envelope/op_env_for_hosp_predictions.csv')
op_env <- fread('/home/j/temp/hospital/2016/data/envelope/op_env_for_hosp_predictions.csv')

## Inpatient enveelope
source("/ihme/code/st_gpr/central/r_functions/utilities/utility.r")
ages <- copy(ages_over1)
ip_envelope <- model_load(47570, 'gpr') %>% setnames('year_id', 'year_start') %>% merge(ages[, c('age_group_id', 'age_start')], by = 'age_group_id')
op_env <- merge(op_env, ip_envelope[, c('location_id', 'year_start', 'age_start', 'sex_id', 'gpr_mean')])
setnames(op_env, 'gpr_mean', 'ip_env')

# HAQ
haqi <- get_covariate_estimates(1099, gbd_round_id = 6, decomp_step = 'step1') %>% setnames(c('year_id', 'mean_value'), c('year_start', 'haqi'))
op_env <- merge(op_env, haqi[, c('location_id', 'year_start', 'haqi')], by = c('location_id', 'year_start'))
setnames(op_env, 'envelope_mean', 'op_env')

preddf <- copy(op_env)
preddf[, age_start := (age_start + age_end)/2] # NOt great, but good inputs
write_csv(preddf, '/share/hospital/scratch/marketscan/mr_brt_df_split/preddf.csv')
