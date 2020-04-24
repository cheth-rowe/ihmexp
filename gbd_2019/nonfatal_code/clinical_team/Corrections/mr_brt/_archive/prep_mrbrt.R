#############################################################
##Author: Nick Roberts
##Date: 3/4/2019
##Purpose: MR-BRT CF3 application
##Notes: First go. Need to make sure we hvae marketscan data by year and state (or just state)
##Updates: Preps and splits it by bundle
#
###########################################################
user <- Sys.info()[7][[1]]
library(data.table)
library(parallel)
library(readstata13)
library(fastDummies, lib.loc = '/home/j/temp/nrober75/packages')
source(paste0('/share/code/hospital/', user, '/clinical_info/clinical_info/Functions/db_utilities.R'))
source('/ihme/cc_resources/libraries/current/r/get_age_metadata.R')
source('/ihme/cc_resources/libraries/current/r/get_covariate_estimates.R')
source('/ihme/cc_resources/libraries/current/r/get_population.R')
source('/ihme/cc_resources/libraries/current/r/get_model_results.R')

ages <- get_age_metadata(12, gbd_round_id = 5) %>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
age_1 <- data.table(age_group_id = 28, age_group_years_start = 0, age_group_years_end = 1)%>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
ages <- ages[age_group_id > 4] %>% rbind(age_1, fill = TRUE) %>% setnames(c('age_group_years_start', 'age_group_years_end'),
                                                                          c('age_start', 'age_end'))

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

### Get outpatient envelope ####
op_env <- get_model_results(gbd_team = 'epi', 
                            gbd_round_id = 5, location_set_id = 35,
                            model_version_id = 246617,
                            age_group_id = c(2:20, 28, 30, 31, 32, 33, 164, 235))

op_env <- op_env[, c('location_id', 'year_id', 'age_group_id', 'sex_id', 'mean')] %>% setnames('mean', 'envelope_mean')
op_env <- interp_under1_env(op_env)
setnames(op_env, 'year_id', 'year_start')
op_env <- op_env[location_id %in% unique(df$location_id)]



ages <- copy(ages_over1)
bundle_map <- bundle_icg()       

# Download and format data
df <- fread('/share/hospital/clinical_runs/run_1/estimates/claims/02_ms_mapped/condensed/bundle_condensed.csv')
df <- dcast(df, age_end + age_start + bundle_id + location_id + sex_id + year_end + year_start ~ estimate_type, value.var = 'val')
df[is.na(df)] <- 0
df[, cf3 := inp_otp_any_adjusted_otp_only_indv_cases/inp_pri_claims_cases]
df[inp_pri_claims_cases == 0, inp_pri_claims_cases_plus1 := 1][inp_pri_claims_cases > 0, inp_pri_claims_cases_plus1 := inp_pri_claims_cases]
df[inp_otp_any_adjusted_otp_only_indv_cases != 0, cf3_no_na := inp_otp_any_adjusted_otp_only_indv_cases/inp_pri_claims_cases_plus1] # Removes infinites, but leaves zero's as NAs to not impute the zero

# Drop US national data, not representative like the states
df <- df[location_id != 102]

# Calculate variance on each point, using wilson's method 
# std_error = sqrt(1 / sample_size * cf * (1 - cf) + 1 / (4 * sample_size^2) * 1.96^2)
# Just need to calculate the total claims to get a cause fraction
df[, total_inp_pri_claims_cases := sum(inp_pri_claims_cases), by = c('age_start', 'location_id', 'sex_id', 'year_start')]
df[, total_inp_pri_claims_cases_plus := sum(inp_pri_claims_cases_plus1), by = c('age_start', 'location_id', 'sex_id', 'year_start')]
df[, cf := inp_pri_claims_cases/total_inp_pri_claims_cases]
df[, cf_plus := inp_pri_claims_cases_plus1/total_inp_pri_claims_cases_plus]

# Get sample size, need Marketscan sample size files and an egeoloc map to get to location id
egeoloc_map <- fread('/home/j/temp/hospital/2016/data/maps/ms_egeoloc_map.csv') # For mergeing
enrolee_files <- Sys.glob('/share/hospital/scratch/marketscan/sample_size/*.dta')
enrolees <- rbindlist(lapply(enrolee_files, read.dta13)) %>% 
  merge(egeoloc_map, by = 'egeoloc') %>% 
  setnames('age_start', 'age') %>% age_binner()
enrolees[, sample_size := sum(sample_size), by = c('age_start', 'sex', 'year', 'location_id')] # Aggregate over 5-year age groups
enrolees <- unique(enrolees[, c('age_start', 'sex', 'year', 'location_id', 'sample_size', 'state')]) %>% setnames(c('sex', 'year'),
                                                                                                                  c('sex_id', 'year_start'))
enrolees[, sex_id := as.numeric(sex_id)]

# Merge and calculate variance
df <- merge(df, enrolees, by = c('age_start', 'sex_id', 'year_start', 'location_id'))
df[, se := sqrt(1/sample_size*cf*(1-cf) + 1/(4*sample_size^2)*1.96^2)] # Should check these order of operations, but this is straigth from Stata
df[, se_plus := sqrt(1/sample_size*cf_plus*(1-cf_plus) + 1/(4*sample_size^2)*1.96^2)] # Should check these order of operations, but this is straigth from Stata

# Merge on covariates of interest (ip_envelope, haqi, sdi, op_envelope)
source("/ihme/code/st_gpr/central/r_functions/utilities/utility.r")

ip_envelope <- model_load(47570, 'gpr') %>% setnames('year_id', 'year_start') %>% merge(ages[, c('age_group_id', 'age_start')], by = 'age_group_id')
df <- merge(df, ip_envelope[, c('location_id', 'year_start', 'age_start', 'sex_id', 'gpr_mean')])
setnames(df, 'gpr_mean', 'ip_env')

# HAQ and SDI
haqi <- get_covariate_estimates(1099, gbd_round_id = 6, decomp_step = 'step1') %>% setnames(c('year_id', 'mean_value'), c('year_start', 'haqi'))
sdi <- get_covariate_estimates(881, gbd_round_id = 6, decomp_step = 'step1') %>% setnames(c('year_id', 'mean_value'), c('year_start', 'sdi'))
df <- merge(df, haqi[, c('location_id', 'year_start', 'haqi')], by = c('location_id', 'year_start'))
df <- merge(df, sdi[, c('location_id', 'year_start', 'sdi')], by = c('location_id', 'year_start'))

## Format outpatient envelope. Need to interpolate to get results for the individual years
source("/ihme/cc_resources/libraries/current/r/interpolate.R")
# env_interp <- interpolate('modelable_entity_id', gbd_id = 19797, version_id = 246617, source = 'epi', measure_id = 19, sex_id = c(1,2), 
#                           gbd_round_id= 5, reporting_year_start = 2010, 
#                           reporting_year_end = 2017)
# env_interp <- env_interp[year_id <= 2015][location_id %in% df$location_id] %>% setnames('year_id', 'year_start')
env_interp <- fread('/home/j/temp/nrober75/USA_interpolated_OP_env.csv') %>% setnames('year_id', 'year_start')

env_interp_long <- melt(env_interp, id.vars = c('age_group_id', 'location_id', 'year_start', 'sex_id'), 
                        measure.vars = patterns('draw_')) 
env_interp_long[, envelope_mean := mean(value), by = c('age_group_id', 'location_id', 'year_start', 'sex_id')] # Get mean of the draws
env_interp_long[, value := NULL][, variable := NULL]
setnames(env_interp_long, 'year_start', 'year_id')
env_interp_long <- env_interp_long[location_id %in% unique(df$location_id)]
env_interp_long <- interp_under1_env(env_interp_long)
env_interp_long <- unique(env_interp_long) %>% setnames('year_id', 'year_start')
env_interp_long <- env_interp_long[age_group_id != 164]

# Merge on 2000 op envelope
op_env <- merge(op_env, ages[, c('age_group_id', 'age_start', 'age_end')])
op_env <- op_env[location_id %in% df$location_id & year_start == 2000]
op_env <- rbind(op_env, env_interp_long)
df <- merge(df, op_env[, c('age_start', 'location_id', 'year_start', 'sex_id', 'envelope_mean')], 
            by = c('age_start', 'location_id', 'year_start', 'sex_id'), all.x = TRUE)
df <- df[location_id > 500] %>% setnames('envelope_mean', 'op_env')
write_csv(df, '/share/hospital/scratch/marketscan/mr_brt_df_split/marketscan_prepped.csv')

##### Prep for MR-BRT analysis by splitting the main dataframe by each bundle #####
# Split into their own csv in /share/hospital/scratch/marketscan/mr_brt_df_split/bundle
bundles <- unique(df$bundle_id)
lapply(bundles, function(x){
  print(x)
  subset <- df[bundle_id == x][!is.na(cf3)]
  subset[, log_mean := log(cf3)]
  subset[, log_se := log(se)]
  subset <- subset[, c('log_mean', 'log_se', 'age_start', 'sex_id', 'haqi', 'ip_env', 'op_env')]
  subset[, intercept := 1]
  write_csv(subset, paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/with_inf/', x, '.csv'))
})

# With dropped infinites
lapply(bundles, function(x){
  print(x)
  subset <- df[bundle_id == x][!is.na(cf3) & cf3 != Inf]
  subset[, log_mean := log(cf3)]
  subset[, log_se := log(se)]
  subset <- subset[, c('log_mean', 'log_se', 'age_start', 'sex_id', 'haqi', 'ip_env', 'op_env')]
  subset[, intercept := 1]
  write_csv(subset, paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', x, '.csv'))
})


# #### Taiwan prep ####
# ## Prep Taiwan too ( From CF3 script)
# taiwan <- fread("/home/j/temp/hospital/2016/data/new sources/taiwan/taiwan_1y.csv") ## Prep uses this one
# taiwan <- taiwan[, c('bundle_id', 'sex', 'age_ihmec', 'correction1', 'correction2', 'correction3')]
# setnames(taiwan, c('correction1', 'correction2', 'correction3'), c('indv_cf', 'incidence', 'prevalence'))
# taiwan <- as.data.table(taiwan)
# colnames(taiwan)[3] <- 'age_end'
# taiwan$age_end <- as.numeric(taiwan$age_end)
# colnames(taiwan)[2] <- 'sex_id'
# taiwan <- taiwan[age_end == 0, age_end := 1][age_end == 95, age_end := 99]
# taiwan$bundle_id <- as.character(taiwan$bundle_id)
# taiwan <- taiwan[bundle_id > 0]
# taiwan$bundle_id <- as.numeric(taiwan$bundle_id)
# 
# # And format a little more
# taiwan[, location_id := 8][, year_start := 2016][age_end != 1 & age_end != 99, age_end := age_end + 1][age_end == 99, age_end := 125]
# taiwan <- merge(taiwan, ages[, c('age_group_id', 'age_start', 'age_end')], by = 'age_end')
# setnames(taiwan, 'prevalence', 'cf3')
# 
# # Use populatoin to get sample size. 99.9% enrolled
# taiwan_pop <- get_population(age_group_id = ages$age_group_id, year_id = 2000:2017, sex_id = c(1,2), gbd_round_id = 6, decomp_step = 'step1', location_id = 8) %>%
#   setnames(c('year_id', 'population'), c('year_start', 'sample_size'))

