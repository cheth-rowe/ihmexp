#############################################################
##Author: Nick Roberts
##Date: 3/4/2019
##Purpose:Prep MR-BRT for ICGs
##Notes: 
##Updates:same as prep_all_cfs.R but at ICG level. TODO: When you work this in the pipeline, make a command to do ICGs or bundles. Very simple to do both.
#
###########################################################
run <- commandArgs()[5]
user <- Sys.info()[7]

library(data.table)
library(readr)
source(paste0('/share/code/hospital/', user, '/clinical_info/clinical_info/Functions/db_utilities.R'))
source("/ihme/code/st_gpr/central/r_functions/utilities/utility.r") # For envelope
source('/ihme/cc_resources/libraries/current/r/get_covariate_estimates.R')
source('/ihme/cc_resources/libraries/current/r/get_age_metadata.R')
source('/ihme/cc_resources/libraries/current/r/get_population.R')
ages <- get_age_metadata(12, gbd_round_id = 5) %>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
age_1 <- data.table(age_group_id = 28, age_group_years_start = 0, age_group_years_end = 1)%>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
ages <- ages[age_group_id > 4] %>% rbind(age_1, fill = TRUE)

bundle_map <- bundle_icg()

# Download and append files
hosp_files <- Sys.glob('/share/hospital/clinical_runs/run_2/estimates/claims/02a_non_ms_mapped/icg/*.csv')
hosp_dt <- rbindlist(lapply(hosp_files, fread))
hosp_dt[is.na(hosp_dt)] <- 0
hosp_dt <- age_binner(hosp_dt)

# Aggregate over ICG just in case they aren't
hosp_dt <- hosp_dt[, lapply(.SD, sum), by = c('location_id', 'year_start', 'age_start', 'icg_id', 'sex_id'), 
                   .SDcols = c('inp_pri_claims_cases','inp_pri_indv_cases','inp_any_claims_cases','inp_any_indv_cases')]
# Aggregate over bundle (annoying duplication because of mapping)
# hosp_dt <- merge(hosp_dt, bundle_map, by = 'icg_id', allow.cartesian = TRUE)
# hosp_dt <- hosp_dt[, lapply(.SD, sum), by = c('location_id', 'year_start', 'age_start', 'bundle_id', 'sex_id'), 
#                    .SDcols = c('inp_pri_claims_cases', 'inp_pri_indv_cases','inp_any_claims_cases','inp_any_indv_cases')]
hosp_dt[, cf1 := inp_pri_indv_cases/inp_pri_claims_cases]
hosp_dt[, cf2 := inp_any_indv_cases/inp_pri_claims_cases]

# Download and merge on haq
haqi <- get_covariate_estimates(1099, location_id = unique(hosp_dt$location_id), sex_id = 3,
                                year_id = unique(hosp_dt$year_start), gbd_round_id = 6, decomp_step = 'step1') %>%
  setnames(c('year_id', 'mean_value'), c('year_start', 'haqi'))
hosp_dt <- merge(hosp_dt, haqi[, c('location_id', 'year_start', 'haqi')], by = c('location_id', 'year_start'))
hosp_dt <- merge(hosp_dt, ages[, c('age_start', 'age_group_id')], by = 'age_start')

# Need to download the envelope and population to get an estimate of inpatient sample size
pops <- get_population(year_id = unique(hosp_dt$year_start), sex_id = c(1,2), location_id = unique(hosp_dt$location_id), age_group_id = unique(hosp_dt$age_group_id),
                       gbd_round_id = 6, decomp_step = 'step1') %>% setnames('year_id', 'year_start')
ip_envelope <- model_load(47570, 'gpr') %>% setnames('year_id', 'year_start') %>% merge(ages[, c('age_group_id', 'age_start')], by = 'age_group_id')

samp_size <- merge(pops[, c('location_id', 'year_start', 'age_group_id', 'sex_id', 'population')], 
                   ip_envelope[, c('age_group_id', 'year_start','sex_id', 'location_id', 'gpr_mean')], by = c('sex_id', 'age_group_id', 'year_start', 'location_id'))
samp_size[, sample_size := gpr_mean*population]
hosp_dt <- merge(hosp_dt, samp_size[, c('sex_id', 'age_group_id', 'year_start', 'location_id', 'sample_size')], by = c('sex_id', 'age_group_id', 'year_start', 'location_id'))

# Calculate sample size
# std_error = sqrt(1 / sample_size * cf * (1 - cf) + 1 / (4 * sample_size^2) * 1.96^2)
# Just need to calculate the total claims to get a cause fraction
hosp_dt[, total_inp_pri_claims_cases := sum(inp_pri_claims_cases), by = c('age_start', 'location_id', 'sex_id', 'year_start')]
hosp_dt[, cf := inp_pri_claims_cases/total_inp_pri_claims_cases]
hosp_dt[, se := sqrt(1/sample_size*cf*(1-cf) + 1/(4*sample_size^2)*1.96^2)] # Should check these order of operations, but this is straigth from Stata

# Drop infinite values and NAs, split by bundle and CF
mod_cf1 <- hosp_dt[!is.na(cf1) & cf1 != Inf]
mod_cf2 <- hosp_dt[!is.na(cf2) & cf2 != Inf]

##### Prep Marketscan inputs ######
# Download and format data
df <- fread('/share/hospital/clinical_runs/run_2/estimates/claims/02_ms_mapped/condensed/icg_condensed.csv')
df <- dcast(df, age_end + age_start + icg_id + location_id + sex_id + year_end + year_start ~ estimate_type, value.var = 'val')
df[is.na(df)] <- 0
df[, cf3 := inp_otp_any_adjusted_otp_only_indv_cases/inp_pri_claims_cases]
df[, cf2 := inp_any_indv_cases/inp_pri_claims_cases]
df[, cf1 := inp_pri_indv_cases/inp_pri_claims_cases]

# Drop US national data, not representative like the states
df <- df[location_id != 102]

# Calculate variance on each point, using wilson's method 
# std_error = sqrt(1 / sample_size * cf * (1 - cf) + 1 / (4 * sample_size^2) * 1.96^2)
# Just need to calculate the total claims to get a cause fraction
df[, total_inp_pri_claims_cases := sum(inp_pri_claims_cases), by = c('age_start', 'location_id', 'sex_id', 'year_start')]
df[, cf := inp_pri_claims_cases/total_inp_pri_claims_cases]

# Get sample size, need Marketscan sample size files and an egeoloc map to get to location id
egeoloc_map <- fread('/home/j/temp/hospital/2016/data/maps/ms_egeoloc_map.csv') # For mergeing
enrolee_files <- Sys.glob('/share/hospital/scratch/marketscan/sample_size/*.dta')
enrolees <- rbindlist(lapply(enrolee_files, read.dta13)) %>% 
  merge(egeoloc_map, by = 'egeoloc') %>% 
  setnames('age_start', 'age') %>% age_binner()
enrolees[, sample_size := sum(sample_size), by = c('age_start', 'sex', 'year', 'location_id')] # Aggregate over 5-year age groups
enrolees <- unique(enrolees[, c('age_start', 'sex', 'year', 'location_id', 'sample_size', 'state')]) %>% setnames(c('sex', 'year'),                                                                                                                  c('sex_id', 'year_start'))
enrolees[, sex_id := as.numeric(sex_id)]

df <- merge(df, enrolees, by = c('age_start', 'sex_id', 'year_start', 'location_id'))
df[, se := sqrt(1/sample_size*cf*(1-cf) + 1/(4*sample_size^2)*1.96^2)] # Should check these order of operations, but this is straigth from Stata

# Merge
haqi <- get_covariate_estimates(1099, gbd_round_id = 6, decomp_step = 'step1') %>% setnames(c('year_id', 'mean_value'), c('year_start', 'haqi'))
df <- merge(df, haqi[, c('location_id', 'year_start', 'haqi')], by = c('location_id', 'year_start'))
df <- merge(df, ages[, c('age_start', 'age_group_id')], by = 'age_start')
#Should we drop Puerto Rico? Maybe don't and see what happens
cf1_cols <- names(mod_cf1)
cf2_cols <- names(mod_cf2)
ms_cf1 <- df[, ..cf1_cols] %>% .[!is.na(cf1) & cf1 != Inf & cf1 != 0]
ms_cf2 <- df[, ..cf2_cols] %>% .[!is.na(cf2) & cf2 != Inf & cf2 != 0]

# Append hospital and marketscna
mod_cf1 <- rbind(mod_cf1, ms_cf1)
mod_cf2 <- rbind(mod_cf2, ms_cf2)

# Weird sample sizes less than 1 causing SE's to be high. SS <1 doens't make sense. Just due to tiny populations
mod_cf1[sample_size < 1, sample_size := 1]
mod_cf1[, logit_se := logit(se)]
mod_cf1 <- mod_cf1[!is.na(logit_se)]

# Get age midpoint
mod_cf1 <- merge(mod_cf1, ages[, c('age_group_id', 'age_end')])
mod_cf2 <- merge(mod_cf2, ages[, c('age_group_id', 'age_end')])
mod_cf1[, age_start := (age_start + age_end)/2]
mod_cf2[, age_start := (age_start + age_end)/2]

# Split
icgs <- unique(mod_cf1$icg_id)
invisible(lapply(icgs, function(x){
  print(x)
  subset <- mod_cf1[icg_id == x & cf1 != 0]
  subset[cf1 == 1, cf1 := 0.999]
  subset[, logit_mean := logit(cf1)]
  subset[, logit_se := logit(se)]
  subset[, log_haqi := log(haqi)]
  
  subset[, intercept := 1]
  write_csv(subset, paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/prep_icg/cf1/', x, '.csv'))
}))

icgs <- unique(mod_cf2$icg_id)
invisible(lapply(icgs, function(x){
  print(x)
  subset <- mod_cf2[icg_id == x & cf2 != 0]
  subset[, log_mean := log(cf2)]
  subset[, log_se := log(se)]
  subset[, log_haqi := log(haqi)]
  
  subset[, intercept := 1]
  write_csv(subset, paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/prep_icg/cf2/', x, '.csv'))
}))

icgs <- unique(df$icg_id)
invisible(lapply(icgs, function(x){
  print(x)
  subset <- df[icg_id == x & cf3 != 0 & !is.na(cf3) & cf3 != Inf]
  subset[, log_mean := log(cf3)]
  subset[, log_se := log(se)]
  subset[, log_haqi := log(haqi)]
  
  subset[, intercept := 1]
  write_csv(subset, paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/prep_icg/cf3/', x, '.csv'))
}))


