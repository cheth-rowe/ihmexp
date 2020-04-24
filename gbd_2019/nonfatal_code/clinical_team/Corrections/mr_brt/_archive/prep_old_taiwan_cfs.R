#############################################################
##Author: Nick Roberts
##Date: 4/3/2019
##Purpose:Compare maps and prep Taiwan in bundles that didn't change across maps
##Notes: 
##Updates: Takes some dB functions to compare. Then move the saved file into prep_all_cfs.R. Did this to keep prep_all_cfs clean. 
# WILL NEED TO RUN THIS WITH NEW MAPS TO CHECK FOR MAP CHANGES (UNTIL TWN DATA COMES IN )
###########################################################
library(readstata13)
#### Compare maps to add Taiwan ####
icd_icg_version <- function(version){
  db_con = fread(paste0(j, "/temp/hospital/2016/code/db_con.csv"))
  
  # Get list of ids and names
  myconn <- RMySQL::dbConnect(RMySQL::MySQL(),
                              host = odbc[[con_def]]$SERVER,
                              username = odbc[[con_def]]$USER,
                              password = odbc[[con_def]]$PASSWORD)
  df <- dbGetQuery(myconn, sprintf("
                                   SELECT
                                   *
                                   FROM
                                   clinical_mapping.cause_code_icg"))
  df$cause_code <- gsub('[.]', '', df$cause_code)
  # max_version <- max(df$map_version)
  df <- as.data.table(df)
  df <- df[map_version == version]
  #WHERE
  #bundle_id IN (%s),
  #paste0(bundles, collapse = ",")))
  dbDisconnect(myconn)
  return(data.table(df))
}

bundle_icg_version <- function(version){
  db_con = fread(paste0(j, "/temp/hospital/2016/code/db_con.csv"))
  
  # Get list of ids and names
  myconn <- RMySQL::dbConnect(RMySQL::MySQL(),
                              host = odbc[[con_def]]$SERVER,
                              username = odbc[[con_def]]$USER,
                              password = odbc[[con_def]]$PASSWORD)
  df <- dbGetQuery(myconn, sprintf("
                                   SELECT
                                   *
                                   FROM
                                   clinical_mapping.icg_bundle"))

  # max_version <- max(df$map_version)
  df <- as.data.table(df)
  df <- df[map_version == version]
  #WHERE
  #bundle_id IN (%s),
  #paste0(bundles, collapse = ",")))
  dbDisconnect(myconn)
  return(data.table(df))
}

new_map <- icd_icg_version(23)
new_bundle_map <- bundle_icg_version(23)
new_icd_map <- merge(new_map, new_bundle_map, by = c('icg_name', 'icg_id'), allow.cartesian = TRUE) %>%
  .[, c('cause_code', 'code_system_id', 'bundle_id')] %>% 
  setnames('bundle_id', 'new_bundle_id')

old_map <- icd_icg_version(20) 
old_bundle_map <- bundle_icg_version(20) %>% setnames('bundle_id', 'old_bundle_id')
old_icd_map <- merge(old_map, old_bundle_map, by = c('icg_name', 'icg_id'), allow.cartesian = TRUE) %>% 
  .[, c('cause_code', 'code_system_id', 'old_bundle_id')]

# Do in a quick for loop
# checked this, every bundle in the old map is present in the new one. 
# So go through all bundles in old map. cause nothing extra in new map will be used and they're all present
#sink('/home/j/temp/nrober75/bundle_comps.txt')
exact_bundles <- c()
for(bundle in unique(old_icd_map$old_bundle_id)){
  new_subset <- new_icd_map[new_bundle_id == bundle & code_system_id %in% c(1,2)]
  new_subset[, num_icd := NULL]
  old_subset <- old_icd_map[old_bundle_id == bundle & code_system_id %in% c(1,2)]
  old_subset[, num_icd := NULL]
  
  if(all.equal(new_subset$cause_code, old_subset$cause_code) == TRUE){
    print('SAME!')
    exact_bundles <- c(exact_bundles, bundle)
  }
  
  #print(bundle)
  # print(setdiff(new_subset$cause_code, old_subset$cause_code))
  # print(setdiff(old_subset$cause_code, new_subset$cause_code))
}
exact_bundles <- data.table(bundle_id = exact_bundles)
#sink()
write_csv(exact_bundles, '/home/j/temp/nrober75/exact_bundles.csv')

ages <- get_age_metadata(12, gbd_round_id = 5) %>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
age_1 <- data.table(age_group_id = 28, age_group_years_start = 0, age_group_years_end = 1)%>% 
  setnames(c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
ages <- ages[age_group_id >= 5] %>% rbind(age_1, fill = TRUE)


taiwan <- fread("/home/j/temp/hospital/2016/data/new sources/taiwan/taiwan_1y.csv")
taiwan <- taiwan[bundle_id != '.']
taiwan <- taiwan[correction3 < 1000000] %>%
  .[, c('bundle_id', 'sex', 'age_ihmec', 'correction1', 'correction2', 'correction3')] %>%
  setnames(c('age_ihmec', 'correction1', 'correction2', 'correction3', 'sex'), c('age_end', 'cf1', 'cf2', 'cf3', 'sex_id'))
taiwan[, age_end := age_end + 1]
taiwan[age_end == 96, age_end := 125]
taiwan <- merge(taiwan, ages[, c('age_start', 'age_end')], by = 'age_end')
taiwan[, age_start := (age_start + age_end)/2]
taiwan <- taiwan[bundle_id %in% exact_bundles]
taiwan[age_start == 110, age_start := 97.5]
## Need to make standard error calculations and get sample size

twn_raw <- read.dta13('/home/j/DATA/TWN/NATIONAL_HEALTH_INSURANCE/CLAIMS/TWN_NATIONAL_HEALTH_INSURANCE_CLAIMS_2016_Y2018M02D14.DTA')
setnames(twn_raw, 'diagnosis_1', 'cause_code')
twn_mapped <- merge(twn_raw, old_map[, c('cause_code', 'icg_name')], by = 'cause_code', all.x = TRUE)
twn_mapped <- as.data.table(twn_mapped)
twn_mapped <- twn_mapped[outcome == 'Discharge']
setnames(twn_mapped, 'age_value', 'age')
twn_mapped[age_unit == 'days', age := 0]
twn_mapped <- age_binner(twn_mapped)
twn_mapped <- twn_mapped[!is.na(age_start)]

# Count numerator and denominator, by admission and not person
twn_mapped[, total_admit := .N, by = c('year', 'sex', 'age_start')]
twn_mapped[, total_admit_icg := .N, by = c('year', 'sex', 'age_start', 'icg_name')]

twn_mapped <- twn_mapped[, c('year', 'sex', 'age_start', 'age_end', 'icg_name', 'total_admit_icg', 'total_admit')]
twn_mapped <- merge(twn_mapped, old_bundle_map[, c("icg_name", "old_bundle_id")], by = 'icg_name', allow.cartesian = TRUE, all.x = TRUE) %>% setnames("old_bundle_id", "bundle_id")
twn_mapped <- unique(twn_mapped[, c('icg_name', 'year', 'sex', 'age_start', 'age_end', 'bundle_id', 'total_admit_icg', 'total_admit')])
twn_mapped[, sum_bundle_admit := sum(total_admit_icg), by = c('year', 'sex', 'age_start', 'bundle_id')]

# Calculate standard error with Wilsons. Using the number of admissions for sample size for now (can't tell which col is the patient identifier or if we have that...)
twn_mapped[, cf := sum_bundle_admit/total_admit]
twn_mapped[, sample_size := total_admit]
twn_mapped[, se := sqrt(1/sample_size*cf*(1-cf) + 1/(4*sample_size^2)*1.96^2)] # Should check these order of operations, but this is straigth from Stata
twn_mapped <- merge(twn_mapped, ages[, c('age_start', 'age_group_id')], by = 'age_start')
twn_mapped[, age_end := NULL]
twn_mapped <- merge(twn_mapped, ages[, c('age_group_id', 'age_end')], by = 'age_group_id')

twn_mapped[, age_start := (age_start + age_end)/2]
twn_mapped[age_start == 110, age_start := 97.5]

# do only 2015 and 2016, no other years
twn_mapped <- twn_mapped[year == 2016]
setnames(twn_mapped, c('sex', 'year'), c('sex_id', 'year_id'))
twn_mapped[, sex_id := as.character(sex_id)][, bundle_id := as.character(bundle_id)]
taiwan[, sex_id := as.character(sex_id)][, bundle_id := as.character(bundle_id)]

taiwan <-  merge(taiwan, twn_mapped[, c('age_group_id', 'age_start', 'sex_id', 'bundle_id', 'cf', 'sample_size', 'se')], by = c('bundle_id', 'age_start', 'sex_id'),all.x = TRUE )
taiwan <- taiwan[!is.na(se)]
taiwan <- taiwan[!is.na(bundle_id)]
taiwan[, year_start := 2016][, location_id := 8]
taiwan <- unique(taiwan[, c('age_group_id', 'age_start', 'year_start', 'sex_id', 'bundle_id', 'cf', 'sample_size', 'se', 'age_end', 'location_id',
                            'cf1', 'cf2', 'cf3')])
taiwan <- taiwan[bundle_id %in% exact_bundles]
write_csv(taiwan, '/share/hospital/clinical_runs/run_2/estimates/claims/taiwan_cf_prep.csv')
