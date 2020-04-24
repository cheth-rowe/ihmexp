#############################################################
##Author: Nick Roberts
##Date: 3/4/2019
##Purpose:Quickly check the percentage of points for CF3 where denominator is > 10
##Notes: For validation of concern that these could be causing high uncertainty (value of dropping vs. the importance of their presence for uncertainty)
##Updates: 
#
###########################################################

df <- fread( '/share/hospital/clinical_runs/run_2/estimates/claims/taiwan_and_ms_prep.csv')
df <- df[inp_pri_claims_cases != 0 & inp_otp_any_adjusted_otp_only_indv_cases != 0] # Subset to represent what's actually used in the model

df[, num_bundle_rows := .N, by = 'bundle_id']

df_10 <- df[inp_pri_claims_cases >= 10]
df_10[, denom_10 := .N, by = 'bundle_id']
df_5 <- df[inp_pri_claims_cases >= 5]
df_5[, denom_5 := .N, by = 'bundle_id']

bundle_ten <- merge(unique(df_10[, c('bundle_id', 'denom_10')]), unique(df[, c('bundle_id', 'num_bundle_rows')]), by = 'bundle_id', all.x = TRUE, all.y = TRUE) %>%
  merge(unique(df_5[, c('bundle_id', 'denom_5')]), by = 'bundle_id', all.x = TRUE, all.y = TRUE)
bundle_ten[is.na(denom_10), denom_10 := 0][is.na(denom_5), denom_5 := 0]

bundle_ten[, prop_dropped_10 := 1-denom_10/num_bundle_rows]
bundle_ten[, prop_dropped_5 := 1-denom_5/num_bundle_rows]

bundle_ten <- merge(bundle_ten, bundle_names[, c('bundle_id', 'bundle_name')], by = 'bundle_id', all.x = TRUE)
setorder(bundle_ten, -denom_10)
bundle_ten <- bundle_ten[, c('bundle_id', 'bundle_name', 'denom_5', 'denom_10', 'num_bundle_rows', 'prop_dropped_5', 'prop_dropped_10')] %>% setorder(prop_dropped_10)

write_csv(bundle_ten, '/home/j/temp/nrober75/bundles_inp_admit_counts.csv')
