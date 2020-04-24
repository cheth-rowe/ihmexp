#############################################################
##Author: Nick Roberts
##Date: 3/8/2019
##Purpose: Scatter the MR-BRT and LMER model predictions for each bundle
##Notes: MR-BRT models are all done and predicted. LMER only one that needs to be re-done
##Updates:
#
###########################################################


library(data.table)
library(tibble)
library(lme4)
library(RMySQL)
library(dplyr)

mround <- function(x,base){ 
  base*round(x/base) 
} 
# Main data for all bundles (that modles were run on)
all_bundles <- fread('/share/hospital/scratch/marketscan/mr_brt_df_split/marketscan_prepped.csv')

# MR-BRT results (needs a bundle ID argument)
results_files <- Sys.glob('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/results/*/model_summaries.csv')


## Get DF and make model
df <- all_bundles[bundle_id == bundle][!is.na(cf3) & cf3 != Inf & cf3 != 0]
df[, log_mean := log(cf3)]
df[, log_se := log(se)]
df <- df[bundle_id == bundle, c('log_mean', 'log_se', 'age_start', 'sex_id', 'haqi', 'ip_env', 'op_env', 'location_id', 'year_start')]

mod <- lmer(log_mean ~ age_start + sex_id + haqi + ip_env + op_env + (1|location_id), data = df)

# Predict on results_files covariates
brt_result <- fread(paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/results/',bundle,'/model_summaries.csv'))

# Round state-level stuff to appropriate level for predictions
pred_df <- copy(df)
pred_df[, op_env := round(op_env, 0)]
pred_df[, ip_env := mround(ip_env, .05)]
pred_df[, haqi := round(haqi, 0)]
pred_df[, pred := predict(mod, newdata = pred_df)]

pred_df <- pred_df %>% 
  mutate_all(as.character)
brt_result <- brt_result %>% mutate_all(as.character)

pred_df_merge <- merge(pred_df, brt_result[, c('Y_mean', 'op_env', 'ip_env', 'haqi', 'sex_id', 'age_start')], 
                 by = c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start'), all.x = TRUE)
pred_df_merge <- as.data.table(pred_df_merge)

pred_df_merge[, pred := as.numeric(pred)]
pred_df_merge[, Y_mean := as.numeric(Y_mean)]

## Now plot
ggplot(pred_df_merge, aes(x = exp(pred), y = exp(Y_mean))) + 
  geom_point() + 
  xlab('LMER Model') + ylab('MR-BRT model') +
  ggtitle('MR-BRT vs. LMER scatter on USA states') + 
  theme_bw() +
  geom_abline(slope = 1, intercept = 0)

# And check out fixed-effects too
pred_df_no_loc <- copy(pred_df) %>% as.data.table()

pred_df_no_loc[, location_id := 1]

mod_lm <- glm(log_mean ~ age_start + sex_id + haqi + ip_env + op_env, data = df)
pred_df_no_loc[, pred_lm := ]
