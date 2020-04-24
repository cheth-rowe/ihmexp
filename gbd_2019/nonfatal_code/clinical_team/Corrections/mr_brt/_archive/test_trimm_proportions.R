#############################################################
##Author: Nick Roberts
##Date: 3/8/2019
##Purpose: Run with Reed's new modeling stuff
##Notes: Run this, should give you an easier way to make predictions on out-of-sample data
##Updates: This is a garbage script. Used just to test things to go into a more real script to run. Don't expect anything from me
#
###########################################################

repo_dir <- "/home/j/temp/reed/prog/projects/run_mr_brt/"
source(paste0(repo_dir, "run_mr_brt_function.R"))
source(paste0(repo_dir, "cov_info_function.R"))
source(paste0(repo_dir, "check_for_outputs_function.R"))
source(paste0(repo_dir, "load_mr_brt_outputs_function.R"))
source(paste0(repo_dir, "predict_mr_brt_function.R"))
source(paste0(repo_dir, "check_for_preds_function.R"))
source(paste0(repo_dir, "load_mr_brt_preds_function.R"))


bundle <- 19
fit1 <- run_mr_brt(
  output_dir = paste0("/share/hospital/scratch/marketscan/mr_brt_df_split/models/run_2/trim_10"), 
  model_label = paste0(bundle, '_run'),
  data = paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', bundle, '.csv') ,
  mean_var = "log_mean",
  se_var = "log_se",
  covs = list( 
    #cov_info("zcov1", "Z"),
    cov_info("op_env", "X"),
    cov_info("ip_env", "X"),
    cov_info("haqi", "X"),
    cov_info("sex_id", "X"),
    cov_info("age_start", "X") ),
  overwrite_previous = TRUE,
  project = 'proj_hospital',
  method = 'trim_maxL', trim_pct = .1
)
#check_for_outputs(fit1, wait_seconds = 10)
results1 <- load_mr_brt_outputs(fit1)

# Make predictions
all_bundles <- fread('/share/hospital/scratch/marketscan/mr_brt_df_split/marketscan_prepped.csv')
df <- all_bundles[bundle_id == 19]
df1 <- fread(paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', bundle, '.csv'))

pred1 <- predict_mr_brt(fit1, newdata = df1)
check_for_preds(pred1)

pred_object <- load_mr_brt_preds(pred1)
trim_20 <- pred_object$model_summaries %>% setnames('Y_mean', 'trim_20_Y_mean')
trim_20 <- as.data.table(trim_20)

pred_trim_10 <- predict_mr_brt(fit1, newdata = df1) # Different fit1 than above. My mistake.
pred_object_10 <- load_mr_brt_preds(pred_trim_10)
trim_10 <- pred_object_10$model_summaries %>% setnames('Y_mean', 'trim_10_Y_mean') %>% as.data.table()
trim_10 <- as.data.table(trim_10)
# Scatter 10 and 20 trim pcts

plot_trim <- merge(unique(trim_10[, c('X_age_start', 'X_sex_id', 'X_op_env', 'X_ip_env', 'X_haqi', 'trim_10_Y_mean')]),
                   unique(trim_20[, c('X_age_start', 'X_sex_id', 'X_op_env', 'X_ip_env', 'X_haqi', 'trim_20_Y_mean')]),
                   by = c('X_age_start', 'X_sex_id', 'X_op_env', 'X_ip_env', 'X_haqi'))

ggplot(plot_trim, aes(x = exp(trim_20_Y_mean), y = exp(trim_10_Y_mean))) + 
  geom_point(color = 'dark green') +
  geom_abline(slope = 1, intercept = 0)

ten <- ggplot(plot_trim, aes(x = X_age_start, y = exp(trim_10_Y_mean))) + 
                geom_point() +
  ggtitle('10% trim model') +
  ylab('10% Trim CF3') +
  theme_bw() +
  geom_smooth()

twenty <- ggplot(plot_trim, aes(x = X_age_start, y = exp(trim_20_Y_mean))) + 
  geom_point() +
  ggtitle('20% trim model') +
  ylab('20% Trim CF3') +
  theme_bw() +
  geom_smooth()
grid.arrange(ten, twenty)              



# Go crazy with a 40% trim
fit_40 <- run_mr_brt(
  output_dir = paste0("/share/hospital/scratch/marketscan/mr_brt_df_split/models/run_2/trim_40"), 
  model_label = paste0(bundle, '_run'),
  data = paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', bundle, '.csv') ,
  mean_var = "log_mean",
  se_var = "log_se",
  covs = list( 
    #cov_info("zcov1", "Z"),
    cov_info("op_env", "X"),
    cov_info("ip_env", "X"),
    cov_info("haqi", "X"),
    cov_info("sex_id", "X"),
    cov_info("age_start", "X") ),
  overwrite_previous = TRUE,
  project = 'proj_hospital',
  method = 'trim_maxL', trim_pct = .4
)

pred40 <- predict_mr_brt(fit_40, newdata = df1)
check_for_preds(pred40)

pred_object <- load_mr_brt_preds(pred40)
trim_40 <- pred_object$model_summaries %>% setnames('Y_mean', 'trim_40_Y_mean')
trim_40 <- as.data.table(trim_40)

forty <- ggplot(trim_40, aes(x = X_age_start, y = exp(trim_40_Y_mean))) + 
  geom_point() +
  ggtitle('40% trim model') +
  ylab('40% Trim CF3') +
  theme_bw() +
  geom_smooth()

grid.arrange(ten, twenty, forty)              

## Now do a lasso
lasso_20 <- run_mr_brt(
  output_dir = paste0("/share/hospital/scratch/marketscan/mr_brt_df_split/models/run_2/lasso_20"), 
  model_label = paste0(bundle, '_run'),
  data = paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', bundle, '.csv') ,
  mean_var = "log_mean",
  se_var = "log_se",
  covs = list( 
    #cov_info("zcov1", "Z"),
    cov_info("op_env", "X"),
    cov_info("ip_env", "X"),
    cov_info("haqi", "X"),
    cov_info("sex_id", "X"),
    cov_info("age_start", "X") ),
  overwrite_previous = TRUE,
  project = 'proj_hospital',
  method = 'trim_maxL', trim_pct = .2, lasso = TRUE
)

pred_lasso_20 <- predict_mr_brt(lasso_20, newdata = df1)
check_for_preds(pred40)

pred_object <- load_mr_brt_preds(pred_lasso_20)
lasso_20 <- pred_object$model_summaries %>% setnames('Y_mean', 'lasso_20_Y_mean')
lasso_20 <- as.data.table(lasso_20)

lasso_20_plot <- ggplot(lasso_20, aes(x = X_age_start, y = exp(lasso_20_Y_mean))) + 
  geom_point() +
  ggtitle('20% lasso model') +
  ylab('20% lasso model') +
  theme_bw() +
  geom_smooth()

system(lasso_20$sh_command)
shQuote(paste0("source /ihme/code/evidence_score/miniconda3/bin/activate mr_brt_env; python /ihme/code/evidence_score/mr_brt_ihme/uw_amo/mr_brt.py \\", lasso_20$sh_command))

## Try to run using Simon's submit script
data_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/'
write_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/models/run_2/simon_trim_20'

submit_mrbrt <- function(bundle){
  
  # Make folder if it doesn't exist
  bun_write_folder <- paste0(write_folder, '/', bundle)
  if(!file.exists(bun_write_folder)){
    dir.create(bun_write_folder)
  }
  
  spec_file <- paste0('no_inf/', bundle,'.csv')
  command<-paste("/bin/bash -c", shQuote(paste0("source /ihme/code/evidence_score/miniconda3/bin/activate mr_brt_env; python /ihme/code/evidence_score/mr_brt_ihme/uw_amo/mr_brt.py \\", ##sy: run script
                                                "--input_dir ", data_folder, " \\",
                                                "--data_file ", spec_file, " \\",
                                                "--cov_metadata_file mr_brt_metadata_linear_with_otp.csv \\",
                                                "--output_dir ", bun_write_folder, " \\",
                                                "--effect_variable log_mean \\",
                                                "--effect_se_variable log_se \\",
                                                "--optimization_method trim_maxL \\",
                                                "--trim_pct 0.2")))
  
  # Run it
  tryCatch(
    {
      system(command)
    },
    error = function(e){
      message("Underlying MR BRT code failed with :")
      message(e)
      message("Checking for model files, if these exist then only plotting has failed")
      
    }
    # finally = {
    #   if(all(file.exists(file_list))){
    #     message(" All model files exist")
    #   }else{
    #     stop("Missing model files")
    #   }
    # }
  )
  
}

