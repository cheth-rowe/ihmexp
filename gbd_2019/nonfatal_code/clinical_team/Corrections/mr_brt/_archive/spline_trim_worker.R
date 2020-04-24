#############################################################
##Author: Nick Roberts
##Date: 3/12/2019
##Purpose: Master script to submit updated MR-BRT models
##Notes: Parallelize your submission of the model
##Updates: 3/12: Run with 20% trim with new code
##Notes2: Use same data as before because that hasn't changed and still has covariates
##Notes3: Update to make predictions for all hospital data places
###########################################################
## Child script to actually run it

# Get arguments from qsub
rm(list = ls())
print(commandArgs())

bundle <- commandArgs()[5]
run <- commandArgs()[6]
#bundle <- 3039
library(data.table)
library(parallel)
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RMySQL)
source('/ihme/cc_resources/libraries/current/r/get_location_metadata.R')

locs <- get_location_metadata(35)
# Libraries and sources from Reed
repo_dir <- "/home/j/temp/reed/prog/projects/run_mr_brt/"
source(paste0(repo_dir, "run_mr_brt_function.R"))
source(paste0(repo_dir, "cov_info_function.R"))
source(paste0(repo_dir, "check_for_outputs_function.R"))
source(paste0(repo_dir, "load_mr_brt_outputs_function.R"))
#source(paste0(repo_dir, "predict_mr_brt_function.R"))
source('/homes/nrober75/code/predict_mr_brt_nrober75.R')
source(paste0(repo_dir, "check_for_preds_function.R"))
source(paste0(repo_dir, "load_mr_brt_preds_function.R"))

loadBundles <- function(bundle_ids) {
  # Get bundles for which we have scalars
  
  db_con = fread(paste0("/home/j/temp/hospital/2016/code/db_con.csv"))
  
  # Get list of ids and names
  con <- dbConnect(dbDriver("MySQL"),
                   username = db_con$username,
                   password = db_con$pass,
                   host = db_con$host)
  df <- dbGetQuery(con, sprintf("
                                SELECT
                                b.bundle_id,
                                b.bundle_name,
                                c.cause_name,
                                c.sort_order
                                FROM
                                bundle.bundle b
                                LEFT JOIN
                                shared.cause_hierarchy_history c ON b.cause_id = c.cause_id AND c.cause_set_version_id = shared.active_cause_set_version(2, 4)
                                WHERE
                                b.bundle_id IN (%s)",
                                paste0(bundle_ids, collapse = ",")))
  dbDisconnect(con)
  df <- data.table(df)
  df <- df[is.na(sort_order), sort_order := 9999]
  df <- df[order(sort_order)][, sort_order := NULL]
  return(df)
}
bun_df <- loadBundles(bundle)
bundle_name <- bun_df[bundle_id == bundle]$bundle_name
# Folders where your stuff is
data_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/'
write_folder <- paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/models/run_', run)
dir.create(write_folder)

# Make bundle-specific write folder if it doesn't already exist
bun_write_folder <- paste0(write_folder, '/', bundle)
if(!file.exists(bun_write_folder)){
  dir.create(bun_write_folder)
}

# Make a folder to store pdfs and trim_pdfs
if(!file.exists(paste0(write_folder, '/pdfs'))){
  dir.create(paste0(write_folder, '/pdfs'))
  dir.create(paste0(write_folder, '/trim_pdfs'))
  
}

# # Run the MR-BRT model
# fit10 <- run_mr_brt(
#   output_dir = bun_write_folder, 
#   model_label = paste0(bundle, '_run_trim_10'),
#   data = paste0(data_folder, bundle, '.csv') ,
#   mean_var = "log_mean",
#   se_var = "log_se",
#   covs = list( 
#     #cov_info("zcov1", "Z"),
#     cov_info("op_env", "X"),
#     cov_info("ip_env", "X"),
#     cov_info("haqi", "X"),
#     cov_info("sex_id", "X"),
#     cov_info("age_start", "X") ),
#   overwrite_previous = TRUE,
#   project = 'proj_hospital',
#   method = 'trim_maxL', trim_pct = .1
# )

fit20 <- run_mr_brt(
  output_dir = bun_write_folder, 
  model_label = paste0(bundle, '_run_trim_20'),
  data = paste0(data_folder, bundle, '.csv') ,
  mean_var = "log_mean",
  se_var = "log_se",
  covs = list( 
    #cov_info("zcov1", "Z"),
    #cov_info("op_env", "X"),
    #cov_info("ip_env", "X"),
    #cov_info("haqi", "X"),
    cov_info("sex_id", "X"),
    cov_info('haqi', 'X',
             prior_mean = 0, 
             prior_var = 'inf',
             bspline_prior_mean = "0, 0, 0", 
             bspline_prior_var = "inf, inf, inf",
             degree = 3,
             n_i_knots = 2, 
             r_linear = TRUE,
             l_linear = TRUE), 
    cov_info('age_start', 'X',
             prior_mean = 0, 
             prior_var = 'inf',
             bspline_prior_mean = "0, 0, 0", 
             bspline_prior_var = "inf, inf, inf",
             degree = 3,
             n_i_knots = 2, 
             r_linear = TRUE,
             l_linear = TRUE) ),
  overwrite_previous = TRUE,
  project = 'proj_hospital',
  method = 'trim_maxL', trim_pct = .25
)

# fit30 <- run_mr_brt(
#   output_dir = bun_write_folder, 
#   model_label = paste0(bundle, '_run_trim_30'),
#   data = paste0(data_folder, bundle, '.csv') ,
#   mean_var = "log_mean",
#   se_var = "log_se",
#   covs = list( 
#     #cov_info("zcov1", "Z"),
#     cov_info("op_env", "X"),
#     cov_info("ip_env", "X"),
#     cov_info("haqi", "X"),
#     cov_info("sex_id", "X"),
#     cov_info("age_start", "X") ),
#   overwrite_previous = TRUE,
#   project = 'proj_hospital',
#   method = 'trim_maxL', trim_pct = .3
# )


# # Check for outputs and continue after. Check 40 first becuase it's the last one, should theoretically finish last
# # Then check for others
# for(i in 1:100){
#   i <- i+1
#   a <- check_for_outputs(fit30, 45)
#   if(a == TRUE){
#     break
#   }
# }

for(i in 1:1000){
  i <- i+1
  print(i)
  a <- check_for_outputs(fit20, 40)
  if(a == TRUE){
    break
  }
}

# for(i in 1:100){
#   i <- i+1
#   a <- check_for_outputs(fit10, 10)
#   if(a == TRUE){
#     break
#   }
# }

#preddf <- fread(paste0(data_folder, bundle, '.csv')) # Training data
#preddf <- fread(paste0('/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/no_inf/', bundle, '.csv'))
preddf <- fread('/share/hospital/scratch/marketscan/mr_brt_df_split/marketscan_prepped.csv')
preddf <- preddf[bundle_id == bundle & cf3 != Inf & !is.na(cf3)]
## Download predictions dataset
#preddf <- fread('/share/hospital/scratch/marketscan/mr_brt_df_split/preddf.csv')
#pred10 <- predict_mr_brt(fit10, newdata = preddf) # This takes ages
pred20 <- predict_mr_brt(fit20, newdata = preddf) # This takes ages
#pred30 <- predict_mr_brt(fit30, newdata = preddf) # This takes ages

# Check for preds, same order as above
# for(i in 1:1000){
#   i <- i+1
#   a <- check_for_preds(pred30, 45)
#   if(a == TRUE){
#     break 
#   }
# }
for(i in 1:1000){
  i <- i+1
  print(i)
  a <- check_for_preds(pred20, 20)
  if(a == TRUE){
    break
  }
}
# for(i in 1:1000){
#   i <- i+1
#   a <- check_for_preds(pred10, 10)
#   if(a == TRUE){
#     break
#   }
# }

#pred_object10 <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_10/prediction_newdata.csv'))
# Now merge and plot
# pred_object10 <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_10/model_summaries.csv')) %>%
#   setnames(c('X_op_env', 'X_ip_env', 'X_haqi', 'X_sex_id', 'X_age_start', 'Y_mean'), c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start', 'Y_mean_10'))

## All covs
# pred_object <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_summaries.csv')) %>%
#   setnames(c('X_op_env', 'X_ip_env', 'X_haqi', 'X_sex_id', 'X_age_start'), c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start'))

## Select covs (fill in from above based on what covariates you are using)
pred_object <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_summaries.csv')) %>%
  setnames(c('X_haqi', 'X_sex_id', 'X_op_env', 'X_age_start'), c('haqi', 'sex_id', 'op_env', 'age_start'))


pred_object <- pred_object %<>% mutate_if(is.numeric, as.character)

preddf <- preddf[, c('age_start', 'location_id', 'year_start', 'sex_id', 'bundle_id', 'cf3', 'se', 'ip_env', 'haqi', 'sdi', 'op_env')]
preddf[, log_mean := log(cf3)]
preddf <- preddf %>% mutate_if(is.numeric, as.character)

pred_results <- merge(preddf, pred_object[, c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')],
                      by = c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start'), all.x = TRUE)

pred_results <- merge(preddf, pred_object[, c('haqi', 'sex_id', 'op_env', 'age_start', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')],
                      by = c('haqi', 'sex_id', 'age_start', 'op_env'), all.x = TRUE)


pred_results <- as.data.table(pred_results)
pred_results[, age_start := as.numeric(age_start)][, Y_mean := as.numeric(Y_mean)]

# Plot
pred_results[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
pred_results[, location_id := as.numeric(location_id)]
pred_results[, haqi := as.numeric(haqi)]
pred_results <- merge(pred_results, locs[, c('location_id', 'super_region_name', 'region_name')])

# Calculate upper and lower at each age for the plot
pred_results[, Y_mean_lo := as.numeric(Y_mean_lo)][, Y_mean_hi := as.numeric(Y_mean_hi)]
pred_results[, lo_lo := min(Y_mean_lo), by = c('sex', 'age_start')]
pred_results[, hi_hi := max(Y_mean_hi), by = c('sex', 'age_start')]
pred_results[location_id >= 523 & location_id <= 573, usa_hi := max(Y_mean_hi), by = c('sex', 'age_start')]

# Get rid of crazy high uncertainty
pred_results[, test_hi := quantile(Y_mean_hi, 0.975), by = c('sex', 'age_start')]

# usa plot
model_df <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/train_data.csv'))
model_df[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
model_df[w >= 0.5, trim := 'Untrimmed'][w <= 0.5, trim := 'Trimmed'][, trim := factor(trim, levels = c('Trimmed', 'Untrimmed'))]

train_df <-fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/train_data.csv'))

# p1 <- ggplot(pred_results[location_id >= 523 & location_id <= 573], aes(x = age_start, y = exp(Y_mean))) +
#   geom_point(data = model_df[w == 1], aes(color = 'Used Training USA data', y = exp(log_mean)), alpha = 0.5) +
#   geom_point(data = model_df[w == 0], aes(color = 'Trimmed USA data', y = exp(log_mean)), alpha = 0.5) +
#   geom_point(aes(color = 'MR-BRT prediction'), alpha = 1) +
#   facet_wrap(~sex) +
#   theme_bw() +
#   ggtitle(paste0('MR-BRT predictions with all data; USA \n', bundle_name)) +
#   scale_color_discrete(name = "") +
#   ylab('CF3 prediction') + xlab('Age')

p2 <- ggplot(model_df, aes(y = exp(log_mean), x = age_start, color = trim)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(name = "", values = c('gold', 'purple')) +
  facet_wrap(~sex) +
  theme_bw() +
  ggtitle(paste0(bundle_name, '\nData used in the model; USA')) +
  ylab('CF3 prediction') + xlab("Age")

p3 <- ggplot(pred_results[location_id >= 523 & location_id <= 573], aes(x = age_start, y = exp(Y_mean))) +
  geom_point(aes(color = 'MR-BRT prediction'), alpha = 0.5) +
  facet_wrap(~sex) +
  theme_bw() +
  ggtitle('USA CF3 MR-BRT results; 25% trimmed')  +
  scale_color_discrete(name = "") +
  ylab('CF3 prediction') + xlab("Age") +
  #geom_ribbon(aes(ymax = exp(usa_hi), ymin = exp(lo_lo), fill = 'MR-BRT prediction'), alpha = 0.2) +
  scale_fill_discrete(name = '')
#coord_cartesian(ylim = c(0, exp(max_val)*1.25))
#ylim(0, y_lim*3)

coefs <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_coefs.csv'))
coefs <- coefs[, 1:3]
names(coefs) <- c('XCov', 'Beta', 'Beta var')
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(coefs, rows = NULL, theme = tt)

# Print
pdf(paste0(write_folder, '/pdfs', '/', bundle, '.pdf'))
grid.arrange(p2, p3, tbl, heights = c(2,2,1.5))
dev.off()



# 
# ten <- ggplot(pred_results, aes(x = age_start, y = exp(Y_mean_10))) + 
#   geom_point() +
#   ggtitle('10% trim model') +
#   ylab('10% Trim CF3') +
#   theme_bw() +
#   geom_smooth()
# twenty <- ggplot(pred_object20, aes(x = age_start, y = exp(Y_mean_20))) + 
#   geom_point() +
#   ggtitle('20% trim model') +
#   ylab('20% Trim CF3') +
#   theme_bw() +
#   geom_smooth()
# 
# thirty <- ggplot(pred_object30, aes(x = age_start, y = exp(Y_mean_30))) + 
#   facet_wrap()
#   geom_point() +
#   ggtitle('30% trim model') +
#   ylab('30% Trim CF3') +
#   theme_bw() +
#   geom_smooth()
# 
# ## Now coefs for each
# coefs10 <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_10/model_coefs.csv'))
# coefs10 <- coefs[, 1:3]
# names(coefs10) <- c('XCov', 'Beta', 'Beta var')
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
# tbl10 <- tableGrob(coefs10, rows = NULL, theme = tt)
# 
# coefs20 <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_coefs.csv'))
# coefs20 <- coefs[, 1:3]
# names(coefs20) <- c('XCov', 'Beta', 'Beta var')
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
# tbl20 <- tableGrob(coefs20, rows = NULL, theme = tt)
# 
# coefs30 <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_30/model_coefs.csv'))
# coefs30 <- coefs[, 1:3]
# names(coefs30) <- c('XCov', 'Beta', 'Beta var')
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
# tbl30 <- tableGrob(coefs30, rows = NULL, theme = tt)
# 
# pdf(paste0(write_folder, '/pdfs', '/', bundle, '.pdf'))
# grid.arrange(ten, tbl10, twenty, tbl20, thirty, tbl30)
# dev.off()

# 
# 
# pdf(paste0(bun_write_folder, '/pdfs/', bundle, '.pdf'))
# all <- grid.arrange(p1, p3, tbl, nrow = 3, as.table = TRUE, heights = c(1.7,1.7,1.3))
# dev.off()
# 
# # Just trim data points
# pdf(paste0(bun_write_folder, '/trim_pdfs/', bundle, '.pdf'))
# all <- grid.arrange(p1, p2, p3, nrow = 3, as.table = TRUE, heights = c(1.5,1.5,1.5))
# dev.off()
# 
# 
# 
# 
# #pred_object <- load_mr_brt_preds(pred1)
# pred_object <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_summaries.csv'))
# setnames(pred_object, c('X_op_env', 'X_ip_env', 'X_haqi', 'X_sex_id', 'X_age_start'), c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start'))
# 
# cols <- c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start')
# 
# pred_object <- pred_object %<>% mutate_if(is.numeric, as.character)
# preddf <- preddf %>% mutate_if(is.numeric, as.character)
# 
# pred_results <- merge(preddf, pred_object[, c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')],
#                       by = c('op_env', 'ip_env', 'haqi', 'sex_id', 'age_start'), all.x = TRUE)
# pred_results <- as.data.table(pred_results)
# pred_results[, age_start := as.numeric(age_start)][, Y_mean := as.numeric(Y_mean)]
# 
# # Plot
# pred_results[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
# pred_results[, location_id := as.numeric(location_id)]
# pred_results[, haqi := as.numeric(haqi)]
# pred_results <- merge(pred_results, locs[, c('location_id', 'super_region_name', 'region_name')])
# 
# # Calculate upper and lower at each age for the plot
# pred_results[, Y_mean_lo := as.numeric(Y_mean_lo)][, Y_mean_hi := as.numeric(Y_mean_hi)]
# pred_results[, lo_lo := min(Y_mean_lo), by = c('sex', 'age_start')]
# pred_results[, hi_hi := max(Y_mean_hi), by = c('sex', 'age_start')]
# pred_results[location_id >= 523 & location_id <= 573, usa_hi := max(Y_mean_hi), by = c('sex', 'age_start')]
# 
# # Get rid of crazy high uncertainty
# pred_results[, test_hi := quantile(Y_mean_hi, 0.975), by = c('sex', 'age_start')]
# 
# # usa plot
# model_df <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/train_data.csv'))
# model_df[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
# 
# 
# 
# # One way to show the full run and trimming, not really necessary
# # p1 <- ggplot(pred_results[location_id >= 523 & location_id <= 573], aes(x = age_start, y = exp(Y_mean))) + 
# #   geom_point(data = model_df[w == 1], aes(color = 'Used Training USA data', y = exp(log_mean)), alpha = 0.5) +
# #   geom_point(data = model_df[w == 0], aes(color = 'Trimmed USA data', y = exp(log_mean)), alpha = 0.5) +
# #   geom_point(aes(color = 'MR-BRT prediction'), alpha = 1) + 
# #   facet_wrap(~sex) + 
# #   theme_bw() + 
# #   ggtitle(paste0('MR-BRT predictions with all data; USA \n', bundle_name)) +
# #   scale_color_discrete(name = "") +
# #   ylab('CF3 prediction') + xlab('Age')
# # p2 <- ggplot(model_df[w == 0], aes(y = exp(log_mean), x = age_start)) + 
# #   scale_color_manual(name = "", values = '#00BA38') +
# #   geom_point(aes(color = 'Trimmed USA data')) +
# #   facet_wrap(~sex) + 
# #   theme_bw() +
# #   ggtitle('Trimmed data from the model; USA') +
# #   ylab('CF3 prediction') + xlab("Age")
# # 
# # max_val <- max(pred_results[location_id >= 523 & location_id <= 573]$Y_mean)
# # 
# # pred_results[location_id >= 523 & location_id <= 573 & exp(hi_hi) > 1000*exp(max_val)]
# # p3 <- ggplot(pred_results[location_id >= 523 & location_id <= 573], aes(x = age_start, y = exp(Y_mean))) + 
# #   geom_point(aes(color = 'MR-BRT prediction'), alpha = 0.5) + 
# #   facet_wrap(~sex) + 
# #   theme_bw() +
# #   ggtitle('USA CF3 MR-BRT results; 20% trimmed')  +
# #   scale_color_discrete(name = "") +
# #   ylab('CF3 prediction') + xlab("Age") +
# #   geom_ribbon(aes(ymax = exp(usa_hi), ymin = exp(lo_lo), fill = 'MR-BRT prediction'), alpha = 0.2) +
# #   scale_fill_discrete(name = '')
# #   coord_cartesian(ylim = c(0, exp(max_val)*1.25))
# #   #ylim(0, y_lim*3)
# #   
# 
# coefs <- fread(paste0(bun_write_folder, '/', bundle, '_run_trim_20/model_coefs.csv'))
# coefs <- coefs[, 1:3]
# names(coefs) <- c('XCov', 'Beta', 'Beta var')
# tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
# tbl <- tableGrob(coefs, rows = NULL, theme = tt)
# 
# 
# pdf(paste0(bun_write_folder, '/pdfs/', bundle, '.pdf'))
# all <- grid.arrange(p1, p3, tbl, nrow = 3, as.table = TRUE, heights = c(1.7,1.7,1.3))
# dev.off()
# 
# # Just trim data points
# pdf(paste0(bun_write_folder, '/trim_pdfs/', bundle, '.pdf'))
# all <- grid.arrange(p1, p2, p3, nrow = 3, as.table = TRUE, heights = c(1.5,1.5,1.5))
# dev.off()
# 
# 
# # 
# # 
# # 
# # trim_20 <- pred_object$model_summaries %>% setnames('Y_mean', 'trim_20_Y_mean')
# # trim_20 <- as.data.table(trim_20)
# # 
# # write_csv(trim_20, paste0(write_folder, ))
