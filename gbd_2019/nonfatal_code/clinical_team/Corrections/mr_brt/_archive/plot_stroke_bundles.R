#############################################################
##Author: Nick Roberts
##Date: 4/1/2019
##Purpose: Compare acute CF2 and CF3 for CVD team
##Notes: 
##Updates: 
#
###########################################################
library(ggplot2)
library(gridExtra)
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

bundles <- c(116, 371, 454, 455)

train_data_files_cf2 <- Sys.glob('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf2/estimates/corrections/models/bundle/cf2/*/mr_brt_mod/train_data.csv')
results_files_cf2 <- Sys.glob('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf2/estimates/corrections/models/bundle/cf2/*/mr_brt_mod/model_summaries.csv')

train_data_files_cf3 <- Sys.glob('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf3/estimates/corrections/models/bundle/cf3/*/mr_brt_mod/train_data.csv')
results_files_cf3 <- Sys.glob('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf3/estimates/corrections/models/bundle/cf3/*/mr_brt_mod/model_summaries.csv')

pdf('/home/j/temp/nrober75/cvd_cf_plots_trim_uncertainty.pdf', onefile = TRUE, compress = TRUE)
for(bundle in bundles){
  train_data<- fread(paste0('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf2/estimates/corrections/models/bundle/cf2/', bundle, '/mr_brt_mod/train_data.csv')) 
  model_sum <- fread(paste0('/share/hospital/clinical_runs/first_run_new_cluster_bundle_cf2/estimates/corrections/models/bundle/cf2/', bundle, '/mr_brt_mod/model_summaries.csv')) %>%
    setnames(c('X_sex_id', 'X_age_start'), c('sex_id', 'age_start'))
  
  
  bun_df <- loadBundles(bundles)
  bundle_name <- bun_df[bundle_id == bundle]$bundle_name
  train_data[w < 0.5, trimmed := 'Trimmed data'][w >= 0.5, trimmed := 'Untrimmed data']
  train_data[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
  
  train_data[, source := factor(source, levels = c('Marketscan', 'HCUP', 'NZL', 'PHL'))]
  data_plot_cf2 <- ggplot(train_data, aes(x = age_start, y = exp(log_mean), color = factor(trimmed))) + 
    geom_point() + theme_bw() +
    geom_line(data = model_sum, aes(y = exp(Y_mean), color = 'MR-BRT predictions'), size = 2) +
    ylab('CF2') + xlab('Age') + 
    ggtitle(paste0(bundle_name, ' CF2 inputs and results')) +
    scale_color_discrete(name = '') +
    facet_wrap(~sex, ncol = 2, scales = 'free') +
    theme(plot.title = element_text(size = 9))
  
  
  if('sex_id' %in% names(model_sum)){
    model_sum[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
  } else{
    model_sum[, sex := 1] # Arbitrary. Should be obvious from the bundle itself
  }
  model_sum[, hi_hi := max(Y_mean_hi), by = c('age_start', 'sex')]
  model_sum[, lo_lo := min(Y_mean_lo), by = c('age_start', 'sex')]
  results_plot_cf2 <- ggplot(unique(model_sum), aes(x = age_start, y = exp(Y_mean))) +
    geom_point() + theme_bw() +
    ylab('CF2 model predictions') + xlab('Age') +
    ggtitle(paste0(bundle_name, ' MR-BRT predictions - CF2')) +
    facet_wrap(~sex) +
    theme(plot.title = element_text(size = 9))
  
  
  
  ## CF3
  train_data_cf3 <- fread(paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/models/bundle/cf3/', bundle, '/mr_brt_mod/train_data.csv'))
  model_sum_cf3 <- fread(paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/models/bundle/cf3/', bundle, '/mr_brt_mod/model_summaries.csv')) %>%
    setnames(c('X_sex_id', 'X_age_start'), c('sex_id', 'age_start'))
  
  
  bun_df <- loadBundles(bundles)
  bundle_name <- bun_df[bundle_id == bundle]$bundle_name
  train_data_cf3[w < 0.5, trimmed := 'Trimmed data'][w >= 0.5, trimmed := 'Untrimmed data']
  train_data_cf3[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']

  data_plot_cf3 <- ggplot(train_data_cf3, aes(x = age_start, y = exp(log_mean), color = factor(trimmed))) + 
    geom_point() + theme_bw() +
    geom_line(data = model_sum_cf3, aes(y = exp(Y_mean), color = 'MR-BRT predictions'), size = 2) +
    ylab('CF3 input data') + xlab('Age') + 
    ggtitle(paste0(bundle_name, ' CF3 inputs and results')) +
    scale_color_discrete(name = '') +
    facet_wrap(~sex, ncol = 2, scales = 'free') +
    theme(plot.title = element_text(size = 9))
  
  if('sex_id' %in% names(model_sum)){
    model_sum_cf3[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
  } else{
    model_sum_cf3[, sex := 1] # Arbitrary. Should be obvious from the bundle itself
  }
  model_sum_cf3[, hi_hi := max(Y_mean_hi), by = c('age_start', 'sex')]
  model_sum_cf3[, lo_lo := min(Y_mean_lo), by = c('age_start', 'sex')]
  model_sum_cf3[, cf := 'CF3']
  model_sum[, cf := 'CF2']
  results_plot_cf3 <- ggplot(unique(model_sum_cf3), aes(x = age_start, y = exp(Y_mean))) +
    geom_line(aes(color = 'CF3 predictions')) + theme_bw() +
    geom_line(data = unique(model_sum), aes(color = 'CF2 predictions')) +
    ylab('Model predictions') + xlab('Age') +
    ggtitle(paste0(bundle_name, ' MR-BRT predictions')) +
    facet_wrap(~sex) +
    theme(plot.title = element_text(size = 9)) +
    geom_ribbon(data = model_sum_cf3, aes(ymin = exp(lo_lo), ymax = exp(hi_hi), fill = 'CF3 prediction'), alpha = 0.5) +
    geom_ribbon(data = model_sum, aes(ymin = exp(lo_lo), ymax = exp(hi_hi), fill = 'CF2 prediction'), alpha = 0.5)
  
  
  # p <- grid.arrange(data_plot_cf2, data_plot_cf3, results_plot_cf3, ncol = 2)
  # print(p)
  p <- results_plot_cf3
  data <- grid.arrange(data_plot_cf2, data_plot_cf3, ncol = 1)
  print(p)
  print(data)
  
}
dev.off()

bun <- str_split(train_data_files[i], '/')[[1]][11]

train_data <- fread(train_data_files[i])
model_sum <- fread(paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/models/bundle/cf2/', bun, '/mr_brt_mod/model_summaries.csv')) 

if('X_sex_id' %in% names(model_sum)){
  setnames(model_sum, c('X_sex_id', 'X_age_start'), c('sex_id', 'age_start'))
} else{
  setnames(model_sum, 'X_age_start', 'age_start')
}

bun <- str_split(train_data_files[i], '/')[[1]][11]
bun_df <- loadBundles(bun)
bundle_name <- bun_df[bundle_id == bun]$bundle_name[1]

# # Weir dformatting
# train_data <- train_data %<>% mutate_if(is.numeric, as.character)
# 
# plot_data <- merge(train_data, model_sum[, c('sex_id', 'age_start', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')],
#                    by = c('sex_id', 'age_start'), all.x = TRUE)

train_data[w < 0.5, trimmed := 'Trimmed data'][w >= 0.5, trimmed := 'Untrimmed data']
train_data[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
train_data <- merge(train_data, locs[, c('location_id', 'parent_id', 'location_name')])
train_data[parent_id == 9, label := 'PHL'][parent_id == 102 | parent_id == 104, label := 'HCUP SIDS'][parent_id == 72, label := 'NZL']
data_plot <- ggplot(train_data, aes(x = age_start, y = exp(log_mean), color = factor(trimmed))) + 
  geom_point() + theme_bw() +
  ylab('CF3 input data') + xlab('Age') + 
  ggtitle(paste0(bundle_name, ' input data')) +
  scale_color_discrete(name = '') +
  facet_wrap(label~sex, ncol = 2, scales = 'free')

if('sex_id' %in% names(model_sum)){
  model_sum[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
} else{
  model_sum[, sex := 1] # Arbitrary. Should be obvious from the bundle itself
}
model_sum[, hi_hi := max(Y_mean_hi), by = c('age_start', 'sex')]
model_sum[, lo_lo := min(Y_mean_lo), by = c('age_start', 'sex')]
results_plot <- ggplot(unique(model_sum), aes(x = age_start, y = exp(Y_mean))) +
  geom_point() + theme_bw() +
  ylab('CF3 model predictions') + xlab('Age') +
  ggtitle(paste0(bundle_name, ' MR-BRT predictions')) +
  facet_wrap(~sex)
