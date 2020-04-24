###############################################################
##Author: Nick Roberts
##Date: 3/4/2019
##Purpose: Make plots for custom correction requests
##Notes: 
##Updates: 
#
###########################################################
library(gridExtra)
library(boot)
preddf <- fread('/share/hospital/clinical_runs/run_1/estimates/corrections/age_sex_preddf.csv')
preddf[age_start == 110, age_start := 97.5]
bun_df <- fread('/share/hospital/scratch/covariates/bundle_names.csv')
restrictions <- get_icg_restrictions()
bundle_map <- bundle_icg() %>% merge(restrictions[map_version == 23], by = c('icg_id', 'icg_name')) %>% .[, c('bundle_id', 'male', 'female', 'yld_age_start', 'yld_age_end')]
bundle_map[, min_age := min(yld_age_start), by = 'bundle_id'][, max_age := max(yld_age_end), by = 'bundle_id']
bundle_map <- unique(bundle_map[, c('bundle_id', 'min_age', 'max_age')])
bundle_map[max_age == 95, max_age := 110]
custom_bundles <- fread('/home/j/temp/nrober75/Custom Correction Requests (4_9) - Form Responses 1 (1).csv')
#custom_bundles <- custom_bundles[`Bundle ID (requests should be submitted for one bundle at a time)` != 4196] # Do separately, shouldn't show up for GBD 2017
#custom_bundles <- custom_bundles$`Bundle ID (requests should be submitted for one bundle at a time)`

cf3_only <- unique(custom_bundles[Figure == 'Age pattern CF3 GBD 2017 vs GBD 2019']$`Bundle ID (requests should be submitted for one bundle at a time)`)
cf2_too <- unique(custom_bundles[Figure == 'Age pattern CF2 fit, CF3 fit, USA CF3 data, TWN CF3 data']$`Bundle ID (requests should be submitted for one bundle at a time)`)
cf2_too <- c(cf2_too, 278) # Guillan Barre too

cfs <- c('cf1', 'cf2', 'cf3')
bun_df <- fread('/share/hospital/scratch/covariates/bundle_names.csv')

pdf('/home/j/temp/nrober75/cf3_cf2_requests_TLs.pdf')
lapply(cf2_too[1:length(cf2_too)], function(bundle){
  
  train_data <- fread(paste0('/share/hospital/clinical_runs/run_2_cf3_with_twn/estimates/corrections/models/bundle/', 'cf3', '/', bundle, '/mr_brt_mod/train_data.csv'))
  train_data[w < 0.5, trimmed := 'Trimmed data'][w >= 0.5, trimmed := 'Untrimmed data']
  train_data[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
  train_data[location_id == 8, location := 'Taiwan'][location_id != 8, location := 'Marketscan']
  min_age <- min(train_data[w == 1]$age_start) 
  max_age <- max(train_data[w == 1]$age_start)
  data_plot <- ggplot(train_data, aes(x = age_start, y = exp(log_mean), color = trimmed)) + 
    geom_point(size = 1.8) + theme_bw() +
    ylab('GBD 2019 CF3 input data') + xlab('Age') + 
    ggtitle('CF3 input data') +
    scale_color_manual(name = '', values = c('#33cc33', '#D354AC')) +
    facet_wrap(~sex, ncol = 4, scales = 'free_y') +
    xlim(c(min_age, max_age))
  

  
  print(bundle)
  restrictions_map <- bundle_map[bundle_id == bundle]
  
  a <- tryCatch(all_cfs <- rbindlist(lapply(cfs, function(cf){
    #print(cf)
    dt <- fread(paste0('/share/hospital/clinical_runs/run_2_cf3_with_twn/estimates/corrections/models/bundle/', cf, '/', bundle, '/mr_brt_mod/model_summaries.csv'))
    dt[, cf_num := cf]
    dt[cf_num == 'cf1', Y_mean := inv.logit(Y_mean)][cf_num == 'cf1', Y_mean_lo := inv.logit(Y_mean_lo)][cf_num == 'cf1', Y_mean_hi := inv.logit(Y_mean_hi)]
    dt[cf_num != 'cf1', Y_mean := exp(Y_mean)][cf_num != 'cf1', Y_mean_lo := exp(Y_mean_lo)][cf_num != 'cf1', Y_mean_hi := exp(Y_mean_hi)]
    dt[, bundle_id := bundle]
    return(dt)
  }), fill = TRUE))
  
  if(bundle == 198){
    all_cfs[is.na(X_sex_id), X_sex_id := 1]
    all_cfs <- all_cfs[X_sex_id != 2]
  }

  all_cfs <- all_cfs[X_age_start >= min_age & X_age_start <= max_age]
  
  # Propagage out ages
  limits <- all_cfs[X_age_start == min_age | X_age_start == max_age]
  limits <- limits[, c('X_sex_id', 'X_age_start', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi', 'cf_num')]
  setnames(limits, 'X_age_start', 'limit')

  plot_preds <- preddf[age_start <= min_age | age_start >= max_age]
  plot_preds[age_start <= min_age, limit := min_age][age_start >= max_age, limit := max_age]
  setnames(plot_preds, 'sex_id', 'X_sex_id')
  extend <- merge(limits, plot_preds, by = c('limit', 'X_sex_id'), allow.cartesian = TRUE)
  setnames(extend, 'age_start', 'X_age_start')

  all_cfs <- rbind(all_cfs, extend, fill = TRUE)
  all_cfs[X_sex_id == 1, sex := 'Males'][X_sex_id == 2, sex := 'Females']
  
  
  rm(past_cf3, past_cf2, past_cf1, past_dt)
  b <- tryCatch({
    past_cf3 <- fread(paste0('/share/hospital/scratch/hospital/cf_draws/cf3_by_bundle/bundle_',bundle, '.csv')) %>%
                  .[location_id == 64] %>% 
                  melt(id.vars = c('location_id', 'sex_id', 'age_start', 'bundle_id', 'model_prediction'), 
                       value.vars = patterns('prevalence_')) %>%
                  .[, Y_mean := mean(value), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, Y_mean_lo := quantile(value, probs = 0.025), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, Y_mean_hi := quantile(value, probs = 0.975), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, c('location_id', 'age_start', 'sex_id', 'model_prediction', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')] %>% unique() %>% 
                  .[sex_id == 1, sex := 'Males'] %>%
                  .[sex_id == 2, sex := 'Females'] %>%
                  .[age_start >= restrictions_map$min_age[1] & age_start <= restrictions_map$max_age] %>%
                  .[, cf_num := 'cf2'] %>%
                  .[, c('location_id', 'age_start', 'sex_id', 'model_prediction', 'Y_mean_lo', 'Y_mean_hi')] %>%
                  setnames('model_prediction', 'Y_mean') %>%
                  .[sex_id == 1, sex := 'Males'] %>%
                  .[sex_id == 2, sex := 'Females'] %>%
                  .[age_start >= restrictions_map$min_age[1] & age_start <= restrictions_map$max_age] %>%
                  .[, cf_num := 'cf3']
                
                past_cf2 <- fread(paste0('/share/hospital/scratch/hospital/cf_draws/cf2_by_bundle/bundle_',bundle, '.csv')) %>%
                  .[location_id == 102] %>%
                  melt(id.vars = c('location_id', 'sex_id', 'age_start', 'bundle_id'), value.vars = patterns('incidence_')) %>%
                  .[, Y_mean := mean(value), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, Y_mean_lo := quantile(value, probs = 0.025), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, Y_mean_hi := quantile(value, probs = 0.975), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                  .[, c('location_id', 'age_start', 'sex_id', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')] %>% unique() %>% 
                  .[sex_id == 1, sex := 'Males'] %>%
                  .[sex_id == 2, sex := 'Females'] %>%
                  .[age_start >= restrictions_map$min_age[1] & age_start <= restrictions_map$max_age] %>%
                  .[, cf_num := 'cf2']
                
                if(bundle != 19){
                  past_cf1 <- fread(paste0('/share/hospital/scratch/hospital/cf_draws/cf1_by_bundle/bundle_',bundle, '.csv')) %>%
                  .[location_id == 102] %>%
                    melt(id.vars = c('location_id', 'sex_id', 'age_start', 'age_end', 'bundle_id'), value.vars = patterns('indv_cf_')) %>%
                    .[, Y_mean := mean(value), by = c('location_id', 'sex_id', 'age_start', 'age_end', 'bundle_id')] %>%
                    .[, Y_mean_lo := quantile(value, probs = 0.025), by = c('location_id', 'sex_id', 'age_start', 'age_end', 'bundle_id')] %>%
                    .[, Y_mean_hi := quantile(value, probs = 0.975), by = c('location_id', 'sex_id', 'age_start', 'age_end', 'bundle_id')] %>%
                    .[, c('location_id', 'age_start', 'sex_id', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')] %>% unique() %>% 
                    .[sex_id == 1, sex := 'Males'] %>%
                    .[sex_id == 2, sex := 'Females'] %>%
                    .[age_start >= restrictions_map$min_age[1] & age_start <= restrictions_map$max_age] %>%
                    .[, cf_num := 'cf1']
                } else{
                  past_cf1 <- fread(paste0('/share/hospital/scratch/hospital/cf_draws/cf1_by_bundle/bundle_',bundle, '.csv')) %>%
                    .[location_id == 102] %>%
                    melt(id.vars = c('location_id', 'sex_id', 'age_start', 'bundle_id'), value.vars = patterns('prevalence_')) %>%
                    .[, Y_mean := mean(value), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                    .[, Y_mean_lo := quantile(value, probs = 0.025), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                    .[, Y_mean_hi := quantile(value, probs = 0.975), by = c('location_id', 'sex_id', 'age_start', 'bundle_id')] %>%
                    .[, c('location_id', 'age_start', 'sex_id', 'Y_mean', 'Y_mean_lo', 'Y_mean_hi')] %>% unique() %>% 
                    .[sex_id == 1, sex := 'Males'] %>%
                    .[sex_id == 2, sex := 'Females'] %>%
                    .[age_start >= restrictions_map$min_age[1] & age_start <= restrictions_map$max_age] %>%
                    .[, cf_num := 'cf1']
                  past_cf1[, Y_mean := 1][, Y_mean_lo := 1][, Y_mean_hi := 1]
                }
               
                
                past_dt <- rbind(past_cf3, past_cf2, past_cf1)
                past_dt[, bundle_id := bundle]
                past_dt <- past_dt[!(bundle_id == 198 & sex == 'Females')]
                past_dt <- past_dt[!(bundle_id == 294 & sex == 'Males')]
                
                },
                error = function(cond){
                  message(paste0(bundle, " broke"))
                  return()
                })
  
  # If BPH don't want males
  all_cfs[bundle_id == 294, sex := 'Females']

  
  #model_sum <- model_sum[age_start > restrictions_map$min_age[1] & age_start < restrictions_map$max_age[1] + 5]
  
  #all_cfs <- all_cfs[X_age_start > restrictions_map$min_age[1] & X_age_start < restrictions_map$max_age[1] + 5]
  if(bundle != 4196){
    p <- ggplot(all_cfs[cf_num == 'cf3' | cf_num == 'cf2'], aes(x = X_age_start, y = Y_mean, color = cf_num)) +
      geom_point() + geom_line(aes(linetype = 'GBD 2019')) + 
      #geom_line(data = past_dt[cf_num == 'cf3' | cf_num == 'cf2'], aes(linetype = 'GBD 2017', x = age_start)) +
      #geom_ribbon(aes(ymin = Y_mean_lo   , ymax = Y_mean_hi, fill = cf_num), alpha = 0.6) +
      #geom_ribbon(data = past_dt[cf_num == 'cf3'], aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num, x = age_start), alpha = 0.6) +
      ylab('Correction factor value') + xlab('Age') +
      facet_wrap(~sex, scales = 'free') +
      theme_bw() +
      ggtitle(bun_df[bundle_id == bundle]$bundle_name) +
      scale_color_discrete(name = '') +
      scale_linetype_discrete(name = '') +
      xlim(c(min_age, max_age))
    
    all <- grid.arrange(p, data_plot)
    print(all)
  } else{
    p <- ggplot(all_cfs[cf_num == 'cf3' | cf_num == 'cf2'], aes(x = X_age_start, y = Y_mean, color = cf_num)) +
      geom_point() + geom_line(aes(linetype = 'GBD 2019')) + 
      #geom_line(data = past_dt[cf_num == 'cf3' | cf_num == 'cf2'], aes(linetype = 'GBD 2017', x = age_start)) +
      #geom_ribbon(aes(ymin = Y_mean_lo   , ymax = Y_mean_hi, fill = cf_num), alpha = 0.6) +
      #geom_ribbon(data = past_dt[cf_num == 'cf3'], aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num, x = age_start), alpha = 0.6) +
      ylab('Correction factor value') + xlab('Age') +
      facet_wrap(~sex, scales = 'free') +
      theme_bw() +
      ggtitle(bun_df[bundle_id == bundle]$bundle_name) +
      scale_color_discrete(name = '') +
      scale_linetype_discrete(name = '') +
      xlim(c(min_age, max_age))
    
    all <- grid.arrange(p, data_plot)
    print(all)
  }
 
  # p_ci <- ggplot(all_cfs, aes(x = X_age_start, y = Y_mean, color = cf_num)) +
  #   geom_point() + geom_line(aes(linetype = 'GBD 2019')) + 
  #   geom_line(data = past_dt, aes(linetype = 'GBD 2017', x = age_start)) +
  #   geom_ribbon(aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num), alpha = 0.6) +
  #   geom_ribbon(data = past_dt, aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num, x = age_start), alpha = 0.6) +
  #   ylab('Correction factor value') + xlab('Age') +
  #   facet_wrap(~sex) +
  #   theme_bw() +
  #   ggtitle(bun_df[bundle_id == bundle]$bundle_name)
  # 
  # all_cfs[, max_by_cf := max(Y_mean), by = 'cf_num']
  # maxes <- unique(all_cfs[, c('cf_num', 'max_by_cf')])
  # p_no_3 <- ggplot(all_cfs[cf_num != 'cf3'], aes(x = X_age_start, y = Y_mean, color = cf_num)) +
  #   geom_point() + geom_line(aes(linetype = 'GBD 2019')) +
  #   geom_line(data = past_dt[cf_num != 'cf3'], aes(linetype = 'GBD 2017', x = age_start)) +
  #   #geom_ribbon(aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num), alpha = 0.6) +
  #   #geom_ribbon(data = past_dt, aes(ymin = Y_mean_lo, ymax = Y_mean_hi, fill = cf_num, x = age_start), alpha = 0.6) +
  #   ylab('Correction factor value') + xlab('Age') +
  #   facet_wrap(~sex) +
  #   theme_bw() +
  #   ggtitle('CF2 only')
    #xlim(c(restrictions_map$min_age[1],restrictions_map$max_age[1]))
  #print(p)
  
  
    #xlim(c(restrictions_map$min_age[1],restrictions_map$max_age[1]))
  
  # if(maxes[cf_num == 'cf3']$max_by_cf > 4*maxes[cf_num == 'cf2']$max_by_cf){
  #   grid.arrange(p, p_no_3, data_plot, nrow = 3)
  # } else{
  #   grid.arrange(p, data_plot, nrow = 2)
  # }

})

dev.off()


