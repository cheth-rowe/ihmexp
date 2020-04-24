# Go through "wild models" and show the number of data points in each mdoel
# Run at different trimming and see the effect on the UI's
# Expect some to be bad, but could alleviate it a bit.
wild_models <- fread('/home/j/temp/nrober75/wild_models.csv')
wild_bundles <- wild_models$index
wild_bun <- loadBundles(wild_bundles)

restrictions <- get_icg_restrictions() %>% setnames(c('yld_age_start', 'yld_age_end'), c('age_start', 'age_end'))
bundle_map <- bundle_icg()

type <- 'bundle'
cf <- 'cf3'
run <- 'wild_models'

for(bundle in bundles){
  preddf <- fread('/share/hospital/clinical_runs/run_1/estimates/corrections/age_sex_preddf.csv')
  bun_df <- fread(paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/prep_bundle/cf3/', bundle, '.csv'))
  rows <- nrow(bun_df)
  
  # Make model
  for(trim in c(0.05,0.1,0.15,0.2)){
    data_folder <- paste0('/share/hospital/clinical_runs/run_1/estimates/corrections/prep_', type, '/', cf, '/')
    write_folder <- paste0('/share/hospital/clinical_runs/', run, '/estimates/corrections/models/', type, '/', cf, '/', trim, '/', bundle, '/')
    dir.create(write_folder, recursive = TRUE)
    
    fit <- run_mr_brt(
      output_dir = write_folder, 
      model_label = 'mr_brt_mod',
      data = paste0(data_folder, bundle, '.csv') ,
      mean_var = y_var,
      se_var = se_var,
      covs = list( 
        #cov_info("zcov1", "Z"),
        cov_info("sex_id", "X"),
        #cov_info("log_haqi", "X"),
        cov_info('age_start', 'X',
                 gprior_mean = 0, 
                 gprior_var = 'inf',
                 bspline_gprior_mean = "0, 0, 0", 
                 bspline_gprior_var = "inf, inf, inf",
                 degree = 3,
                 n_i_knots = 2, 
                 r_linear = TRUE,
                 l_linear = TRUE) ),
      overwrite_previous = TRUE,
      #project = 'proj_hospital',
      method = 'trim_maxL', trim_pct = trim)
    pred <- predict_mr_brt(fit, newdata = preddf, write_draws = FALSE) 
  }
  
  all_preds <- data.table()
  all_input <- data.table()
  for(trim in c(0.05,0.1,0.15,0.2)){
    
    preds <- fread(paste0('/share/hospital/clinical_runs/', run, '/estimates/corrections/models/', type, '/', cf, '/', trim, '/', bundle, '/mr_brt_mod/model_summaries.csv'))
    dir <- preds$working_dir
    preds <- as.data.table(preds)
    
    preds[, trim_pct := trim]
    all_preds <- rbind(all_preds, preds)
    
    train_data <- fread(paste0('/share/hospital/clinical_runs/', run, '/estimates/corrections/models/', type, '/', cf, '/', trim, '/', bundle, '/mr_brt_mod/train_data.csv'))
    train_data[, trim_pct := trim]
    all_input <- rbind(all_input, train_data)
  }
   
  all_input[w <= 0.5, Trimmed := 'Trimmed'][w > 0.5, Trimmed := 'Not trimmed']
  all_input[Trimmed == 'Trimmed', trim_at_pct := paste0(Trimmed, ' at ', trim_pct)]
  all_input[Trimmed == 'Not trimmed', trim_at_pct := 'Not trimmed']
  all_input[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
  
    # Now plot it
    bundle_name <- unique(wild_bun[bundle_id == bundle]$bundle_name)
    train_data <- fread(paste0(dir, 'train_data.csv'))
    train_data[sex_id == 1, sex := 'Males'][sex_id == 2, sex := 'Females']
    pdf('/home/j/temp/nrober75/test_trims.pdf')
    for(t in c(0.05,0.1,0.15,0.2)){
      subset <- all_input[trim_pct == t]
      train_p <- ggplot(all_input, aes(x = age_start, y = exp(log_mean), color = trim_at_pct)) +
        geom_point(aes(shape = sex), size = 2.5) +
        theme_bw() + 
        #facet_wrap(~sex, scales = 'free') +
        scale_color_discrete(name = '') +
        ylab('CF3 input data') +
        xlab('Age') +
        ggtitle(paste0(bundle_name, ' - Trimmed: ', trim, '; N = ',rows))
    }
    train_p <- ggplot(all_input, aes(x = age_start, y = exp(log_mean), color = trim_at_pct)) +
      geom_point(aes(shape = sex), size = 2.5) +
      theme_bw() + 
      #facet_wrap(~sex, scales = 'free') +
      scale_color_discrete(name = '') +
      ylab('CF3 input data') +
      xlab('Age') +
      ggtitle(paste0(bundle_name, ' - Trimmed: ', trim, '; N = ',rows))
    
    preds[X_sex_id == 1, sex := 'Males'][X_sex_id == 2, sex := 'Females']
    results_p <- ggplot(all_preds, aes(x = X_age_start, y = exp(Y_mean))) +
      geom_point() + 
      theme_bw() + 
      ylab('MR-BRT CF3 prediction') +
      xlab('Age') +
      facet_wrap(trim_pct~sex, scales = 'y_free')
    
    results_p_ui <- ggplot(preds, aes(x = X_age_start, y = exp(Y_mean))) +
      geom_point() + 
      theme_bw() + 
      ylab('MR-BRT CF3 prediction') +
      xlab('Age') +
      facet_wrap(~sex, scales = 'free') +
      geom_ribbon(aes(ymax = exp(Y_mean_hi), ymin = exp(Y_mean_lo)), alpha = 0.3)
    
    pdf(paste0('/home/j/temp/nrober75/wild_models_run/', bundle, '_', trim, '.pdf'))
    all <- grid.arrange(train_p, results_p, results_p_ui, nrow = 3)
    dev.off()
  }
  
  

}
