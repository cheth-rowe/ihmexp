#############################################################
##Author: Nick Roberts
##Date: 3/8/2019
##Purpose: Run lmer mixed-effects GBD 2017 models on each bundle with USA states/years
##Notes: Probably parallelize and use qsub if it's taking a long time
##Updates:
#
###########################################################

library(data.table)
library(tibble)
library(lme4)
library(RMySQL)

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
# Multiplied it by 1000 draws to merge on
copy_draws <- function(data){
  append_draws <- list()
  for(draws in 0:999){
    append_draws[[paste0(draws)]] <- copy(data)
    append_draws[[paste0(draws)]]$draw <- draws
  }
  append_draws <- rbindlist(append_draws)
  return(append_draws)
}


data_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/'
spec_file <- paste0('no_inf/', bundle,'.csv')

df <- fread(paste0(data_folder, spec_file))

all_bundles <- fread('/share/hospital/scratch/marketscan/mr_brt_df_split/marketscan_prepped.csv')


#### Format and make model ####
bundles <- unique(all_bundles$bundle_id)
bun_df <- loadBundles(bundles)

mclapply(bundles, function(bundle){
  print(bundle)
  pdf(paste0('/home/j/temp/hospital/2016/data/mr_brt/bundle_runs/lmer_bundles/bootmer/random_effects_lmer_bundle_results_', bundle, '.pdf'), onefile = TRUE, compress = TRUE)
  df <- all_bundles[bundle_id == bundle][!is.na(cf3) & cf3 != Inf & cf3 != 0]
  df[, log_mean := log(cf3)]
  df[, log_se := log(se)]
  
  df <- df[bundle_id == bundle, c('log_mean', 'log_se', 'age_start', 'sex_id', 'haqi', 'ip_env', 'op_env', 'location_id', 'year_start')]
  
  mod <- lmer(log_mean ~ age_start + sex_id + haqi + ip_env + op_env + (1|location_id), data = df)
  
  #### Get draws, merge on to dataset ####
  betas <- fixef(mod)
  beta_names <- names(betas)
  cov_matrix <- as.matrix(vcov(mod))
  
  mvrnorm <- as.data.table(MASS::mvrnorm(n = 1000, mu = betas, cov_matrix))
  setnames(mvrnorm, '(Intercept)', 'intercept')
  colnames(mvrnorm) <- paste(colnames(mvrnorm), "beta", sep = "_")
  mvrnorm[, draw := .I-1]
  
  # Use simulate
  # mod_sim <- simulate(mod, 1000, newdata = df) %>% as.data.table()
  # mod_sim[, id := .I]
  # 
  # df <- cbind(df, mod_sim)
  # df <- melt(df, id.vars = c(names(df)[1:9], 'id'),  value.var = patterns('sim_'))
  # 
  # Use bootmer too
  test_draws <- bootMer(mod, FUN = predict, nsim = 1000)
  test_df <- cbind(df, t(test_draws$t))
  test_df[, id := .I]
  test_df[, pred := predict(mod, newdata = test_df)]
  
  test_df_long <- melt(test_df, id.vars = c(names(test_df)[1:9], 'id', 'pred'), value.var = patterns('V'))
  test_df_long[, draw_lo := quantile(value, 0.025), by = 'id']
  test_df_long[, draw_hi := quantile(value, 0.975), by = 'id']
  test_df_long[, min_lo := min(draw_lo), by = c('age_start', 'sex_id')]
  test_df_long[, max_hi := max(draw_hi), by = c('age_start', 'sex_id')]
  plot_df <- unique(test_df_long[, c('age_start', 'sex_id', 'pred', 'min_lo', 'max_hi')])
  
  bun_name <- bun_df[bundle_id == bundle]$bundle_name
  p1 <- ggplot(plot_df, aes(x = age_start, y = exp(pred))) + 
    geom_ribbon(aes(ymin = exp(min_lo), ymax = exp(max_hi)), alpha = 0.3, fill = 'dark green') + 
    geom_point(aes(color = 'Predictions'), color = 'dark green') + 
    facet_wrap(~sex_id) + 
    theme_bw() +
    ylab('Random effects prediction') + xlab('Age') + 
    ggtitle(paste0(bun_name, ' predictions on USA training data'))
  
  coefs <- summary(mod)
  coefs <- coefs$coefficients
  coefs <- as.data.table(coefs)
  coefs[, xcov := c('intercept', 'op_env', 'ip_env', 'haqi', 'sex_id', 'age_start')]
  coefs <- coefs[, c('xcov', 'Estimate', 'Std. Error', 't value')]
  
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  tbl <- tableGrob(coefs, rows = NULL, theme = tt)
  
  all <- grid.arrange(p1, tbl,
                      nrow = 2, as.table = TRUE, heights = c(3,1))
  print(all)
  dev.off()
  
  #### This was original code without predicting draws ofrandom effects
  # Make actual prediction
  # df[, pred := predict(mod, newdata = df)]
  # df[, draw_lo := quantile(value, 0.025), by = 'id']
  # df[, draw_hi := quantile(value, 0.975), by = 'id']
  # 
  # df[, min_lo := min(draw_lo), by = c('age_start', 'sex_id')]
  # df[, max_hi := max(draw_hi), by = c('age_start', 'sex_id')]
  # plot_df <- unique(df[, c('age_start', 'sex_id', 'pred', 'min_lo', 'max_hi')])
  # bun_name <- bun_df[bundle_id == bundle]$bundle_name
  # 
  # p1 <- ggplot(plot_df, aes(x = age_start, y = exp(pred))) + 
  #   geom_ribbon(aes(ymin = exp(min_lo), ymax = exp(max_hi)), alpha = 0.3, fill = 'dark green') + 
  #   geom_point(aes(color = 'Predictions'), color = 'dark green') + 
  #   facet_wrap(~sex_id) + 
  #   theme_bw() +
  #   ylab('Random effects prediction') + xlab('Age') + 
  #   ggtitle(paste0(bun_name, ' predictions on USA training data'))
  #   #scale_y_continuous(limits=c(0,150))
  # 
  # df <- copy_draws(df)
  # df <- merge(df, mvrnorm, by = 'draw')
  
  # ## Get draws of random effects
  # re = ranef(mod, condVar = TRUE)$location_id
  # 
  # fe_d
  # re_draws = NULL
  # for(i in seq(dim(attr(re, 'postVar'))[3])){
  #   re_draws = cbind(re_draws, mvrnorm(n = 1000, mu = as.matrix(re[i,]), 
  #                                      Sigma = attr(re, 'postVar')[,,i]))
  # }
  # 
  
  # Manually make predictions. Do you need to add random effects? Do I need to get draws from random effects
  # df[, draw_pred := intercept_beta + age_start*age_start_beta + sex_id*sex_id_beta + haqi*haqi_beta + ip_env*ip_env_beta + op_env*op_env_beta]
  # 
  # df[, mean_draw := mean(draw_pred), by = c('age_start', 'sex_id', 'location_id', 'year_start')]
  # df[, low_draw := quantile(draw_pred, probs = 0.025), by = c('age_start', 'sex_id', 'location_id', 'year_start')]
  # df[, hi_draw := quantile(draw_pred, probs = 0.975), by = c('age_start', 'sex_id', 'location_id', 'year_start')]
  # df[, min_lo := min(low_draw), by = c('age_start', 'sex_id')]
  # df[, max_hi := max(hi_draw), by = c('age_start', 'sex_id')]
  # 
  # plot_df <- unique(df[, c('age_start', 'sex_id', 'pred', 'mean_draw', 'min_lo', 'max_hi')])
  # p1 <- ggplot(plot_df, aes(x = age_start, y = exp(pred))) + 
  #   geom_ribbon(aes(ymin = exp(min_lo), ymax = exp(max_hi)), alpha = 0.3, fill = 'dark green') + 
  #   geom_point(aes(color = 'Predictions'), color = 'dark green') + 
  #   facet_wrap(~sex_id) + 
  #   theme_bw() +
  #   ylab('Random effects prediction') + xlab('Age') + 
  #   ggtitle(paste0(bun_name, ' predictions on USA training data'))
  # 
  # 
  # coefs <- summary(mod)
  # coefs <- coefs$coefficients
  # coefs <- as.data.table(coefs)
  # coefs[, xcov := c('intercept', 'op_env', 'ip_env', 'haqi', 'sex_id', 'age_start')]
  # coefs <- coefs[, c('xcov', 'Estimate', 'Std. Error', 't value')]
  # 
  # tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  # tbl <- tableGrob(coefs, rows = NULL, theme = tt)
  # 
  # all <- grid.arrange(p1, tbl,
  #                     nrow = 2, as.table = TRUE, heights = c(3,1))
  # print(all)
  # dev.off()
}, mc.cores = 30)
  




# ## And make the same predictions as you did with MR-BRT
# 
# intercept <- 1
# op_env <- seq(1,13,0.5)
# ip_env <- seq(0,1,0.05)
# age_start <- c(seq(0,95,5), 1)
# haqi <- seq(60,90,1)
# sex_id <- c(1,2)
# test <- tidyr::crossing(intercept, op_env, ip_env, haqi, sex_id) %>% as.data.table()
# 
# test[, pred := predict(mod, newdata = test)]
# 
# 
# 
# er <- ranef(mod)
# er$location_id[df$location_id, "(Intercept)"]
# er$Beach[dataF$Beach, "(Intercept)"]