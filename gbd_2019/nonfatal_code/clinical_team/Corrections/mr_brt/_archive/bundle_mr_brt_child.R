#############################################################
##Author: Nick Roberts
##Date: 3/4/2019
##Purpose: Launch MR-BRT for every bundle
##Notes: Parallelize your submission of the model
##Updates:
#
###########################################################

# Run the script
print(commandArgs())

bundle <- commandArgs()[5]

data_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/'
write_folder <- '/share/hospital/scratch/marketscan/mr_brt_df_split/bundle/results'

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
                                                "--optimization_method trim_remL \\",
                                                "--trim_pct 0.2")))
  print(command)
  
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

submit_mrbrt(bundle)

