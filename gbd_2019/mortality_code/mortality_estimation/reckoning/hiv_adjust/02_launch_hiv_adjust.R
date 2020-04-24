
## Run Ensemble process

## Outputs: HIV adjusted death numbers



###############################################################################################################
## Set up settings
  rm(list=ls())
  library(foreign); library(RMySQL); library(dplyr); library(data.table); library(rhdf5); library(assertable)
  
  if (Sys.info()[1]=="Windows") {
    root <- "J:/" 
    user <- Sys.getenv("USERNAME")
  } else {
    root <- "/home/j/"
    user <- Sys.getenv("USER")
  }
  
## Grab functions
  library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
  library(mortcore, lib = "/home/j/WORK/02_mortality/shared/r")
  source(paste0("/share/cc_resources/libraries/current/r/get_population.R"))

  singularity_shell <- paste0("/homes/", user, "/hiv_gbd2019/singR_shell.sh")


###############################################################################################################
## Set start and end options
## Note: These are categorized toggles rather than numeric because of the number of possibilities that can be run in parallel
## First, do you want to run the ensemble process?
  run_ensemble <- T
  finish_reckoning <- F
  
## Run 01b prep script
  run_prep <- F

## Set other run options
  spec_name <- "190630_rhino_combined"# Update this with the new Spectrum run results each time
  test <- F # Test submission of everything 
  file_del <- F
  run_comment <- paste0("HIV Decomp Step 4 ", Sys.Date())
  
## Set years that you are running for
  start_year <- 1950
  gbd_year <- 2019

  gbd_round <- ifelse(gbd_year == 2019, 6, gbd_year - 2012)
  years <- c(start_year:gbd_year) 

###############################################################################################################
## Get version number, create process lineages, make directories
 # new_hiv_adjust_version <- gen_new_version(model_name="hiv adjustment", model_type="estimate",  comment=run_comment)
  
  #print(paste0("Creating new version ", new_hiv_adjust_version))
  
  #new_hiv_adjust_version <- 205 #Original decomp 4
  new_hiv_adjust_version <- 277 #Final, Final GBD 19 version
  #new_hiv_adjust_version <- 278 #Test version using ST-GPR final for group 2B
  
  # Get parent processes of the version, and update version mapping appropriately
  run_parents <- get_best_versions(c("mlt life table estimate", "mlt death number estimate","population estimate"))

  #mlt_lt_version <- run_parents[["mlt life table estimate"]]  
  mlt_lt_version <- 240
  
  ## TODO: This call is throwing an error
  #gen_parent_child(child_process = "hiv adjustment estimate", child_id = new_hiv_adjust_version, parent_runs = run_parents)  

  ## Set code_dir
  ## TODO: Keep this in /ihme/code
  code_dir <- paste0("/homes/",user,"/hiv_gbd2019/04_prep_all_results/")
  setwd(code_dir)
  
  # Tag Git commit with the current version number
  # git_tag(code_dir=code_dir, name=paste0("v",new_hiv_adjust_version), message=run_comment)
  
## Set directories
  #new_hiv_adjust_version <- "test_84"
  master_dir <- paste0("/ihme/mortality/hiv_adjust/",new_hiv_adjust_version)
  input_dir <- paste0(master_dir, "/inputs")
  out_dir <- paste0(master_dir,"/envelope_whiv")
  out_dir_hiv_free <- paste0(master_dir,"/envelope_hivdel")
  out_dir_hiv <- paste0(master_dir,"/hiv_output")
  out_dir_products <- paste0(root,"/WORK/02_mortality/03_models/5_lifetables/results/post_reckoning_products")
  out_dir_adjust_summary <- paste0(master_dir,"/results")
  
  ## Create folder tree for reckoning runs (new version = new tree)
  for(dir in c(master_dir, input_dir, out_dir, out_dir_hiv_free,
               out_dir_adjust_summary, out_dir_hiv)
  ) {
    system(paste0("mkdir ", dir))
  }

## Grab locations
  codes <- data.table(get_locations(level = "lowest", gbd_year = gbd_year))
  run_countries <- unique(codes$ihme_loc_id)
  print(run_countries)
  write.csv(codes,paste0(input_dir,"/locations_lowest.csv"), row.names = F)

  hiv_groups <- data.table(get_locations(hiv_metadata = T, gbd_year = gbd_year))
  mort_locs <- as.data.table(hiv_groups)[, list(location_id, location_name, ihme_loc_id)]
  write.csv(mort_locs, paste0(input_dir, "/locations_mort.csv"), row.names = F)

  all_locs <- data.table(get_locations(level = "all", gbd_year = gbd_year))
  write.csv(all_locs, paste0(input_dir, "/locations_all.csv"))
  
  
## Save spectrum version in a flat file, until HIV is included in DB tracking for Mortality DB
  hiv_output <- data.table(used_spec_run = spec_name)
  write.csv(hiv_output,paste0(out_dir_hiv,"/used_hiv_run.csv"),row.names=F)
  
## Create draw maps to scramble draws so that they are not correlated over time
  create_draws <- function(location_id) {
    new_seed <- location_id * 100 + 121 # Just to avoid any unintended correlation with other seeds that rely on location_id
    set.seed(new_seed)
    data <- data.table(expand.grid(location_id=location_id,old_draw=c(0:999)))
    data[,draw_sort:=rnorm(1000)]
    data <- data[order(location_id,draw_sort),]
    data[,new_draw:=seq_along(draw_sort)-1,by = location_id] # Creates a new variable with the ordering based on the values of draw_sort 
    data[,draw_sort:=NULL]
  }
  draw_map <- rbindlist(lapply(unique(codes$location_id),create_draws))
  
  write.csv(draw_map,paste0(input_dir, "/draw_map.csv"), row.names = F)

  
###############################################################################################################
## Create HIV type file
hiv_groups <- as.data.table(get_locations(hiv_metadata = T, gbd_year = gbd_year))

## Create groups list
master_types <- hiv_groups[,list(location_id, ihme_loc_id, group)]

write.csv(master_types,paste0(input_dir,"/ensemble_groups.csv"), row.names = F)

## Create population file
age_map_all <- data.table(get_age_map(type = "all"))
age_map_all <- age_map_all[(age_group_id != 1 & age_group_id <= 20) | age_group_id == 28 | (age_group_id >=30 & age_group_id <=32) | age_group_id == 235, list(age_group_name_short,age_group_id)]
write.csv(age_map_all, paste0(input_dir,"/age_map_01.csv"), row.names = F)

age_map_lt <- data.table(get_age_map(type="lifetable"))
age_map_lt <- age_map_lt[,list(age_group_id,age_group_name_short)]
setnames(age_map_lt,"age_group_name_short","age")
write.csv(age_map_lt, paste0(input_dir,"/age_map_02.csv"), row.names = F)

age_map_mort <- data.table(get_age_map())
age_map_mort <- age_map_mort[,list(age_group_id,age_group_name)]
write.csv(age_map_mort, paste0(input_dir,"/age_map_mort.csv"), row.names = F)

pop_age_ids <- unique(c(age_map_all$age_group_id, age_map_lt$age_group_id, age_map_mort$age_group_id), 2:4)
pop_age_ids <- pop_age_ids[(pop_age_ids != 1 & pop_age_ids <= 20) | pop_age_ids == 22 | pop_age_ids == 28 | (pop_age_ids >=30 & pop_age_ids <=32) | pop_age_ids == 235]
pop <- get_population( location_id = -1, location_set_id = 21, year_id = years, sex_id = c(1:3), age_group_id = pop_age_ids, gbd_round_id = gbd_round, decomp_step="step4")
pop <- pop[location_id %in% unique(hiv_groups$location_id), list(age_group_id, location_id, year_id, sex_id, population)]

#Assert population values
assert_values(pop, "population", "gt", 0)
assert_values(pop, "population", "not_na")
id_vars <- list(year_id = years, age_group_id = pop_age_ids, sex_id = c(1:3), location_id = unique(hiv_groups$location_id))
assertable::assert_ids(pop, id_vars)

write.csv(pop,paste0(input_dir, "/pop.csv"), row.names = F)

###############################################################################################################
## Remove Existing Files
if (file_del == T) {
  if (run_ensemble==T) {
    print("Deleting 01_ensemble output -- should take 5-10 minutes (maybe more)")
    system(paste0("perl -e 'unlink <",out_dir,"/en_summ_*dta>' ")) # Not produced by my code, but should be deleted to not be out of date
    system(paste0("perl -e 'unlink <",out_dir,"/scalars_*.dta>' "))
    system(paste0("perl -e 'unlink <",out_dir,"/env_*.h5>' "))
    system(paste0("perl -e 'unlink <",out_dir,"/comparison*.csv>' "))
    system(paste0("perl -e 'unlink <",out_dir_hiv_free,"/env_*.h5>' "))
    system(paste0("perl -e 'unlink <",out_dir_hiv,"/hiv_death_*.csv>' "))
    system(paste0("perl -e 'unlink <",out_dir_hiv,"/reckon_reporting_*.csv>' "))
    system(paste0("perl -e 'unlink <",out_dir_adjust_summary,"/results*.csv>' "))
    system(paste0("perl -e 'unlink <",out_dir_adjust_summary,"/result/results_hiv_adjust_v",new_hiv_adjust_version,".csv>' "))
  }
 }

###############################################################################################################
## Submit Jobs
  
## NOTE: MAKE SURE THAT IN HIV, PREP_SPEC_RESULTS HAS BEEN RUN ON ALL INPUT DATA

#These CHINA locs are misclassified in the database as 2A
special.china <- c("CHN_493","CHN_495","CHN_521")



## Run ensemble process
   if (run_ensemble==T) {  
    for(country in run_countries) {
           group <- paste0(unique(master_types$group[master_types$ihme_loc_id==country]))
           
           if(country %in% special.china){
             group <- "2B"
           }
           
          ensemble.string <- paste0("qsub -l m_mem_free=50.0G -l fthread=2 -l h_rt=03:00:00 -l archive=True -q all.q ",
                      "-e /share/temp/sgeoutput/",user,"/errors/ ",
                      "-o /share/temp/sgeoutput/", user, "/output ",
                      "-N ensemble_01_",country," ",
                      "-cwd -P proj_hiv ",
                        singularity_shell," ",
                        code_dir,"/run_ensemble.R ",
                        country," ", group," ",
                        spec_name," ", new_hiv_adjust_version," ", mlt_lt_version, " ",start_year," ", gbd_year)
          print(ensemble.string)      
          system(ensemble.string)
          Sys.sleep(0.2)
       
           
        }
     
    

    print("Waiting 20 minutes before checking results for 01_ensemble")
    Sys.sleep(60*20)
    

    ## Check for output
    check_loc_results(run_countries,out_dir, prefix="env_",postfix=".h5")
    check_loc_results(run_countries,out_dir_hiv_free, prefix="env_",postfix=".h5")
    check_loc_results(run_countries,out_dir_hiv, prefix="hiv_death_",postfix=".csv")
   }    


