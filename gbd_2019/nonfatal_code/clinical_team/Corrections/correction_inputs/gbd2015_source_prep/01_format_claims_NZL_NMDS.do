// *********************************************************************************************************************************************************************
// *********************************************************************************************************************************************************************
// Author: 		
// Date: 			September 2015
// Modified:		--
// Project:		GBD
// Purpose:		Template for formatting NZL_NMDS for the combined hospital database


** **************************************************************************
** CONFIGURATION
** **************************************************************************
	** ****************************************************************
	** Prepare STATA for use
	**
	** This section sets the application perferences.  The local applications
	**	preferences include memory allocation, variables limits, color scheme,
	**	 and setting a local for the date.
	**
	** ****************************************************************
		// Set application preferences
			// Clear memory and set memory and variable limits
				clear all
				set mem 10G
				

			// Set to run all selected code without pausing
				set more off

			// Set graph output color scheme
				set scheme s1color

			// Get date
				local today = date(c(current_date), "DMY")
				local year = year(`today')
				local month = string(month(`today'),"%02.0f")
				local day = string(day(`today'),"%02.0f")
				local today = "`year'_`month'_`day'"


		** ****************************************************************
		** SET LOCALS
		**
		** Set data_name local and create associated folder structure for
		**	formatting prep.
		**
		** ****************************************************************
			// Data Source Name
				local data_name "NZL_NMDS"
			// Original data folder
				local input_folder "FILEPATH
			// Log folder
				local log_folder FILEPATH
				capture mkdir "`log_folder'"
			// Output folder
				local output_folder FILEPATH
				local archive_folder "`output_folder'/_archive"
				capture mkdir "`output_folder'"
				capture mkdir "`archive_folder'"


		** ****************************************************************
		** CREATE LOG
		** ****************************************************************
			capture log close
			log using "`log_folder'/00_format_`today'.log", replace

	// IHME written mata functions which speed up collapse and egen
	//do "FILEPATH" 

** **************************************************************************
** RUN FORMAT PROGRAGM
** **************************************************************************
	// GET DATA
		clear
		gen y = .
		forvalues y = 2000/2014 {
			append using "`input_folder'/NZL_NMDS_PUS8929_EVENTS_`y'_Y2015M10D08.DTA"
			replace y = `y' if y == .
		}
		// now add 2015 which follows a different naming format
		local year = 2015
		foreach month in JAN_JUN JUL_DEC {
			append using "`input_folder'/NZL_NMDS_2015_EVENTS_`month'_Y2017M04D21.DTA"
			replace y = `year' if y == .
		}
		
		tempfile EVENTS
		save `EVENTS', replace
		
		clear
		forvalues y = 2000/2014 {
			append using "`input_folder'/NZL_NMDS_PUS8929_DIAGS_`y'_Y2015M10D08.DTA"
		}
		// add 2015 data for diagnoses
		foreach month in JAN_JUN JUL_DEC {
			append using "`input_folder'/NZL_NMDS_2015_DIAGS_`month'_Y2017M04D21.DTA"
		}
		
		drop OP_ACDTE 	// Operation/Procedure date
		drop if DIAG_TYP == "O"	// Operation/Procedures
		
		// icd_vers (string): ICD version - "ICD10", "ICD9_detail"
			gen icd_vers = "ICD10" if inlist(CLIN_SYS, "10", "11", "12", "13", "14")
			replace icd_vers = "ICD9_detail" if CLIN_SYS == "01" | CLIN_SYS == "02" | CLIN_SYS == "06"
			drop CLIN_SYS
		
		tempfile DIAGS
		save `DIAGS', replace
		
		// DIAGNOSIS & ECODE VARIABLES
			// Primary dx
				keep if DIAG_TYP == "A"
				rename CLIN_CD dx_
				isid EVENT_ID
				replace DIAG_SEQ = 1
				drop DIAG_TYP
				reshape wide dx_, i(EVENT_ID icd_vers) j(DIAG_SEQ)
				sort EVENT_ID
				tempfile primary
				save `primary', replace
				
			// Secondary dx's. Only keep first 15 diagnoses
				use `DIAGS', clear
				keep if DIAG_TYP == "B"
				rename CLIN_CD dx_
				sort EVENT_ID DIAG_SEQ
				by EVENT_ID : gen num = _n
				drop DIAG_SEQ DIAG_TYP
				drop if num > 15
				reshape wide dx_, i(EVENT_ID icd_vers) j(num)
				sort EVENT_ID
				tempfile secondary
				save `secondary', replace
			
			// Ecodes
				use `DIAGS', clear
				keep if DIAG_TYP == "E"
				rename CLIN_CD ecode_
				sort EVENT_ID DIAG_SEQ
				by EVENT_ID : gen num = _n
				drop DIAG_SEQ DIAG_TYP
				reshape wide ecode_, i(EVENT_ID icd_vers) j(num)
				sort EVENT_ID
				
				merge 1:1 EVENT_ID using `primary', nogen
				merge 1:1 EVENT_ID using `secondary', nogen
				
		
		merge 1:1 EVENT_ID using `EVENTS', keep(3) nogen
		
		

	// ENSURE ALL VARIABLES ARE PRESENT
		// source (string): source name
			gen source = "`data_name'"
		// year (numeric)
			gen year = string(evendate, "%tdCCYY")
			destring year, replace
			drop evendate
		// NID (numeric)
			gen NID = .
			replace NID = 220786 if year == 2000
			replace NID = 220787 if year == 2001
			replace NID = 220788 if year == 2002
			replace NID = 220789 if year == 2003
			replace NID = 220790 if year == 2004
			replace NID = 220791 if year == 2005
			replace NID = 220792 if year == 2006
			replace NID = 220793 if year == 2007
			replace NID = 220794 if year == 2008
			replace NID = 220795 if year == 2009
			replace NID = 220796 if year == 2010
			replace NID = 220797 if year == 2011
			replace NID = 220798 if year == 2012
			replace NID = 220799 if year == 2013
			replace NID = 220800 if year == 2014
			replace NID = 293984 if year == 2015
			drop y
		// iso3 (string)
			gen iso3 = "NZL"
		// subdiv (string)
			gen subdiv = ""
		// location_id (numeric)
		// set location_id based on ethnicgp
		// See NZL govt's ethnic codes in the link below.  The first digit (level 1)
		// denotes the larger group and the second digit denotes more detail (level 2).
		// http://www.health.govt.nz/nz-health-statistics/data-references/code-tables/common-code-tables/ethnicity-code-tables
			gen location_id = .
			replace ethnicgp = substr(ethnicgp, 1, 1)  // keep first char (of a string)
			replace location_id = 44850 if ethnicgp == "2" // 2 for Maori
			replace location_id = 44851 if ethnicgp != "2"  // if not 2 then nor Maori
		// national (numeric): 0 = no, 1 = yes
			gen national = 1
		// admission date
			rename evstdate adm_date
		// age (numeric)
			rename AGE_AT_DISCHARGE age
			// gen age = .
			// replace age = 0 if AGE_AT_DISCHARGE == 0
			// replace age = 1 if AGE_AT_DISCHARGE >=1 & AGE_AT_DISCHARGE <=4
			// forvalues a = 5(5)90 {
			// 	replace age = `a' if AGE_AT_DISCHARGE >= `a' & AGE_AT_DISCHARGE <= (`a' + 4)
			// }
			// replace age = 95 if AGE_AT_DISCHARGE >= 95
			
			// this method was not correct
			// gen age = .
			// replace age = 0 if AGE_AT_DISCHARGE == 0
			// replace age = floor(AGE_AT_DISCHARGE / 5) * 5
			// replace age = 95 if age > 95
			// drop AGE_AT_DISCHARGE
		// frmat (numeric): find the WHO format here "Age formats documentation.xlsx"
			gen frmat = .
		// im_frmat (numeric): from the same file as above
			gen im_frmat = .
		// sex (numeric): 1=male 2=female 9=missing
			gen sex = .
			replace sex = 1 if gender == "M"
			replace sex = 2 if gender == "F"
			replace sex = 9 if sex == .
		// platform (string): "Inpatient", "Outpatient"
			gen platform = "Inpatient"
		// patient_id (string)
			// Var already present as patient_id
		// icd_vers (string): ICD version - "ICD10", "ICD9_detail"
			
		// dx_* (string): diagnoses
			
		// ecode_* (string): variable if E codes are specifically mentioned
			
		// Inpatient variables
			// discharges (numeric)
				gen metric_discharges = 1
			// day cases (numeric)
				gen metric_day_cases = 1 if real(los) == 0
			// bed_days (numeric)
				gen metric_bed_days = real(los)
			// deaths (numeric)
				gen metric_deaths = 0
				replace metric_deaths = 1 if END_TYPE == "DD"
		// Keep only inpatient 24 hour stays ie drop day cases
			keep if real(los) != 0

	// VARIABLE CHECK
		// If any of the variables in our template are missing, create them now (even if they are empty)
		// All of the following variables should be present
			#delimit;
			order
			adm_date
			iso3 subdiv location_id national
			source NID
			year
			age frmat im_frmat
			sex platform patient_id
			icd_vers dx_* ecode_*
			metric_*;
		// Drop any variables not in our template of variables to keep
			keep
			adm_date
			iso3 subdiv location_id national
			source NID
			year
			age frmat im_frmat
			sex platform patient_id
			icd_vers dx_* ecode_*
			metric_*;
			#delimit cr

	// COLLAPSE DATA
		//collapse (sum) metric_*, by(iso3 subdiv location_id national source NID year age frmat im_frmat sex platform patient_id icd_vers dx_* ecode_*) fast

** **************************************************************************
** RUN SELECT PRIMARY PROGRAGM
** **************************************************************************
// drop bed days of zero
		keep if metric_bed_days > 0
	// Select primary GBD code
	** NOTE: use the primary GBD code unless there is an E code (per  2013_12_23)
		// gen yld_cause_primary = dx_1_yld_cause
		// gen cause_primary = dx_1
		// // Prioritize External Cause of Injury & Service Codes over Nature of Injury
		// 	// Diagnosis codes
		// 		local dx_count = 0
		// 		foreach var of varlist dx_* {
		// 			local dx_count = `dx_count' + 1
		// 		}
		// 		local dx_count = `dx_count'
		// 		// Prioritizing ICD-10 V, W, X, Y codes and ICD-9 E codes
		// 			forvalues dx = `dx_count'(-1)1 {
		// 				display "Prioritizing ICD-10 V, W, X, Y codes and ICD-9 E codes in dx_`dx'"
		// 				//replace yld_cause_primary = dx_`dx'_yld_cause if (inlist(substr(dx_`dx'_mapped,1,1),"V","W","X","Y") & icd_vers == "ICD10") | (inlist(substr(dx_`dx'_mapped,1,1),"E") & icd_vers == "ICD9_detail")
		// 				replace cause_primary = dx_`dx' if (inlist(substr(dx_`dx',1,1),"V","W","X","Y") & icd_vers == "ICD10") | (inlist(substr(dx_`dx',1,1),"E") & icd_vers == "ICD9_detail")
		// 			}
		// 	// External cause set
		// 		local ecode_count = 0
		// 		foreach var of varlist ecode_* {
		// 			local ecode_count = `ecode_count' + 1
		// 		}
		// 		local ecode_count = `ecode_count'
		// 		// Prioritizing ICD-10 V, W, X, Y codes and ICD-9 E codes
		// 			forvalues dx = `ecode_count'(-1)1 {
		// 				display "Prioritizing ICD-10 V, W, X, Y codes and ICD-9 E codes in ecode_`dx'"
		// 				//replace yld_cause_primary = ecode_`dx'_yld_cause if (inlist(substr(ecode_`dx'_mapped,1,1),"V","W","X","Y") & icd_vers == "ICD10") | (inlist(substr(ecode_`dx'_mapped,1,1),"E") & icd_vers == "ICD9_detail")
		// 				replace cause_primary = ecode_`dx' if (inlist(substr(ecode_`dx',1,1),"V","W","X","Y") & icd_vers == "ICD10") | (inlist(substr(ecode_`dx',1,1),"E") & icd_vers == "ICD9_detail")
		// 			}
	// Reshape long
		//drop dx_*_original dx_*_yld_cause ecode_*_yld_cause ecode_*_original
		// drop dx_1
		// rename cause_primary dx_1

		// drop Ecodes since we've already replaced dx_1 with an ecode where appropriate
		drop ecode_*
		//rename dx_*_mapped dx_mapped_*
		//rename ecode_*_mapped ecode_mapped_*
		//reshape long dx_mapped_ ecode_mapped_, i(iso3 subdiv location_id national source NID year age frmat im_frmat sex platform patient_id icd_vers yld_cause_primary cause_primary) j(dx_num)
		// stata kept saying we specified too many vars for reshape or that smaller # weren't unique. so we explicitly created an index
		. gen long id = _n
		reshape long dx_, i(id) j(dx_ecode_id)
		replace dx_ecode_id = 2 if dx_ecode_id > 2
		drop if missing(dx_)
		//drop id

	// Collapse
		//keep iso3 subdiv location_id national source NID year age frmat im_frmat sex platform patient_id icd_vers dx_ecode_id metric_* dx_
		//collapse (sum) metric_*, by(iso3 subdiv location_id national source NID year age frmat im_frmat sex platform patient_id icd_vers dx_ecode_id dx_) fast

		//rename dx_ dx_mapped_
	
** **************************************************************************
** RUN EPI COMPILE CLUSTER PROGRAGM
** **************************************************************************

	// STANDARDIZE PLATFORMS
		replace platform = "Outpatient" if platform == "Emergency"
		replace platform = "1" if platform == "Inpatient" | platform == "Inpatient 1"
		drop if platform == "Inpatient 2"
		replace platform = "2" if platform == "Outpatient"
		destring platform, replace
		assert platform != .
		
	// STANDARDIZE METRICS
	// 	gen cases = .
	// 		capture confirm var metric_discharges
	// 		if !_rc {
	// 			replace cases = metric_discharges if platform == 1
	// 		}
	// 		capture confirm var metric_discharges_weighted
	// 		if !_rc {
	// 			replace cases = metric_discharges_weighted if platform == 1
	// 		}
	// 		capture confirm var metric_visits
	// 		if !_rc {
	// 			replace cases = metric_visits if platform == 2
	// 		}
	// 		capture confirm var metric_visits_weighted
	// 		if !_rc {
	// 			replace cases = metric_visits_weighted if platform == 2
	// 		}
	// 	gen deaths = .
	// 		capture confirm var metric_deaths
	// 		if !_rc {
	// 			replace deaths = metric_deaths
	// 		}
	// 		capture confirm var metric_deaths_weighted
	// 		if !_rc {
	// 			replace deaths = metric_deaths_weighted
	// 		}

	// // STANDARDIZE AGES
	// 	drop if age == .
	// 	replace age = 95 if age > 95
	// 	replace age = 0 if age < 1

	// 	rename age age_start
	// 	gen age_end = age_start + 4
	// 	replace age_end = 1 if age_start == 0
	// 	replace age_end = 99 if age_start == 95
	// 	replace age_end = 4 if age_start == 1

	// STANDARDIZE SEX
		//drop if sex != 1 & sex != 2

	// LOCATION_ID
		// rename location_id location_id_orig
		// merge m:1 iso3 using "location_ids.dta"
		// 	assert _m != 1
		// 	drop if _m == 2
		// 	drop _m
		// replace location_id = location_id_orig if location_id_orig != .
		// drop location_id_orig

	// COLLAPSE CASES
		// fastcollapse cases deaths, type(sum) by(iso3 subdiv location_id source year age_* sex platform icd_vers cause_primary)
		// fastcollapse discharges deaths, type(sum) by(iso3 subdiv location_id source national year age_* sex platform icd_vers dx_ dx_ecode_id NID)
		//collapse (sum) cases deaths, by(iso3 subdiv location_id source national year age_* sex platform icd_vers dx_ dx_ecode_id NID) fast

	// modify columns to fit our structure
		rename dx_ cause_code
		rename dx_ecode_id diagnosis_id

	// SAVE
		compress
		save "`output_folder'/NZL_NMDS.dta", replace
		save "`archive_folder'/NZL_NMDS_`today'.dta", replace

	capture log close


// *********************************************************************************************************************************************************************
// *********************************************************************************************************************************************************************
