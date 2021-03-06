# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 15:52:35 2017

Subsets incidence to data-dense geographies only
Calculates combined standard error given two input means + SEs

This program must be run in the GBD Computing Environment

@author: USERNAME
updated: Feb 08 2019
"""

import numpy as np
from os.path import join
import pandas as pd
from db_queries import get_location_metadata, get_model_results, get_population
from db_tools.ezfuncs import query
#import new_emr

# Subset a df to all results in a particular category, and all results NOT
#   in that category
def split_on_cat(in_df, column, category):
    all_in_cat = in_df.loc[in_df[column] == category]
    not_in_cat = in_df.loc[in_df[column] != category]
    return(all_in_cat, not_in_cat)

# Find the new standard error generated by dividing two distributions
def divided_se(num_mean, num_se, denom_mean, denom_se):
    new_se = 0
    if num_mean != 0:
        # Formula for calculating the combined SE
        new_se = np.sqrt((num_mean**2 / denom_mean**2) *
                         ((num_se**2 / num_mean**2) +
                         (denom_se**2 / denom_mean**2)))
    return new_se


# Takes a df, splits on incidence, subsets incidence to ONLY data dense
#   countries, concatenates this with prevalence for all countries
def subset_incidence_only(in_df,split_col = 'measure_combined',
                                split_val = 'incidence',
                             location_col = 'location_id'):
    inc_all, prev_all = split_on_cat(in_df,split_col,split_val)
    inc_dense_only = new_emr.all_fourplus_locs(inc_all,location_col)
    # Combine data dense incidence with all prevalence
    combined_df = pd.concat([prev_all,inc_dense_only],ignore_index = True, sort = True)
    return combined_df
    

# Subsets incidence to data-dense only
# Generates a new SE based on some specified columns
def combined_transform(in_df, num_mean, num_se, denom_mean, denom_se,
                       split_col = "measure_combined",
                       split_val = "incidence",
                    location_col = "location_id"):
    combined_df = subset_incidence_only(in_df)
    combined_df['mean_divided'] = combined_df.mean_allforms / combined_df.mean_latent
    combined_df['se_divided'] = combined_df.apply(lambda x:divided_se(num_mean = x[num_mean],
                                                                        num_se = x[num_se],
                                                                    denom_mean = x[denom_mean],
                                                                      denom_se = x[denom_se]),
                                                                          axis = 1)
    return combined_df


# Mean across columns    
def df_mean(df, mean_col_name, input_cols):
    df[mean_col_name] = df[input_cols].mean(axis = 1)
    return df


# Round a column to the base digits
def df_round(df, round_col_name, input_col, base = 1):

    def rounder(x, base = base):
        return int(base * round(float(x) / base))

    df[round_col_name] = df[input_col].apply(rounder, args = (base,))
    return df
   

# Given sex-specific model results, creates a both-sex dataframe with the same
#  information
def both_sex_model_results(by_sex_df):
    # Get the set of locations, age groups, and years in the input df
    in_age_groups = by_sex_df.age_group_id.unique().tolist()
    in_years = by_sex_df.year_id.unique().tolist()
    in_locs = by_sex_df.location_id.unique().tolist()
    # Get populations for both male and female for each of these age groups,
    #  years, and locations
    pops = get_population(age_group_id = in_age_groups,
                           location_id = in_locs,
                               year_id = in_years,
                                sex_id = [1,2],
                           decomp_step = 'step2')
    # Drop non-useful columns from the populations df
    # pops = pops.drop(labels=['process_version_map_id'],axis=1)
    # Merge the input dataframe on location, age-group, year, sex
    merge_on_cols = ['location_id','year_id','age_group_id','sex_id']
    merged = by_sex_df.merge(pops, on = merge_on_cols, how = 'inner')
    # Now, multiply the 'mean', 'lower', and 'upper' columns by the population
    #   to get a total COUNT for each location-year-age-sex group
    merged['mean_count'] = merged['mean'] * merged['population']
    merged['lower_count'] = merged['lower'] * merged['population']
    merged['upper_count'] = merged['upper'] * merged['population']
    # Drop the old rate columns and the sex column
    merged = merged.drop(labels = ['mean','lower','upper','sex_id'], axis = 1)
    # Group by all columns except mean_count, lower_count, upper_count,
    #   and population
    dont_group_by = ['mean_count','lower_count','upper_count','population']
    group_by_these = [i for i in merged.columns.tolist() if i not in dont_group_by]
    # Get the sum of the counts AND the total population
    # This is the both-sex df
    both_sex = merged.groupby(group_by_these).sum().reset_index()
    both_sex['sex_id'] = 3
    # Calculate the rate as count/population
    both_sex['mean'] = both_sex['mean_count'] / both_sex['population']
    both_sex['lower'] = both_sex['lower_count'] / both_sex['population']
    both_sex['upper'] = both_sex['upper_count'] / both_sex['population']
    # Drop the count columns
    both_sex = both_sex.drop(labels = ['mean_count','lower_count',
                                     'upper_count','population'],
                               axis = 1)
    # Return the both sexes dataframe
    return both_sex


def both_sex_model_results_figures(by_sex_df):
    # Get the set of locations, age groups, and years in the input df
    in_age_groups = by_sex_df.age_group_id.unique().tolist()
    in_years = by_sex_df.year_id.unique().tolist()
    in_locs = by_sex_df.location_id.unique().tolist()
    # Get populations for both male and female for each of these age groups,
    #  years, and locations
    pops = get_population(age_group_id = in_age_groups,
                           location_id = in_locs,
                               year_id = in_years,
                                sex_id = [1,2],
                           decomp_step = 'step2')
    # Merge the input dataframe on location, age-group, year, sex
    merge_on_cols = ['location_id','year_id','age_group_id','sex_id']
    merged = by_sex_df.merge(pops, on = merge_on_cols, how = 'inner')
    # Now, multiply the 'mean', 'lower', and 'upper' columns by the population
    #   to get a total COUNT for each location-year-age-sex group
    merged['mean_count'] = merged['mean'] * merged['population']
    merged['lower_count'] = merged['lower'] * merged['population']
    merged['upper_count'] = merged['upper'] * merged['population']
    # Drop the old rate columns and the sex column
    merged = merged.drop(labels = ['mean','lower','upper','sex_id'], axis = 1)
 
    # Group by all columns except mean_count, lower_count, upper_count,
    #   and population
    dont_group_by = ['mean_count','lower_count','upper_count','population']
    group_by_these = [i for i in merged.columns.tolist() if i not in dont_group_by]
    # Get the sum of the counts AND the total population
    # This is the both-sex df
    both_sex = merged.groupby(group_by_these).sum().reset_index()
    both_sex['sex_id'] = 3
    return both_sex


# This function returns results for sex ids 1, 2, and 3, (male, female, AND both)
# Inputs: model_version_id and measure_id (for get_model_results)
def enhanced_get_model_results(model_version_id, measure_id = 5):
    get_age_groups = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                      22,30,31,32,33,235]
    # Get the sex specific results
    results_by_sex = get_model_results('epi',
                                       model_version_id = model_version_id,
                                             measure_id = measure_id,
                                           age_group_id = get_age_groups,
                                        location_set_id = 35,
                                            decomp_step = 'step2')
    # Get the both-sex results
    both_sex_results = both_sex_model_results(results_by_sex)
    results_combined = pd.concat([results_by_sex,both_sex_results], sort = True)
    return results_combined    


def LTBI_data_adjust(df):
    # Prepare sheet for upload
    df = df.drop(['standard_error',
                  'mean_allforms','lower_allforms','upper_allforms',
                  'mean_latent','lower_latent','upper_latent','measure_latent',
                  'se_allforms','se_latent',
                  'measure_id','model_version_id'], axis = 1)
    df['lower'] = ''
    df['upper'] = ''
    df['crosswalk_parent_seq'] = df['seq']
    df['seq'] = np.nan
    df = df.rename(columns = {'mean_combined':'mean',
                              'measure_allforms':'measure',
                              'se_combined':'standard_error'})
    print('(1)mean max, (2)mean min, (3)std.error max, (4)std.error min.')
    print('If any are less than 0 or greater than 1, you have a problem')
    print(df['mean'].max(),
          df['mean'].min(),
          df['standard_error'].max(),
          df['standard_error'].min())
    return(df)
 

# Joins an in_df to results from a Epi Model Version
def join_to_latent(in_df, model_version_id):
    latent_df = enhanced_get_model_results(model_version_id = model_version_id,
                                                 measure_id = 5)
    if 'year_id' not in in_df:
        in_df = df_mean(in_df, mean_col_name = 'years_average',
                                  input_cols = ['year_start','year_end'])
        in_df = df_round(in_df,round_col_name= 'year_id',
                                    input_col= 'years_average',
                                         base= 5)
    # Join latent to in_df
    join_on = ['year_id','location_id','sex_id','age_group_id']
    join_df = in_df.merge(latent_df,on = join_on, 
                                   how = 'inner',
                              suffixes = ('_allforms','_latent'))

    # A function that either uses existing standard error, or else creates the new 
    #  standard error using known upper and lower values
    def get_new_se(in_se, in_lower, in_upper, in_mean, in_sample_size=None):
        if np.isnan(in_se):
            if np.isnan(in_lower):
                return np.sqrt(in_mean * (1-in_mean) / in_sample_size)
            else:
                return (in_upper - in_lower) / (2*1.96)
        else:
            return in_se

    join_df['se_allforms'] = join_df.apply(lambda x: get_new_se(x['standard_error'],
                                                                x['lower_allforms'],
                                                                x['upper_allforms'],
                                                                x['mean_allforms'],
                                                                x['sample_size']),
                                                                axis = 1)

    join_df['se_latent'] = np.nan
    join_df['se_latent'] = join_df.apply(lambda x: get_new_se(x['se_latent'],
                                                              x['lower_latent'],
                                                              x['upper_latent'],
                                                              x['mean_latent']),
                                                              axis = 1)

    join_df['mean_combined'] = join_df['mean_allforms'] / join_df['mean_latent']
    join_df['se_combined'] = join_df.apply(lambda x: divided_se(num_mean = x['mean_allforms'],
                                                                  num_se = x['se_allforms'],
                                                              denom_mean = x['mean_latent'],
                                                                denom_se = x['se_latent']),
                                                                    axis = 1)
    join_clean_df = LTBI_data_adjust(join_df)
    return join_clean_df


    #filepath shoulde be LTBI dismod outputs adjusted into new custom age groups.
    #Under-5 age-group data
def join_to_latent_custom_ages(in_df,
     filepath = "FILEPATH/ltbi_custom_ages_493040.csv"):
    latent_df = pd.read_csv(filepath)
    if 'year_id' not in in_df:
        in_df = df_mean(in_df, mean_col_name = 'years_average',
                                  input_cols = ['year_start','year_end'])
        in_df = df_round(in_df,round_col_name= 'year_id',
                                    input_col= 'years_average',
                                         base= 5)
    
    # Join latent to in_df
    join_on = ['year_id','location_id','sex_id','age_group_id']
    join_df = in_df.merge(latent_df,on = join_on,
                                   how = 'inner',
                              suffixes = ('_allforms','_latent'))

    # A function that either uses existing standard error, or else creates the new 
    #  standard error using known upper and lower values
    #  If no upper and lower exist, calculate standard error from mean and sample size
    def get_new_se(in_se, in_lower, in_upper, in_mean, in_sample_size=None):
        if np.isnan(in_se):
            if np.isnan(in_lower):
                return np.sqrt(in_mean * (1 - in_mean) / in_sample_size)
            else:
                return (in_upper - in_lower) / (2 * 1.96)
        else:
            return in_se

    join_df['se_allforms'] = join_df.apply(lambda x: get_new_se(x['standard_error'],
                                                                x['lower_allforms'],
                                                                x['upper_allforms'],
                                                                x['mean_allforms'],
                                                                x['sample_size']),
                                                                axis = 1)

    join_df['se_latent'] = join_df.apply(lambda x: get_new_se(x['se_latent'],
                                                              x['lower_latent'],
                                                              x['upper_latent'],
                                                              x['mean_latent']),
                                                                axis = 1)

    join_df['mean_combined'] = join_df['mean_allforms'] / join_df['mean_latent']
    join_df['se_combined'] = join_df.apply(lambda x: divided_se(num_mean = x['mean_allforms'],
                                                                  num_se = x['se_allforms'],
                                                              denom_mean = x['mean_latent'],
                                                                denom_se = x['se_latent']),
                                                                    axis = 1)
    join_clean_df = LTBI_data_adjust(join_df)
    return join_clean_df

####################################
####################################
# Find the allforms standard error from uploaded and latent data
def backcalculate_allforms_se(num_new, num_mean, denom_mean, denom_se):
    new_se = 0
    if num_mean != 0:
        # Formula for calculating the combined SE
        new_se = np.sqrt((num_new**2) * 
                         (denom_mean**2) - 
                         ((num_mean**2) * (denom_se**2) / (denom_mean**2)))
    return new_se

#####################################

def multiply_by_latent(in_df, model_version_id):
    latent_df = enhanced_get_model_results(model_version_id = model_version_id,
                                                 measure_id = 5)
    if 'year_id' not in in_df:
        in_df = df_mean(in_df, mean_col_name = 'years_average',
                                  input_cols = ['year_start','year_end'])
        in_df = df_round(in_df,round_col_name= 'year_id',
                                    input_col= 'years_average',
                                         base= 5)
    
    # Join latent to in_df
    join_on = ['year_id','location_id','sex_id','age_group_id']
    join_df = in_df.merge(latent_df,on = join_on,
                                   how = 'inner',
                              suffixes = ('_uploaded','_latent'))

    # A function that either uses existing standard error, or else creates the new 
    #  standard error using known upper and lower values
    def get_new_se(in_se, in_lower, in_upper):
        if np.isnan(in_se):
            return (in_upper - in_lower) / (2*1.96)
        else:
            return in_se

    join_df['se_latent'] = (join_df['upper_latent'] - join_df['lower_latent']) / (2*1.96)
    join_df['mean_allforms'] = join_df['mean_uploaded'] * join_df['mean_latent']
    join_df['se_allforms'] = join_df.apply(lambda x: backcalculate_allforms_se(num_new = x['standard_error'],
                                                      num_mean = x['mean_allforms'],
                                                    denom_mean = x['mean_latent'],
                                                      denom_se = x['se_latent']),
                                                          axis = 1)
    return join_df

    # Find the new standard error generated by dividing two distributions



