"""
map the ICG level data to bundle

aggregate to 5 year bands

"""
import sys
import time
import itertools
import multiprocessing
import functools
import warnings
import os
import pandas as pd
import numpy as np
from db_queries import get_population, get_covariate_estimates, get_cause_metadata
from db_tools.ezfuncs import query

from clinical_info.Mapping import clinical_mapping
from clinical_info.Functions import hosp_prep, data_structure_utils as dsu

def expand_bundles(df, drop_null_bundles=True):
    """
    This Function maps baby sequelae to Bundles.
    When our data is at the baby seq level there are no duplicates, one icd code
    goes to one baby sequela. At the bundle level there are multiple bundles
    going to the same ICD code so we need to duplicate rows to process every
    bundle. This function maps bundle ID onto baby sequela and duplicates rows.
    Then it can drop rows without bundle ID.  It does this by default.

    Parameters:
        df: Pandas DataFrame
            Must have nonfatal_cause_name column
        drop_null_bundles: Boolean
            If true, will drop the rows of data with baby sequelae that do not
            map to a bundle.
    """

    assert "bundle_id" not in df.columns, "bundle_id has already been attached"
    assert "icg_id" in df.columns,\
        "'icg_id' must be a column"

    df = clinical_mapping.map_to_gbd_cause(df,
                                           input_type='icg',
                                           output_type='bundle',
                                           write_unmapped=False,
                                           truncate_cause_codes=False,
                                           extract_pri_dx=False,
                                           prod=True)

    if drop_null_bundles:
        # drop rows without bundle id
        df = df[df.bundle_id.notnull()]
    assert df.shape[0] > 0, "All the data has been lost in this function"
    return df

def drop_cols(df, propogate_2017=True):
    """
    icg id and name are dropped because the data has been mapped to bundle
    sample size and cases are dropped because they will be re-made when
    aggregating to 5 years. They also contain Null values which lose rows
    in a groupby
    """
    if propogate_2017:
        master_drop = ['icg_id', 'icg_name', 'sample_size', 'cases']
    else:
        master_drop = ['icg_id', 'icg_name', 'sample_size', 'cases', 'lower', 'upper']

    to_drop = [c for c in master_drop if c in df.columns]
    df.drop(to_drop, axis=1, inplace=True)
    assert df.shape[0] > 0, "All the data has been lost in this function"
    return df

def agg_to_bundle(df, propogate_2017=True):

    df = expand_bundles(df)

    df = drop_cols(df, propogate_2017=propogate_2017)

    if propogate_2017:
        sum_cols = ['mean', 'upper', 'lower']
    else:
        sum_cols = ['mean']

    grouper = df.columns.drop(sum_cols).tolist()

    sum_dict = dict(list(zip(sum_cols, ['sum'] * len(sum_cols))))

    df = df.groupby(grouper).agg(sum_dict).reset_index()
    assert df.shape[0] > 0, "All the data has been lost in this function"
    return df

def make_square(df, try_to_break=False):
    """
    Function that inserts zeros for demographic and etiolgy groups that were not
    present in the source.  A zero is inserted for every age/sex/etiology
    combination for each location and for only the years available for that
    location. If a location has data from 2004-2005 we will create explicit
    zeroes for those years but not for any others.

    Age and sex restrictions are applied after the data is made square.

    Parameters:
        df: Pandas DataFrame
            Must be aggregated and collapsed to the bundle level.
    """
    start_square = time.time()
    assert "bundle_id" in df.columns, "'bundle_id' must exist"
    assert "icg_id" not in df.columns, ("Data must not be at the",
                                        " icg_id level")

    # create a series of sorted, non-zero mean values to make sure
    # the func doesn't alter anything
    check_mean = df.loc[df['mean'] > 0, 'mean'].sort_values().\
        reset_index(drop=True)
    if try_to_break:
        # confirmed this breaks the assert below.. which means the assert is 
        # working as we'd expect. A better assert fail output would be good
        check_mean.iloc[3] = .030214

    # square the dataset
    # do we need to make it all square? or just 1 col? Try to do it a faster way below
    df = make_zeros(df, etiology='bundle_id', cols_to_square=['mean'], n_pools=15)

    # assert the sorted means are identical
    assert (check_mean == df.loc[df['mean'] > 0, 'mean'].sort_values().\
        reset_index(drop=True)).all()
    print("Data squared in {} min".format((time.time() - start_square)/60))

    # re-apply age sex restrictions after making data square
    df = clinical_mapping.apply_restrictions(df, age_set='age_group_id', cause_type='bundle')

    # assign upper and lower to 0 where mean is 0
    df.loc[df['mean'] == 0, ['upper', 'lower']] = 0

    print("Data made square and age/sex restrictions applied in {} min".\
        format((time.time() - start_square)/60))

    return df

def merge_population(df, gbd_round_id, decomp_step):
    """
    Function that attaches population info to the DataFrame.  Checks that there
    are no nulls in the population columns.  This has to be ran after the data
    has been made square!

    Parameters:
        df: Pandas DataFrame
    """

    assert (df['mean'] == 0).any(), """There are no rows with zeros, implying
        that the data has not been made square.  This function should be ran
        after the data is square"""

    # create age/year/location lists to use for pulling population
    age_list = list(df.age_group_id.unique())
    loc_list = list(df.location_id.unique())
    year_list = list(df.year_id.unique())

    # pull population
    pop = get_population(age_group_id=age_list, location_id=loc_list,
                         sex_id=[1, 2], year_id=year_list,
                         gbd_round_id=gbd_round_id, decomp_step=decomp_step)

    demography = ['location_id', 'year_id', 'sex_id',
                  'age_group_id']

    pre_shape = df.shape[0]  # store for before comparison
    # then merge population onto the hospital data

    df = df.merge(pop, how='left', on=demography)  # attach pop info to hosp
    assert pre_shape == df.shape[0], "number of rows don't match after merge"

    # assert that there are no nulls in population column:
    hosp_prep.report_if_merge_fail(df, check_col="population",
                                   id_cols=demography, store=True,
                                   filename="population_merge_failure")

    return df

def rate_to_count(df):
    """
    Converts information in ratespace to countspace by mutliplying rates by
    population.

    Parameters:
        df: Pandas DataFrame
            Must already have population information.
    """

    assert "population" in df.columns, "'population' has not been attached yet."
    assert {'mean', 'upper', 'lower'}.issubset(df.columns)

    # we want to take every single rate column into count space
    cols_to_count = df.filter(regex="^mean|^upper|^lower").columns

    # do the actual conversion
    for col in cols_to_count:
        df[col + "_count"] = df[col] * df['population']


    # Drop the rate columns
    df.drop(cols_to_count, axis=1, inplace=True)

    return df

def expandgrid(*itrs):
    # create a template df with every possible combination of
    #  age/sex/year/location to merge results onto
    # define a function to expand a template with the cartesian product
    product = list(itertools.product(*itrs))
    return({'Var{}'.format(i+1):[x[i] for x in product] for i in range(len(itrs))})

def pooled_square_builder(loc, est_types, ages, sexes, loc_bundle, loc_years, eti_est_df, etiology):

    cause_type = loc_bundle.loc[loc_bundle.location_id == loc, etiology].unique().tolist()

    dat = pd.DataFrame(expandgrid(ages, sexes, [loc],
                                  loc_years[loc_years.location_id == loc].
                                  year_id.unique(),
                                  est_types,
                                  cause_type))
    dat.columns = ['age_group_id', 'sex_id', 'location_id', 'year_id', 'estimate_id', etiology]
    dat = dat.merge(eti_est_df, how='left', on=[etiology, 'estimate_id'])
    dat = dat[dat['keep'] == 1]
    dat.drop('keep', axis=1, inplace=True)

    return dat

def make_zeros(df, cols_to_square, etiology='bundle_id', n_pools=15):
    """
    takes a dataframe and returns the square of every
    age/sex/cause type(bundle or baby seq) which
    exists in the given dataframe but only the years
    available for each location id
    """
    if type(cols_to_square) == str:
        cols_to_square = [cols_to_square]

    # badsqr_df = []
    ages = df.age_group_id.unique()
    sexes = df.sex_id.unique()
    est_types = df['estimate_id'].unique().tolist()

    # get all the unique bundles by source and location_id
    src_loc = df[['source', 'location_id']].drop_duplicates()
    src_bundle = df[['source', etiology]].drop_duplicates()
    loc_bundle = src_bundle.merge(src_loc, how='outer', on='source')
    loc_years = df[['year_id', 'location_id']].drop_duplicates()

    eti_est_df = df[[etiology, 'estimate_id']].drop_duplicates()
    eti_est_df['keep'] = 1

    partial_builder = functools.partial(pooled_square_builder, ages=ages,
                        sexes=sexes, loc_years=loc_years, loc_bundle=loc_bundle,
                        eti_est_df=eti_est_df, etiology=etiology, est_types=est_types)
    locs = df.location_id.unique()

    print("booting up multiproc squaring...")
    warnings.warn("This script uses multiprocessing with {} pools. If ran on fair cluster then it needs to be ran with at least {} threads.".format(n_pools))
    p = multiprocessing.Pool(n_pools)
    sqr_df = p.map(partial_builder, locs)
    sqr_df = pd.concat(sqr_df, ignore_index=True, sort=False)
    print("done with that")

    # get a key to fill in missing values for newly created rows
    missing_col_key = df[['location_id', 'year_id', 'source_type_id', 'nid', 'representative_id', 'source']].copy()
    missing_col_key.drop_duplicates(inplace=True)

    # inner merge sqr and missing col key to get all the column info we lost
    sqr_df = sqr_df.merge(missing_col_key, how='inner', on=['location_id', 'year_id'])

    # left merge on our built out template df to create zero rows
    sqr_df = sqr_df.merge(df, how='left', on=['age_group_id', 'sex_id', 'location_id', 'year_id',
                                                  'estimate_id', etiology, 'source_type_id', 'nid',
                                                  'representative_id', 'source'])

    # fill missing values of the col of interest with zero
    for col in cols_to_square:
        sqr_df[col] = sqr_df[col].fillna(0)

    return sqr_df

def split_covered_sources(df, full_coverage_sources):
    """
    Function that splits the data into two pieces: one for sources that
    fully cover their population, and another for sources that do not.
    The data is split based on a list of sources that is hard coded into
    this function.

    RETURNS TWO DATAFRAMES, the SECOND will be the dataframe with fully
    covered populations

    Example call:
        df, covered_df = split_covered_sources(df)

    Parameters:
        df: Pandas DataFrame
    """
    # make condition mask that indicates rows that have full coverage
    has_full_coverage = df.source.isin(full_coverage_sources)

    # make dataframe that only contains fully covered sources
    covered_df = df[has_full_coverage].copy()

    # drop this data from the main dataframe
    df = df[~has_full_coverage]

    return df, covered_df

def merge_denominator(df, run_id):
    """
    Merges denominator stored from when cause fractions were run.  Returns df
    with a column named "denominator".  Meant to be run after data has been
    made square.  There WILL be null values of denominator after the merge.
    This would happen whenever there was a demographic that literally had no
    admissions whatsoever in a data source.  At the moment it is an open
    question whether or not we should drop these null rows, or give them a
    'denominator' value of zero.

    This is a function susceptible to problems. If create_cause_fractions was
    ran a while ago, the information may not match or be up to date anymore.
    Making the denominator should happen right before this function is used.

    The purpose of this is so that the rows that are inserted as zero to make
    the data square can have information for uncertainty for the inserted rows.
    Later, 'denominator' will be renamed to sample_size.

    Parameters:
        df: Pandas DataFrame
    """
    denom_path = FILEPATH.format(run_id)
    warnings.warn("""

                  ENSURE THAT THE DENOMINATORS FILE IS UP TO DATE.
                  the file was last edited at {}

                  """.format(time.strftime('%Y-%m-%d %H:%M:%S',
                  time.localtime(os.path.getmtime(denom_path
                  )))))
    denoms = pd.read_csv(denom_path)
    denoms = dsu.year_id_switcher(denoms)

    pre = df.shape[0]
    df = df.merge(denoms,
                  how='left', on=["age_group_id", "sex_id",
                                  "year_id", "location_id"])

    assert pre == df.shape[0], "number of rows changed after merge"

    # check that data is zero when denominator is null
    assert (df.loc[df.denominator.isnull(), 'mean_count'] == 0).all(), ("mean_count"
        " should be 0")
    assert (df.loc[df.denominator.isnull(), 'lower_count'] == 0).all(), ("lower_count"
        " should be 0")
    assert (df.loc[df.denominator.isnull(), 'upper_count'] == 0).all(), ("upper_count"
        " should be 0")

    # This removes rows that were inserted for demographic groups that never
    # existed in hospital data sources.
    print("pre null denom drop shape", df.shape)
    df = df[df.denominator.notnull()]
    print("post null denom", df.shape)

    return df

def aggregate_to_dismod_years(df):
    """
    Function to collapse data into 5-year bands.  This is done in order to
    reduce the strain on on dismod by reducing the number of rows of data. Data
    must already be square! Data must be in count space!

    Before we have years like 1990, 1991, ..., 2014, 2015.  After we will have:
        1988-1992
        1993-1997
        1998-2002
        2003-2007
        2008-2012
        2013-2017

    Parameters:
        df: Pandas DataFrame
            Contains data at the bundle level that has been made square.
    """

    assert (df.mean_count == 0).any(), """There are no rows with zeros, implying
        that the data has not been made square.  This function should be ran
        after the data is square"""

    assert {"mean_count", "upper_count", "lower_count"}.\
        issubset(df.columns), """The columns mean_count, upper_count,
            lower_count, are not all present, implying that the data has
            not been coverted into count space."""

    if df['upper_count'].max() > df['population'].max():
        warnings.warn("Some CFs are insanely large so we're matching them to 5 times population")
        df.loc[df['upper_count'] > df['population'] * 5, 'upper_count'] =\
            df.loc[df['upper_count'] > df['population'] * 5, 'population'] * 5

    # we need to switch back to year start end
    df.rename(columns={'year_id': 'year_start'}, inplace=True)
    df['year_end'] = df['year_start']

    # first change years to 5-year bands
    df = hosp_prep.year_binner(df)

    # we are mixing data from different NIDs together.  There are special NIDs
    # for aggregated data that the Data Indexers 
    # makes for us.
    warnings.warn("update this to use the new source table in the database")
    df = hosp_prep.apply_merged_nids(df, assert_no_nulls=False, fillna=True)

    if "sample_size" in df.columns:
        # then df already has sample size
        df.drop("sample_size", axis=1, inplace=True)

    # rename "denominator" to "sample_size" in df (not covered_df)
    df.rename(columns={"denominator": "sample_size"}, inplace=True)

    cols_to_sum = df.filter(regex="count$").columns.tolist()
    cols_to_sum = cols_to_sum + ['population', 'sample_size']
    pre_cases = 0
    for col in cols_to_sum:
        pre_cases += df[col].sum()

    groups = ['location_id', 'year_start', 'year_end', 'age_group_id', 'sex_id',
              'nid', 'representative_id', 'bundle_id', 'estimate_id']

    sum_dict = dict(list(zip(cols_to_sum, ["sum"] * len(cols_to_sum))))
    df = df.groupby(groups).agg(sum_dict).reset_index()

    post_cases = 0
    for col in cols_to_sum:
        post_cases += df[col].sum()

    accepted_loss = 300
    lost_cases = abs(pre_cases - post_cases)
    print("The acceptable loss is {} cases. {} cases were lost to floating point error".format(accepted_loss, lost_cases))
    assert lost_cases < accepted_loss,\
        ("some cases were lost. "
         "From {} to {}. A difference of {}".format(pre_cases, post_cases, abs(pre_cases - post_cases)))

    # set sample size to np.nan when mean/upper/lower are greater than 0
    df.loc[(df['mean_count'] > 0) &
           (df['lower_count'] > 0) &
           (df['upper_count'] > 0), 'sample_size'] = np.nan

    return df

def aggregate_covered_sources(covered_df):
    """
    Aggregate sources with full coverage into 5 year bands. Before we have
    years like 1990, 1991, ..., 2014, 2015.  After we will have:
        1988-1992
        1993-1997
        1998-2002
        2003-2007
        2008-2012
        2013-2017

    Parameters:
        covered_df: Pandas DataFrame
            df containing the covered sources.
    """
    assert "mean_count" in covered_df.columns, """The column 'mean_count' is
        not present in the columns of covered_df, implying that the df hasn't
        been transformed into count space."""

    # we need to switch back to year start end
    covered_df.rename(columns={'year_id': 'year_start'}, inplace=True)
    covered_df['year_end'] = covered_df['year_start']

    # first change years to 5-year bands
    covered_df = hosp_prep.year_binner(covered_df)

    # we are mixing data from different NIDs together.  There are special NIDs
    # for aggregated data that the Data Indexers
    # makes for us.
    warnings.warn("update this to use the new source table in the database for UTLAs")
    covered_df = hosp_prep.apply_merged_nids(covered_df)

    # These are the same groups as for the other df.
    groups = ['location_id', 'year_start', 'year_end', 'age_group_id', 'sex_id',
              'nid', 'representative_id', 'bundle_id', 'estimate_id']

    pre_cases = covered_df['mean_count'].sum()

    # sample_size has some null values from being made square, but it was just
    # population, so we're using pop instead. so remember,
    # population == sample_size for covered_df
    covered_df = covered_df.groupby(groups)\
        .agg({'mean_count':'sum',
              'population': 'sum'}).reset_index()

    assert round(pre_cases, 0) == round(covered_df['mean_count'].sum(), 0),\
        ("some cases were lost. "
         "From {} to {}".format(pre_cases, covered_df['mean_count'].sum()))

    return covered_df

def count_to_rate(df, covered_df, review):
    """
    Function that transforms from count space to rate space by dividing by
    dividing by the (aggregated) population.  Returns two dataframes.  This
    performs the transformation for each dataframe separatly.
    Pass in both dataframes!  The order of the returne dataframes is important!
    The main part of hosital data, i.e., sources without full coverage, is
    returned first.  Second item returned is the fully covered sources.

    Example call:
        not_fully_covered, full_covered = count_to_rate(df, covered_df)

    Parameters:
        df: Pandas DataFrame
            df containing most hospital data, the sources that are not fully
            covered.
        covered_df: Pandas DataFrame
            df containing the covered sources.
    """

    # get a list of the count columns
    cnt_cols = df.filter(regex="count$").columns

    for col in cnt_cols:
        new_col_name = col[:-6]
        print(new_col_name, col)
        df[new_col_name] = df[col] / df['population']

    if not review:
        # drop the count columns
        df.drop(cnt_cols, axis=1, inplace=True)

    # for covered_df
    covered_df['mean'] = covered_df['mean_count'] / covered_df['population']
    covered_df.drop('mean_count', axis=1, inplace=True)

    return df, covered_df

def concat_covered_sources(df, covered_df, review):
    """
    Reattaches the data that was split up into fully covered and not fully
    covered sources.  Both data frames should be in rate space!

    Parameters:
        df: Pandas DataFrame
            df containing most hospital data, the sources that are not fully
            covered.
        covered_df: Pandas DataFrame
            df containing the covered sources.
    """
    if not review:
        # this drops the population that was merged on for coverting to counts.
        df.drop('population', axis=1, inplace=True)  # don't need pop anymore

    covered_df.rename(columns={"population": "sample_size"}, inplace=True)

    # bring covered_df and df together.
    # where we have full coverage, lower and upper should be null
    # mean will never be null
    df = pd.concat([df, covered_df], ignore_index=True, sort=False)

    return df


def get_parent_injuries(df):
    """
    roll up the child causes into parents add parent injuries

    Parameters:
        df: Pandas DataFrame
    """
    parent_path = FILEPATH.format(parent_path))
    pc_injuries = pd.read_csv(parent_path)

    # take only parent causes from map
    # pc_injuries = pc_injuries[pc_injuries.parent==1]
    # create BID and ME ID dicts for parent causes
    # TODO do we need meid?
    me_dict = dict(list(zip(pc_injuries.e_code, pc_injuries['level1_meid'])))
    bundle_dict = dict(list(zip(pc_injuries.e_code, pc_injuries['Level1-Bundle ID'])))

    # loop over ecodes pulling every parent name
    df_list = []
    for parent in pc_injuries.loc[pc_injuries['parent']==1, 'e_code']:
        inj_df =\
         pc_injuries[(pc_injuries['baby sequela'].str.contains(parent)) &\
                      (pc_injuries['parent'] != 1)]
        inj_df = inj_df[inj_df['Level1-Bundle ID'].notnull()]
        inj_df['me_id'] = me_dict[parent]
        inj_df['bundle_id'] = bundle_dict[parent]
        # inj_df['nonfatal_cause_name'] = parent
        df_list.append(inj_df)

    parent_df = pd.concat(df_list, sort=False)

    # number of child injuries in csv should match new df
    assert pc_injuries.child.sum() == parent_df.child.sum(),\
        "sum of child causes doesn't match sum of parent causes"

    parent_df.drop(['baby sequela', 'ME name level 1', 'level 1 measure',
                    'e_code', 'parent', 'child'], axis=1,
                    inplace=True)

    # TODO remove MEID parts of this?
    parent_df.rename(columns={'me_id': 'parent_me_id',
                              'level1_meid': 'modelable_entity_id',# child me_id
                              'bundle_id': 'parent_bundle_id',
                              'Level1-Bundle ID': 'bundle_id'},# child bundle_id
                     inplace=True)

    # reorder cuz it's hard to look at
    col_before = parent_df.columns
    parent_df = parent_df[['modelable_entity_id', 'parent_me_id', 'bundle_id',
                           'parent_bundle_id']].copy()
    assert set(col_before) == set(parent_df.columns),\
        'you dropped a column while reordering columns'

    # TODO change to bundle ID
    # duplicate data for parent injuries
    parent_df = parent_df.merge(df, how='left', on='bundle_id')

    # TODO remove meids?
    # drop child ME and bundle IDs
    parent_df.drop(['modelable_entity_id', 'bundle_id'], axis=1, inplace=True)
    # add parent bundle and ME IDs
    parent_df.rename(columns={'parent_me_id': 'modelable_entity_id',
                              'parent_bundle_id': 'bundle_id'}, inplace=True)

    # don't need modelable_entity_id, it was just along for the ride
    parent_df.drop('modelable_entity_id', axis=1, inplace=True)

    assert set(parent_df.columns).symmetric_difference(set(df.columns)) == set()

    groups = ['location_id', 'year_start', 'year_end',
              'age_group_id', 'sex_id', 'nid',
              'representative_id', 'estimate_id',
              'bundle_id'] # re-collapsing on nonfatal cause will double count

    cols_to_sum = df.filter(regex="^mean|^upper|^lower").columns.tolist()
    cols_to_sum = cols_to_sum + ['sample_size']
    # sum parent rows of data together
    parent_df = parent_df.groupby(groups).agg(\
        dict(list(zip(cols_to_sum, ['sum'] * len(cols_to_sum))))).reset_index() # keep sample size

    # set sample size to np.nan when upper exists and isn't zero
    parent_df.loc[(parent_df.upper != 0) &\
                (parent_df.upper.notnull()) &\
                (parent_df.sample_size.notnull()), 'sample_size'] = np.nan
    # rbind duped parent data back onto hospital data
    df = pd.concat([df, parent_df], sort=False, ignore_index=True).reset_index(drop=True)
    return(df)

def merge_measure(df, drop_mult_measures=True):
    cm = clinical_mapping.get_clinical_process_data("cause_code_icg", prod=True)
    cm = cm[['icg_id', 'icg_measure']].drop_duplicates()
    
    bun = clinical_mapping.get_clinical_process_data('icg_bundle', prod=True)
    measures = bun.merge(cm, how='left', on='icg_id')
    measures = measures[['bundle_id', 'icg_measure']].drop_duplicates()
    if drop_mult_measures:
        doubles = measures.bundle_id[measures.bundle_id.duplicated(keep=False)].unique()
        warnings.warn("We will be dropping distinct measures for bundle ids {}".format(doubles))
        measures.drop_duplicates(subset=['bundle_id'], inplace=True)
    measures.rename(columns={'icg_measure': 'measure'}, inplace=True)
    
    pre = df.shape[0]
    df = df.merge(measures, how='left', on='bundle_id')
    assert pre == df.shape[0], "Row counts have changed, not good"
    
    parent_inj = [264., 269., 272., 276, 270, 275]
    df.loc[df.bundle_id.isin(parent_inj), 'measure'] = 'inc'
    # assert df.measure.isnull().sum() == 0, "There shouldn't be any null measures"
    if df.measure.isnull().sum() != 0:
        no_measure_bundles = df.loc[df.measure.isnull(), 'bundle_id'].unique()
        warnings.warn("There are null measures for bundles {}".format(no_measure_bundles))
    return df

def align_uncertainty(df):
    """
    We're using different forms of uncertainty, depending on estimate type and
    data source so set the values we'd expect and test them
    """
    # A bunch of uncertainty checks
    # When mean_raw is zero, then sample_size should not be null.
    assert df.loc[df['mean'] == 0, "sample_size"].notnull().all(),\
        "Imputed zeros are missing sample_size in some rows."

    # Make upper and lower Null when sample size is not null
    est_types = ['inp-primary-cf1_modeled', 'inp-primary-cf2_modeled',
               'inp-primary-cf3_modeled']
    df.loc[(df['sample_size'].notnull()) &\
           (df['estimate_id'].isin(est_types)),
           ['upper', 'lower']] = np.nan

    return df

def get_5_year_haqi_cf(gbd_round_id, decomp_step, min_treat=0.1, max_treat=0.75):
    """
    A function to get the health access quality covariates data which we'll
    use to divide our mean_raw values by to adjust our estimates up

    Parameters:
        min_treat: float
            minimum access. Sets a floor for the CF. If 0.1 then the lowest possible CF will be 0.1,
            in practice this is a 10x increase in the estimate
        max_treat: float or int
            maximum acess. Sets a cap for the CF. If 75 then any loc/year with a covariate above 75
            will have a CF of 1 and the data will be unchanged
    """
    # get a dataframe of haqi covariate estimates
    df = get_covariate_estimates(covariate_id=1099, gbd_round_id=gbd_round_id, decomp_step=decomp_step)

    df.rename(columns={'year_id': 'year_start'}, inplace=True)

    if df.mean_value.mean() > 1 and max_treat < 1:
        warn_msg = """Increasing max_treat variable 100X. Mean of the HAQi column is larger 
        than 1. We assume this means the range is from 0 to 100. Summary stats for the 
        mean_value column in the haqi covar are \n {}""".format(df.mean_value.describe())
        warnings.warn(warn_msg)
        max_treat = max_treat * 100

    # set the max value
    df.loc[df.mean_value > max_treat, 'mean_value'] = max_treat

    # get min df present in the data
    # Note, should this just be 0.1?
    min_df = df.mean_value.min()

    # make the correction
    df['haqi_cf'] = \
        min_treat + (1 - min_treat) * ((df['mean_value'] - min_df) / (max_treat - min_df))

    # drop the years outside of hosp_data so year binner doesn't break
    df['year_end'] = df['year_start']
    warnings.warn("Currently dropping HAQi values before 1988 and after 2017")
    df = df[(df.year_start > 1987) & (df.year_start < 2018)].copy()
    df = hosp_prep.year_binner(df)

    # Take the average of each 5 year band
    df = df.groupby(['location_id', 'year_start', 'year_end']).agg({'haqi_cf': 'mean'}).reset_index()

    assert df.haqi_cf.max() <= 1, "The largest haqi CF is too big"
    assert df.haqi_cf.min() >= min_treat, "The smallest haqi CF is too small"

    return df

def apply_haqi_corrections(df, gbd_round_id, decomp_step):
    """
    merge the haqi correction (averaged over 5 years) onto the hospital data
    """
    haqi = get_5_year_haqi_cf(gbd_round_id, decomp_step)

    pre = df.shape
    df = df.merge(haqi, how='left', on=['location_id', 'year_start', 'year_end'])
    assert pre[0] == df.shape[0],\
        "DF row data is different. That's not acceptable. Pre shape {}. Post shape {}".\
        format(pre, df.shape)
    assert df.haqi_cf.isnull().sum() == 0,\
        "There are rows with a null haqi value. \n {}".format(\
            df[df.haqi_cf.isnull()])

    return df


def create_source_map(df):
    # create a map from location and year to source
    source_map = df[['location_id', 'year_id', 'source']].drop_duplicates()
    source_map.rename(columns={'year_id': 'year_start'}, inplace=True)
    source_map['year_end'] = source_map['year_start']
    # we have two sources for the same 5 year bin, this creates some problems of duping data
    # for now just overwrite shanghai with NHSIRS, wont' affect underlying data
    source_map.loc[source_map['source'] == "CHN_SHD", 'source'] = "CHN_NHSIRS"
    # bin to 5 year groups so it matches with the data after agg to dismod
    source_map = hosp_prep.year_binner(source_map).drop_duplicates()
    return source_map

def write_estimates(df, run_id):
    # convert mixed type cols to numeric
    df = hosp_prep.fix_col_dtypes(df, errors='raise')
    df['measure'] = df['measure'].astype(str)
    print("Writing the df file...")
    file_path = FILEPATH.format(run_id)
    hosp_prep.write_hosp_file(df, file_path, backup=False)
    return df

def create_bundle_estimates(df, gbd_round_id, decomp_step, run_id):
    source_map = create_source_map(df)

    print("begin create bundle estimates method agg to bundle")
    df = agg_to_bundle(df, propogate_2017=True)
    print("begin cbe.make square")
    df = make_square(df)
    print("begin cbe.merge population")
    df = merge_population(df, gbd_round_id=gbd_round_id, decomp_step=decomp_step)
    print("begin cbe.rate to count")
    df = rate_to_count(df)

    print("splitting covered and uncovered sources")
    df, covered_df = split_covered_sources(df, full_coverage_sources=["UK_HOSPITAL_STATISTICS"])
    print("merge on the denominator for zero obs")
    df = merge_denominator(df, run_id)

    # code is breaking b/c cases are being lost determine if it's floating point error (it was)
    df.loc[df['upper_count'] > df['population'], 'upper_count'] = df.loc[df['upper_count'] > df['population'], 'population']
    print("aggregating to 5 year groups")
    df = aggregate_to_dismod_years(df)
    covered_df = aggregate_covered_sources(covered_df)
    print("going back to rate space")
    df, covered_df = count_to_rate(df, covered_df, review=False)

    # bring the two streams back into one.
    df = concat_covered_sources(df, covered_df, review=False)

    print("adding parent injury bundles")
    df = get_parent_injuries(df)
    for col in ['mean', 'lower', 'upper']:
        hosp_prep.check_parent_injuries(df, col_to_sum=col, verbose=True)

    # merge measures on using icg measures
    df = merge_measure(df)

    # set nulls
    df = align_uncertainty(df)

    # apply the haqi correction
    df = apply_haqi_corrections(df, gbd_round_id=gbd_round_id, decomp_step=decomp_step)

    # bring data sources back
    pre = df.shape[0]
    df = df.merge(source_map, how='left', on=['location_id', 'year_start', 'year_end'])
    assert pre == df.shape[0], "rows changed unexpectedly"
    assert df.source.isnull().sum() == 0, "We are missing some sources"

    # apply the 5 year inj corrections
    df = hosp_prep.apply_inj_corrections(df, run_id)

    # apply e code proportion cutoff, removing rows under cutoff
    df = hosp_prep.remove_injuries_under_cutoff(df, run_id)

    df = write_estimates(df, run_id)

    return df
