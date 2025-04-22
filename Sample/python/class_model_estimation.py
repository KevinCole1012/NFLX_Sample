import pandas as pd
import numpy as np
from datetime import datetime
import pickle
import model.model_common as ecm
import weight.weighted as wgt
from warehouse.utilities.connection import get_connection_list_sa, get_engine_list_sa
from warehouse.utilities.update import upsert_table_sa
from dateutil.relativedelta import relativedelta
from sqlalchemy import text
from tbox.utilities.workflow import workflow_set_state
from tbox.statistics.data import clean_missing, fill_if_all_nan
from tbox.statistics.standardize import standardize
import logging
from tbox.logging import setup_logging
from argparse import ArgumentParser
from tbox.statistics.classification import xgboost_classifier_tune, logistic_classifier_tune, lightgbm_classifier_tune
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, log_loss
from tbox.statistics.correlation import correlation
from tbox.statistics.performance import compute_quantile_returns
from tbox.statistics.regression import linear_regression
from tbox.statistics.neutralize import neutralize
from tbox.statistics.winsorize import winsorize
from model.model_common import fetch_and_pivot_class_numerical_features, fetch_and_pivot_model_exposures
from configuration import (construct_factor_configs, construct_alpha_model_configs, get_alpha_model_alpha_factors,
                           get_alpha_model_risk_factors, get_alpha_model_inv_univ_portfolio_ids,
                           get_alpha_model_factor_params)

setup_logging('debug')

py_config_active = 'dev'

sql_connections = get_connection_list_sa('dev')
sql_connections_dev = get_connection_list_sa('dev')  # TODO: Change hyperparam query to prod when ready
sql_connections_prod = get_connection_list_sa('prod')  # TODO: Change hyperparam query to prod when ready
sql_engines = get_engine_list_sa(py_config_active)

# TODO: Need to get process setup for pulling file paths from 'Equity/Common/src/resources/system/config.yml'

if py_config_active == 'prod':
    class_path = ('###PATH###'
                  '\\Global\\Database\\Classification')
if py_config_active == 'dev':
    class_path = ('###PATH###'
                  '\\Database\\Classification')

workflow_set_state(engine=sql_engines['investment'],
                   context='equity.model.classification.python.estimate',
                   state_id='completed',
                   state='0')

# This script runs the classification modeling portion of the classification signal estimation. This has a dependency
# that the R training datasets have already been persisted to model_class_features in the models DB

# Overwrite should be set to False if you are attempting to generate new data and append it to previously existing data

overwrite = True

# Designate whether permutation feature impact & accumulated local effects analysis should be run

run_feature_impact = False

# Define model_id and version that will be used to import/label our output files

# Input list of model_ids to compute classification models
# in general the list model_id_list = ['global','us','emerging'] and routine iterates through models

# parser = ArgumentParser(
#     prog = 'getAutosysArgs',
#     description = 'Retrieves the arguments from autosys job'
# )
#
# parser.add_argument('-models', '--models')
# parser.add_argument('-modelType', '--modelType')
#
# args = parser.parse_args()
#
# string = args.models
# delimiters = [","]
#
# for delimiter in delimiters:
#     string = " ".join(string.split(delimiter))

# model_id_list = string.split()

print(model_id_list)

model_id_list = ['us']
version_suffix = 'ex_base_signal_starmlag'

# Define begin & end date if not just running for the current periods

begin_date = '2024-1-1'
end_date = '2024-12-1'

# Set training & validation lengths in terms of the number of months for creating our period metadata

training_months = 180
min_training_months = 36
halflife = None
identifier_vars = ['period_dt', 'identifier', 'exchange_country']
positive_classes = ['long', 'short']
max_hyperparam_age = 12
model_list = ['ridge', 'lightgbm', 'xgboost']
n_jobs = -1
compute_performance = True
n_trials = 20
lambda_values = np.logspace(-3, 3, 16)
version = f'{training_months}m_{version_suffix}'

# Logical operations to set variables to production values if not defined specifically above

# Set begin & end dates

if begin_date:
    pass
    if end_date:
        pass
    else:
        end_date = str(datetime.today().date().replace(day=1))
else:
    begin_date = str(datetime.today().date().replace(day=1))
    end_date = str(datetime.today().date().replace(day=1))

begin_datetime = datetime.strptime(begin_date, '%Y-%m-%d')
begin_query_date = begin_datetime + relativedelta(months=-training_months-1)
end_query_date = datetime.strptime(end_date, '%Y-%m-%d')

for model_id in model_id_list:

    factor_config = construct_factor_configs()
    model_config = construct_alpha_model_configs(factors_config=factor_config)

    risk_factors = get_alpha_model_risk_factors(model_config[model_id])
    inv_universe = get_alpha_model_inv_univ_portfolio_ids(model_config[model_id])

    # Grab alpha factors metadata and make relevant modifications
    alpha_factors = get_alpha_model_alpha_factors(model_config[model_id])
    alpha_factors_dummy = get_alpha_model_alpha_factors(model_config[model_id], dummy=True)
    alpha_factor_params = get_alpha_model_factor_params(model_config[model_id], factor_ids=alpha_factors)
    alpha_factors_lagged_for_training = alpha_factor_params.loc[
        alpha_factor_params.lag == 1, 'factor_id'
    ].unique().tolist()
    alpha_factors = ['model_factor_value_' + model_id + '_' + x for x in alpha_factors]
    alpha_factors_dummy = ['model_factor_value_' + model_id + '_' + x for x in alpha_factors_dummy]
    alpha_factors_lagged_for_training = [
        'model_factor_value_' + model_id + '_' + x for x in alpha_factors_lagged_for_training]
    alpha_factors_continuous = [x for x in alpha_factors if x not in alpha_factors_dummy]

    class_model_ids = model_config[model_id]['class_model_ids']
    riskweight_var = model_config[model_id]['alpha_risk_weight']

    return_variables = [
        'market_totalreturn_usd', 'market_totalreturn_usd_win', 'market_totalreturn_usd_win_rfres_rfres',
    ]

    if len(class_model_ids) > 1:
        class_where_clause = "class_model_id in " + str(tuple(class_model_ids))
    else:
        class_where_clause = "class_model_id = '" + str(class_model_ids[0]) + "'"

    # Pull feature names for categoricals & numericals to explicitly assign in the modeling process

    class_categorical_sql = """
                                SELECT * FROM model_class_feature_metadata
                                WHERE ({model_where_clause} and type = 'categorical')
                            """.format(model_where_clause=class_where_clause)

    categorical_features_df = pd.read_sql(text(class_categorical_sql), sql_connections['models'])
    categorical_features = categorical_features_df.feature_id.unique().tolist()

    # Numerical
    class_numerical_sql = """
                        SELECT * FROM model_class_feature_metadata
                        WHERE ({model_where_clause}
                        AND type = 'continuous')
                    """.format(model_where_clause=class_where_clause)
    numerical_features = pd.read_sql(text(class_numerical_sql),
                                     sql_connections['models'])
    numerical_features = numerical_features.feature_id.unique().tolist()

    # Grab class target variable from the metadata
    class_target_sql = """
                                    SELECT * FROM model_class_feature_metadata
                                    WHERE ({model_where_clause} and type = 'target' and category = 'return_side')
                                """.format(model_where_clause=class_where_clause)

    class_target_df = pd.read_sql(text(class_target_sql), sql_connections['models'])
    class_target = class_target_df.feature_id.unique().tolist()[0]

    # Pull numerical data
    class_data_num = fetch_and_pivot_class_numerical_features(
        sql_connection=sql_connections['models'],
        class_where_clause=class_where_clause,
        begin_query_date=str(begin_query_date),
        end_query_date=str(end_query_date),
        numerical_features=numerical_features,
        chunk_size=10
    )

    # Pull categorical data
    class_data_sql_cat = """
                            SELECT * FROM model_class_features_character
                            WHERE ({model_where_clause} and period_dt >= '{query_date_start}'
                            and period_dt <= '{query_date_end}')
                         """.format(model_where_clause=class_where_clause,
                                    query_date_start=begin_query_date,
                                    query_date_end=end_query_date)
    class_data_cat = pd.read_sql(text(class_data_sql_cat), sql_connections['models'])

    # Pivot tall SQL training data to wide format
    class_data_cat = class_data_cat.pivot(index=['class_model_id', 'period_dt', 'identifier', 'exchange_country']
                                          , columns='feature_id'
                                          , values='value').reset_index()

    # Merge numerical data on our categorical data
    class_data = pd.merge(class_data_cat, class_data_num
                          , how='left'
                          , on=['class_model_id', 'period_dt', 'identifier', 'exchange_country']).reset_index(drop=True)

    del class_data_num
    del class_data_cat

    # Convert datetype of numerical features to float32 from string object type
    for feature in numerical_features:
        class_data[feature] = class_data[feature].astype(np.float32)

    # Fill missing by period_dt & GICS6 to get better peer grouping for filling missing
    gics6_hist_sql = """
                            SELECT period_dt
                             ,identifier
                             ,exchange_country
                             ,substring(GICS8_hist, 1, 6) as gics6_hist 
                            FROM cmmn_data_gics_hist
                            WHERE (identifier in {identifier_list} and period_dt >= '{query_date_start}'
                            and period_dt <= '{query_date_end}')
                         """.format(identifier_list=tuple(class_data['identifier'].unique().tolist()),
                                    query_date_start=begin_query_date,
                                    query_date_end=end_query_date)
    gics6_hist = pd.read_sql(text(gics6_hist_sql), sql_connections_prod['investment']).drop_duplicates(
        subset=['period_dt', 'identifier', 'exchange_country'])

    class_data = pd.merge(class_data,
                          gics6_hist,
                          how='left',
                          on=['period_dt', 'identifier', 'exchange_country'])
    class_data.fillna({'gics6_hist': 'missing'}, inplace=True)

    class_data = clean_missing(data=class_data,
                               columns=numerical_features,
                               strategy='median',
                               groups=['period_dt', 'class_model_id', 'gics6_hist'])

    standardized_features = [x for x in numerical_features if x not in alpha_factors_dummy]
    class_data = standardize(data=class_data,
                             columns=standardized_features,
                             groups=['period_dt', 'class_model_id'],
                             normal=False,
                             suffix=None,
                             missing_strategy=None)

    class_data = fill_if_all_nan(data=class_data,
                                 groups=['period_dt', 'class_model_id'],
                                 columns=numerical_features,
                                 value=0)

    # Fill any missing values for categorical variables with missing
    class_data[categorical_features] = class_data[categorical_features].fillna('missing')
    class_data.drop_duplicates(subset=['period_dt', 'identifier', 'exchange_country'], inplace=True)

    # Create categorical features that will work in tree based (ordinal) & linear models (dummy)
    cat_features_ordinal = []
    cat_features_dummy = []
    cat_values = {}

    for cat_feature in categorical_features:
        logging.info('Creating ordinal category features for %s.', cat_feature)
        cat_values[cat_feature] = pd.DataFrame()
        cat_values[cat_feature][cat_feature] = class_data[cat_feature].unique()
        cat_values[cat_feature][cat_feature + '_ord'] = np.arange(len(cat_values[cat_feature])) + 1
        cat_features_ordinal.append(cat_feature + '_ord')

        class_data = pd.merge(class_data, cat_values[cat_feature], how='left', on=[cat_feature])

        logging.info('Creating ordinal category features for %s.', cat_feature)
        cat_features_dummy = cat_features_dummy + class_data[cat_feature].unique().tolist()
        cat_features_temp = class_data.loc[:, identifier_vars + [cat_feature]]
        cat_features_temp['dummy'] = 1
        cat_features_temp = cat_features_temp.pivot(index=identifier_vars,
                                                    columns=cat_feature,
                                                    values='dummy').reset_index()

        class_data = pd.merge(class_data, cat_features_temp, how='left', on=identifier_vars)
        class_data.loc[:, cat_features_dummy] = class_data[cat_features_dummy].fillna(0)

    # Ensure that we don't have any duplicate observations

    class_data.drop_duplicates(subset=['period_dt', 'identifier', 'exchange_country'], inplace=True)
    class_data.sort_values(by=['period_dt', 'class_model_id', 'identifier', 'exchange_country'], inplace=True)

    # Compute exponential weights and normalize
    if halflife is not None:
        weights = wgt.exp_weights(length=training_months,
                                  halflife=halflife,
                                  sort_desc=False,
                                  weight_sum=1)

        weights['weight'] = weights.weight / np.mean(weights.weight)

    # Pull risk factors for return neutralization
    risk_factors = pd.read_sql(text(
        f"""
        SELECT distinct factor 
        FROM cmmn_model_risk_factors
        WHERE model_id = '{model_id}'
        """), sql_connections['investment'])['factor'].unique().tolist()

    # Need to remove the residual return risk factors for exposure neutralization if model = 'us'
    if model_id == 'us':
        values_to_remove = ['u_f_nresret_1m', 'u_f_nresret_2m1']
        risk_factors_exp = [item for item in risk_factors if item not in values_to_remove]
    else:
        risk_factors_exp = risk_factors

    # Pull all risk factors and alpha factors
    model_exp_var_list = risk_factors + [riskweight_var]
    model_factor_exposures = fetch_and_pivot_model_exposures(sql_connection=sql_connections_prod['investment'],
                                                             base_model=model_id,
                                                             model_exp_var_list=model_exp_var_list,
                                                             begin_query_date=str(begin_query_date),
                                                             end_query_date=str(end_query_date),
                                                             identifier_vars=identifier_vars,
                                                             chunk_size=10)

    # Make sure that all risk factors are included in exposures and fill missing dummy vars with 0
    risk_factors = [x for x in risk_factors if x in model_factor_exposures.columns.tolist()]
    risk_factors_exp = [x for x in risk_factors_exp if x in model_factor_exposures.columns.tolist()]

    if len(inv_universe) > 1:
        inv_universe_where_statement = f"portfolio_id in {tuple(inv_universe)}"
    else:
        inv_universe_where_statement = f"portfolio_id = '{inv_universe[0]}'"

    # Grab our universe data for relevant investment universe IDs to compute investable signals & diagnostics
    investment_universe_query = text(
        f"""
        SELECT period_dt, 
            portfolio_id, 
            identifier, 
            excntry as exchange_country
        FROM v_cmmn_benchmarks_monthly 
        WHERE {inv_universe_where_statement}
        AND period_dt >= '{begin_query_date}' 
        AND period_dt <= '{end_date}'
        """
    )
    universe_holdings = pd.read_sql(investment_universe_query, sql_connections_prod["investment"])
    universe_holdings['portfolio_id'] = universe_holdings['portfolio_id'].str.rstrip('_p')
    universe_holdings['portfolio_id'] = universe_holdings['portfolio_id'].str.rstrip('_ai')

    model_factor_exposures = pd.merge(model_factor_exposures, universe_holdings, on=identifier_vars
                                      , how='inner')  # This is set to inner join for shrunk method

    identifier_list = model_factor_exposures.identifier.unique().tolist()

    returns_query = text(
        f"""
        SELECT period_dt,
            identifier,
            exchange_country,
            totalreturn_usd   
        FROM cmmn_data_returns_monthly 
        WHERE identifier IN {tuple(identifier_list)} 
        AND period_dt >= '{begin_query_date}' 
        AND period_dt <= '{end_date}'"""
    )
    returns_data = pd.read_sql(returns_query, sql_connections_prod['investment'])

    # Merge with combined data
    merged_model_data = pd.merge(model_factor_exposures,
                                 returns_data,
                                 how='inner',
                                 on=identifier_vars)

    class_data = pd.merge(class_data,
                          merged_model_data,
                          how='left',
                          on=identifier_vars)

    class_data.drop_duplicates(subset=identifier_vars, inplace=True)

    del model_factor_exposures, returns_data, universe_holdings, merged_model_data

    class_data.loc[:, risk_factors] = class_data.loc[:, risk_factors].fillna(0)
    class_data.fillna({'portfolio_id': 'residual'}, inplace=True)

    # Fill missing returns and riskweight variable by gics6 if possible
    class_data = clean_missing(data=class_data,
                               columns=['totalreturn_usd', riskweight_var],
                               groups=['period_dt', 'class_model_id', 'gics6_hist'],
                               strategy='median')

    # Fill missing returns and riskweight variable by class_model_id
    class_data = clean_missing(data=class_data,
                               columns=['totalreturn_usd', riskweight_var],
                               groups=['period_dt', 'class_model_id'],
                               strategy='median')

    # Winsorize returns
    class_data = winsorize(data=class_data,
                           columns='totalreturn_usd',
                           lower_percentile=0.01,
                           upper_percentile=0.99,
                           groups=['period_dt', 'class_model_id'],
                           suffix='_win')

    risk_factors = [x for x in risk_factors if x in class_data.columns]
    risk_factors_exp = [x for x in risk_factors_exp if x in class_data.columns]
    model_exp_var_list = [x for x in model_exp_var_list if x in class_data.columns]

    # Create new neutralized return across the estimation model
    class_data = neutralize(data=class_data,
                            target_vars='totalreturn_usd_win',
                            neutralization_vars=risk_factors,
                            groups=['period_dt', 'class_model_id'],
                            weight=riskweight_var,
                            intercept=True,
                            missing_strategy='median',
                            residual_suffix='_rfres')

    # Create lagged training features
    class_data.sort_values(by=['identifier', 'exchange_country', 'period_dt'], inplace=True)
    lagged_features_map = {}
    lagged_features = []
    for lagged_feature in alpha_factors_lagged_for_training:
        class_data[lagged_feature + '_1'] = class_data.groupby(
            ['identifier','exchange_country'])[lagged_feature].shift(1)
        lagged_features_map[lagged_feature] = lagged_feature + '_1'
        lagged_features.append(lagged_feature + '_1')

    # Ensure we fill any missing and inf values with 0
    class_data[model_exp_var_list + numerical_features + lagged_features] = class_data[
        model_exp_var_list + numerical_features + lagged_features].fillna(0)
    class_data = class_data.replace([np.inf, -np.inf], 0)

    # Establish the first date of the current month to filter target returns below. We don't know returns for the
    # current month.
    return_cutoff = str(datetime.today().date().replace(day=1))

    # Create target return decile var and cv_groups
    class_data.loc[class_data.period_dt < return_cutoff, 'target_10'] = class_data.loc[
        class_data.period_dt < return_cutoff].groupby(
        ['period_dt', 'class_model_id'])['totalreturn_usd_win_rfres'].transform(
        lambda x: pd.qcut(x, q=10, labels=range(1, 11), duplicates='drop'))
    cv_group_vars = ['period_dt', 'portfolio_id', 'target_10']

    # Grab min period_dt from our data for training start and computing first prediction period

    min_period = class_data.period_dt.min()
    prediction_start_period = min_period + relativedelta(months=min_training_months)

    forecast_periods = pd.DataFrame({'period_dt': class_data['period_dt'].unique()})
    forecast_periods['training_start'] = forecast_periods['period_dt'].shift(training_months)
    forecast_periods['training_start'] = forecast_periods['training_start'].fillna(min_period)
    forecast_periods['training_end'] = forecast_periods['period_dt'].shift(1)
    forecast_periods = forecast_periods.loc[forecast_periods.period_dt >= prediction_start_period]

    # Drop periods where we don't have enough history to satisfy our training window
    forecast_periods = forecast_periods.dropna(subset=['training_start'])

    # Get list of period_dt values to loop over before setting period_dt as the index
    prediction_periods = forecast_periods.loc[
        (forecast_periods['period_dt'] >= begin_date) & (
                forecast_periods['period_dt'] <= end_date), 'period_dt'].to_list()

    # Set period_dt as the date metadata index to loop across when setting up our projects
    forecast_periods = forecast_periods.set_index(['period_dt'])

    # Setup various feature sets for training & prediction
    linear_features = numerical_features
    nonlinear_features = numerical_features + cat_features_ordinal

    class_data.sort_values(by=['period_dt', 'class_model_id', 'identifier', 'exchange_country'], inplace=True)

    for period in prediction_periods:

        logging.info('Running classification models for %s for the holdout period %s.', model_id, period)

        seed = np.random.randint(1, 1000)

        # Output datasets for appending
        (combined_hyperparams,
         combined_summary,
         combined_prob,
         combined_coef,
         combined_trials,
         combined_models) = (pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                             {})

        training_start = forecast_periods.loc[period, 'training_start']
        training_end = forecast_periods.loc[period, 'training_end']

        # Training data for current period
        train_data = class_data.loc[
            (class_data['period_dt'] <= training_end) &
            (class_data['period_dt'] >= training_start)].sort_values(by=identifier_vars)

        # Here we have the added step of overwriting the contemp with lagged features for training purposes where
        # appropriate. We need to do this so the dataframes will align with the prediction set.
        for key, value in lagged_features_map.items():
            train_data[key] = train_data[value]

        # Prediction data for current period
        prediction_data = (class_data.loc[
                               class_data.period_dt == period]
                           .sort_values(by=identifier_vars).reset_index(drop=True))

        if halflife is not None:
            training_weights = pd.DataFrame()
            training_weights['period_dt'] = train_data.period_dt.unique()
            training_weights = pd.merge(training_weights.sort_values(by=['period_dt']).reset_index(drop=True),
                                        weights,
                                        left_index=True,
                                        right_index=True)
            weight_column = 'weight'
            train_data = pd.merge(train_data,
                                  training_weights,
                                  how='left',
                                  on='period_dt')

        else:
            weight_column = None

        min_hyperparam_date = period + relativedelta(months=-max_hyperparam_age)

        for target_class in positive_classes:

            # The dates below run off the class_data_periods metadata that was established earlier

            logging.info('Generating %s side class model for period %s.', target_class, period)

            train_data.loc[train_data[class_target] == target_class, target_class] = 1
            train_data.loc[train_data[class_target] != target_class, target_class] = 0

            # Downsampling of majority class which we know to be "none" by grabbing the proper number we'd like to
            # sample and use scaler term for scale_pos_weight parameter
            class_counts = train_data.loc[:, [target_class, 'identifier']].groupby(
                [target_class], as_index=False).count()
            major_class_count = int(
                round(class_counts.loc[class_counts[target_class] == 0, 'identifier'].iloc[0]))
            minor_class_count = int(
                round(class_counts.loc[class_counts[target_class] == 1, 'identifier'].iloc[0]))

            # Create our scaling variable to balance the training process in xgboost algorithm
            pos_scaler = int(round(major_class_count / minor_class_count))
            cv_groups = ['period_dt', 'class_model_id', target_class]
            # cv_groups = cv_group_vars

            if 'xgboost' in model_list:

                pid = 'xgboost_' + str(period.strftime('%Y_%m_%d')) + '_' + target_class
                hyperparams_periods_sql = text(
                    f"""
                    SELECT distinct period_dt 
                    FROM model_class_hyperparameters 
                    WHERE period_dt < '{period}'
                    AND period_dt > '{min_hyperparam_date}' 
                    AND model_type = 'xgboost'
                    AND version = '{version}'
                    AND target_variable = '{target_class}'
                    """
                )
                hyperparams_periods = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])

                if len(hyperparams_periods) == 0:
                    logging.info('Tuning xgboost classifier hyperparameters')
                    xgboost_tune = xgboost_classifier_tune(data=train_data,
                                                           y_var=target_class,
                                                           x_vars=nonlinear_features,
                                                           weight_var=weight_column,
                                                           scale_pos_weight=pos_scaler,
                                                           cv_column=None,
                                                           cv_groups=cv_groups,
                                                           cv_cuts=3,
                                                           metric='auc',
                                                           n_trials=n_trials,
                                                           n_jobs=n_jobs,
                                                           refit=True,
                                                           early_stopping_rounds=50,
                                                           seed=seed)

                    best_model = xgboost_tune['best_model']

                    # Store trials data
                    xgboost_trials = xgboost_tune['trials']
                    xgboost_trials['period_dt'] = period
                    xgboost_trials['model_id'] = model_id
                    xgboost_trials['model_type'] = 'xgboost'
                    xgboost_trials['version'] = version
                    xgboost_trials['target_variable'] = target_class
                    combined_trials = pd.concat([combined_trials, xgboost_trials], ignore_index=True)

                    # Store hyperparams
                    xgboost_params = xgboost_tune["best_params"]
                    xgboost_params_df = pd.DataFrame({'period_dt': period,
                                                      'model_id': model_id,
                                                      'model_type': 'xgboost',
                                                      'version': version,
                                                      'target_variable': target_class,
                                                      'hyperparameter': xgboost_params.keys(),
                                                      'estimate': xgboost_params.values()})
                    combined_hyperparams = pd.concat([combined_hyperparams, xgboost_params_df],
                                                     ignore_index=True)

                else:
                    logging.info('Using previously turned xgboost classifier hyperparameters')
                    # Grab best hyperparams available
                    hyperparams_periods_max = max(hyperparams_periods)
                    hyperparams_periods_sql = text(
                        f"""
                        SELECT * 
                        FROM model_class_hyperparameters 
                        WHERE period_dt < '{period}'
                        AND period_dt > '{min_hyperparam_date}' 
                        AND model_type = 'xgboost'
                        AND version = '{version}'
                        AND target_variable = '{target_class}'
                        """
                    )
                    xgboost_params_df = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])
                    xgboost_params = dict(zip(xgboost_params_df['hyperparameter'],
                                              xgboost_params_df['estimate']))
                    # List of keys to convert to integers
                    keys_to_convert = ['n_estimators', 'max_depth', 'min_child_weight', 'scale_pos_weight',
                                       'random_state']
                    # Convert specified keys to integers
                    xgboost_params.update({key: int(xgboost_params[key]) for key in keys_to_convert})

                    best_model = XGBClassifier(
                        **xgboost_params,
                        n_jobs=n_jobs,
                        use_label_encoder=False)
                    best_model.fit(
                        train_data.loc[:, temp_features],
                        train_data.loc[:, target_class],
                        sample_weight=train_data.loc[:, weight_column] if weight_column else None
                    )

                # Add best model to output
                combined_models.update({pid: best_model})

                # Training prob and model fit
                train_pred = best_model.predict_proba(train_data.loc[:, temp_features])[:, 1]
                auc = roc_auc_score(train_data.loc[:, target_class], train_pred)
                logloss = log_loss(train_data.loc[:, target_class], train_pred)

                logging.info('AUC for %s: %s', pid, auc)
                logging.info('logloss for %s: %s', pid, logloss)

                # Make prediction with the best estimator that was refit across our data
                prob_temp = pd.Series(
                    best_model.predict_proba(prediction_data.loc[:, temp_features])[:, 1], name='positive_probability')
                prob_temp = pd.concat([prediction_data.loc[:, identifier_vars + ['class_model_id']],
                                       prob_temp], axis=1)
                prob_temp.loc[:, 'model_id'] = model_id
                prob_temp.loc[:, 'version'] = version
                prob_temp.loc[:, 'model_type'] = 'xgboost'
                prob_temp.loc[:, 'target_variable'] = target_class
                combined_prob = pd.concat([combined_prob, prob_temp])

                # Model summary
                temp_summary = pd.DataFrame({
                    'name': pid,
                    'model_type': 'xgboost',
                    'model_class': 'tree',
                    'version': version,
                    'class_model_id': model_id,
                    'target_variable': target_class,
                    'holdout_date': period,
                    'training_start': training_start,
                    'training_end': training_end,
                    'create_time': datetime.now(),
                    'final_logloss': logloss,
                    'final_auc': auc,
                }, index=[0])
                combined_summary = pd.concat([combined_summary, temp_summary], ignore_index=True)

            if 'lightgbm' in model_list:

                pid = 'lightgbm_' + str(period.strftime('%Y_%m_%d')) + '_' + target_class
                hyperparams_periods_sql = text(
                    f"""
                    SELECT distinct period_dt 
                    FROM model_class_hyperparameters 
                    WHERE period_dt < '{period}'
                    AND period_dt > '{min_hyperparam_date}' 
                    AND model_type = 'lightgbm'
                    AND version = '{version}'
                    AND target_variable = '{target_class}'
                    """
                )
                hyperparams_periods = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])

                temp_features = feature_dictionary[version_suffix]['nonlinear']

                if len(hyperparams_periods) == 0:
                    logging.info('Tuning lightgbm classifier hyperparameters')
                    lightgbm_tune = lightgbm_classifier_tune(data=train_data,
                                                             y_var=target_class,
                                                             x_vars=temp_features,
                                                             weight_var=weight_column,
                                                             cv_column=None,
                                                             cv_groups=cv_groups,
                                                             cv_cuts=3,
                                                             metric='auc',
                                                             n_trials=n_trials,
                                                             n_jobs=n_jobs,
                                                             refit=True,
                                                             early_stopping_rounds=50,
                                                             seed=seed)

                    best_model = lightgbm_tune['best_model']

                    # Store trials data
                    lightgbm_trials = lightgbm_tune['trials']
                    lightgbm_trials['period_dt'] = period
                    lightgbm_trials['model_id'] = model_id
                    lightgbm_trials['model_type'] = 'lightgbm'
                    lightgbm_trials['version'] = version
                    lightgbm_trials['target_variable'] = target_class
                    combined_trials = pd.concat([combined_trials, lightgbm_trials], ignore_index=True)

                    # Store hyperparams
                    lightgbm_params = lightgbm_tune["best_params"]
                    lightgbm_params_df = pd.DataFrame({'period_dt': period,
                                                       'model_id': model_id,
                                                       'model_type': 'lightgbm',
                                                       'version': version,
                                                       'target_variable': target_class,
                                                       'hyperparameter': lightgbm_params.keys(),
                                                       'estimate': lightgbm_params.values()})
                    combined_hyperparams = pd.concat([combined_hyperparams, lightgbm_params_df],
                                                     ignore_index=True)

                else:
                    logging.info('Using previously turned lightgbm classifier hyperparameters')
                    # Grab best hyperparams available
                    hyperparams_periods_max = max(hyperparams_periods)
                    hyperparams_periods_sql = text(
                        f"""
                        SELECT * 
                        FROM model_class_hyperparameters 
                        WHERE period_dt < '{period}'
                        AND period_dt > '{min_hyperparam_date}' 
                        AND model_type = 'lightgbm'
                        AND version = '{version}'
                        AND target_variable = '{target_class}'
                        """
                    )
                    lightgbm_params_df = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])
                    lightgbm_params = dict(zip(lightgbm_params_df['hyperparameter'],
                                               lightgbm_params_df['estimate']))
                    # List of keys to convert to integers
                    keys_to_convert = ['n_estimators', 'num_leaves', 'bagging_freq', 'min_data_in_leaf', 'random_state']
                    # Convert specified keys to integers
                    lightgbm_params.update({key: int(lightgbm_params[key]) for key in keys_to_convert})

                    best_model = LGBMClassifier(
                        **lightgbm_params,
                        n_jobs=n_jobs,
                        verbose=-1,
                        class_weight='balanced',
                        objective='binary')
                    best_model.fit(
                        train_data.loc[:, temp_features],
                        train_data.loc[:, target_class],
                        sample_weight=train_data.loc[:, weight_column] if weight_column else None
                    )

                # Add best model to output
                combined_models.update({pid: best_model})

                # Training prob and model fit
                train_pred = best_model.predict_proba(train_data.loc[:, temp_features])[:, 1]
                auc = roc_auc_score(train_data.loc[:, target_class], train_pred)
                logloss = log_loss(train_data.loc[:, target_class], train_pred)

                logging.info('AUC for %s: %s', pid, auc)
                logging.info('logloss for %s: %s', pid, logloss)

                # Make prediction with the best estimator that was refit across our data
                prob_temp = pd.Series(
                    best_model.predict_proba(prediction_data.loc[:, temp_features])[:, 1], name='positive_probability')
                prob_temp = pd.concat([prediction_data.loc[:, identifier_vars + ['class_model_id']],
                                       prob_temp], axis=1)
                prob_temp.loc[:, 'model_id'] = model_id
                prob_temp.loc[:, 'version'] = version
                prob_temp.loc[:, 'model_type'] = 'lightgbm'
                prob_temp.loc[:, 'target_variable'] = target_class
                combined_prob = pd.concat([combined_prob, prob_temp])

                # Model summary
                temp_summary = pd.DataFrame({
                    'name': pid,
                    'model_type': 'lightgbm',
                    'model_class': 'tree',
                    'version': version,
                    'class_model_id': model_id,
                    'target_variable': target_class,
                    'holdout_date': period,
                    'training_start': training_start,
                    'training_end': training_end,
                    'create_time': datetime.now(),
                    'final_logloss': logloss,
                    'final_auc': auc,
                }, index=[0])
                combined_summary = pd.concat([combined_summary, temp_summary], ignore_index=True)

            if 'ridge' in model_list:

                pid = 'ridge_' + str(period.strftime('%Y_%m_%d')) + '_' + target_class
                hyperparams_periods_sql = text(
                    f"""
                    SELECT distinct period_dt 
                    FROM model_class_hyperparameters 
                    WHERE period_dt < '{period}'
                    AND period_dt > '{min_hyperparam_date}' 
                    AND model_type = 'ridge'
                    AND version = '{version}'
                    AND target_variable = '{target_class}'
                    """
                )
                hyperparams_periods = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])

                temp_features = feature_dictionary[version_suffix]['linear']

                if len(hyperparams_periods) == 0:
                    logging.info('Tuning ridge classifier hyperparameters')
                    ridge_tune = logistic_classifier_tune(
                        data=train_data,
                        y_var=target_class,
                        x_vars=temp_features,
                        penalty='l2',
                        weight_var=weight_column,
                        cv_column=None,
                        cv_groups=cv_groups,
                        cv_cuts=3,
                        metric='logloss',  # Default metric for Ridge classifier
                        lambda_values=lambda_values,
                        n_jobs=n_jobs,
                        refit=True,
                        seed=seed
                    )

                    best_model = ridge_tune['best_model']

                    # Store hyperparams
                    ridge_params = ridge_tune["best_params"]
                    ridge_params_df = pd.DataFrame({
                        'period_dt': [period],
                        'model_id': [model_id],
                        'model_type': ['ridge'],
                        'version': [version],
                        'target_variable': [target_class],
                        'hyperparameter': ridge_params.keys(),
                        'estimate': ridge_params.values()
                    })
                    combined_hyperparams = pd.concat([combined_hyperparams, ridge_params_df], ignore_index=True)

                else:
                    logging.info('Using previously tuned ridge classifier hyperparameters')
                    # Grab best hyperparams available
                    hyperparams_periods_max = max(hyperparams_periods)
                    hyperparams_periods_sql = text(
                        f"""
                        SELECT * 
                        FROM model_class_hyperparameters 
                        WHERE period_dt < '{period}'
                        AND period_dt > '{min_hyperparam_date}' 
                        AND model_type = 'ridge'
                        AND version = '{version}'
                        AND target_variable = '{target_class}'
                        """
                    )
                    ridge_params_df = pd.read_sql(hyperparams_periods_sql, sql_connections_dev['models'])
                    ridge_params = dict(zip(ridge_params_df['hyperparameter'], ridge_params_df['estimate']))

                    best_model = LogisticRegression(
                        penalty='l2',
                        fit_intercept=True,
                        solver='saga',
                        random_state=seed,
                        max_iter=1000,
                        class_weight='balanced',
                        **ridge_params
                    )
                    best_model.fit(
                        train_data.loc[:, temp_features],
                        train_data.loc[:, target_class],
                        sample_weight=train_data.loc[:, weight_column] if weight_column else None
                    )

                # Add best model to output
                combined_models.update({pid: best_model})

                # Training predictions and model evaluation
                train_pred = best_model.decision_function(train_data.loc[:, temp_features])
                auc = roc_auc_score(train_data.loc[:, target_class], train_pred)
                logloss = log_loss(train_data.loc[:, target_class],
                                   best_model.predict_proba(train_data.loc[:, temp_features]))

                logging.info('AUC for %s: %s', pid, auc)
                logging.info('logloss for %s: %s', pid, logloss)

                # Make prediction with the best estimator that was refit across our data
                prob_temp = pd.Series(
                    best_model.predict_proba(prediction_data.loc[:, temp_features])[:, 1], name='positive_probability')
                prob_temp = pd.concat([prediction_data.loc[:, identifier_vars + ['class_model_id']],
                                       prob_temp], axis=1)
                prob_temp.loc[:, 'model_id'] = model_id
                prob_temp.loc[:, 'version'] = version
                prob_temp.loc[:, 'model_type'] = 'ridge'
                prob_temp.loc[:, 'target_variable'] = target_class
                combined_prob = pd.concat([combined_prob, prob_temp])

                # Model summary
                temp_summary = pd.DataFrame({
                    'name': pid,
                    'model_type': 'ridge',
                    'model_class': 'logistic',
                    'version': version,
                    'class_model_id': model_id,
                    'target_variable': target_class,
                    'holdout_date': period,
                    'training_start': training_start,
                    'training_end': training_end,
                    'create_time': datetime.now(),
                    'final_logloss': logloss,
                    'final_auc': auc,
                }, index=[0])
                combined_summary = pd.concat([combined_summary, temp_summary], ignore_index=True)

        if compute_performance:
            combined_prob_w = combined_prob.pivot(index=identifier_vars + ['model_id', 'version', 'model_type'],
                                                  columns='target_variable',
                                                  values='positive_probability').reset_index()
            combined_signal_data = pd.merge(combined_prob_w,
                                            prediction_data,
                                            how='inner',
                                            on=identifier_vars)

            # Standardize new signal expected returns
            combined_signal_data = standardize(data=combined_signal_data,
                                               columns=['long', 'short'],
                                               normal=False,
                                               groups=['period_dt', 'portfolio_id', 'model_type'],
                                               suffix='_std',
                                               missing_strategy='median')

            # Create net
            combined_signal_data['net_std'] = combined_signal_data['long_std'] - combined_signal_data['short_std']

            combined_signal_data = standardize(data=combined_signal_data,
                                               columns=['net_std'],
                                               normal=False,
                                               groups=['period_dt', 'portfolio_id', 'model_type'],
                                               suffix=None,
                                               missing_strategy='median')

            # Create new neutralized signal without weighting
            combined_signal_data = neutralize(data=combined_signal_data,
                                              target_vars=['net_std'],
                                              neutralization_vars=risk_factors_exp,
                                              groups=['period_dt', 'portfolio_id', 'model_type'],
                                              intercept=True,
                                              missing_strategy=None,
                                              residual_suffix='_rfres')

            # Standardize out final signal variable
            combined_signal_data = standardize(data=combined_signal_data,
                                               columns='net_std_rfres',
                                               normal=False,
                                               groups=['period_dt', 'portfolio_id', 'model_type'],
                                               suffix='_std',
                                               missing_strategy=None)

            # Drop risk factors from our data
            combined_signal_data.drop(columns=risk_factors, inplace=True)

            # Compute payoffs
            signal_payoffs = linear_regression(data=combined_signal_data,
                                               y_var='totalreturn_usd_win_rfres',
                                               x_vars='net_std_rfres_std',
                                               groups=['period_dt', 'model_id', 'model_type', 'portfolio_id',
                                                       'version'],
                                               intercept=True,
                                               weight=riskweight_var
                                               )

            signal_payoffs_est = signal_payoffs['coeffs'].rename(columns={'estimate_id': 'factor_id',
                                                                          'estimate': 'payoff'})

            signal_ic = correlation(data=combined_signal_data,
                                    corr_vars='net_std_rfres_std',
                                    with_vars='totalreturn_usd_win_rfres',
                                    groups=['period_dt', 'model_id', 'model_type', 'portfolio_id',
                                            'version'],
                                    weight=riskweight_var
                                    )

            signal_quantiles_data = compute_quantile_returns(data=combined_signal_data,
                                                             quantile_variable='net_std_rfres_std',
                                                             n_quantiles=5,
                                                             return_variables='totalreturn_usd_win_rfres',
                                                             groups=['period_dt', 'model_id', 'model_type',
                                                                     'portfolio_id', 'version'],
                                                             weight=None,
                                                             compute_spread=True)

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=signal_quantiles_data['quantiles'],
                            table_name='class_signal_quantile_returns',
                            keys=['period_dt', 'model_id', 'model_type', 'version',
                                  'portfolio_id', 'quantile_group_id', 'quantile_variable',
                                  'return_variable', 'weight'],
                            create_if_missing=True,
                            var_char_max=40
                            )

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=signal_quantiles_data['spreads'],
                            table_name='class_signal_quantile_spread_returns',
                            keys=['period_dt', 'model_id', 'model_type', 'version',
                                  'portfolio_id', 'quantile_variable', 'return_variable', 'weight'],
                            create_if_missing=True,
                            var_char_max=40
                            )

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=signal_ic,
                            table_name='class_signal_ic',
                            keys=['period_dt', 'model_id', 'model_type', 'version',
                                  'portfolio_id', 'corr_var', 'with_var', 'corr_type'],
                            create_if_missing=True,
                            var_char_max=40
                            )

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=signal_payoffs_est,
                            table_name='class_signal_payoffs',
                            keys=['period_dt', 'model_id', 'model_type', 'portfolio_id',
                                  'version', 'factor_id'],
                            create_if_missing=True,
                            var_char_max=40
                            )

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=signal_payoffs['summary'],
                            table_name='class_signal_pay_reg_summary',
                            keys=['period_dt', 'model_id', 'model_type', 'portfolio_id',
                                  'version'],
                            create_if_missing=True,
                            var_char_max=40
                            )

            # Output signal data
            combined_signal_data = combined_signal_data.loc[:, ['period_dt', 'identifier', 'exchange_country',
                                                                'model_id', 'version', 'model_type', 'portfolio_id',
                                                                riskweight_var,
                                                                'totalreturn_usd_win_rfres', 'net_std',
                                                                'net_std_rfres_std'
                                                                ]].rename(
                columns={riskweight_var: 'riskweight'}).reset_index(drop=True)

            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=combined_signal_data,
                            table_name='class_signal_data',
                            keys=['period_dt', 'model_id', 'model_type', 'portfolio_id',
                                  'version', 'identifier', 'exchange_country'],
                            create_if_missing=True,
                            var_char_max=40
                            )

        # Format/clean and export prediction data, feature impacts, best hyperparameters, and other metadata

        upsert_table_sa(engine=sql_engines['models'],
                        data_frame=combined_summary,
                        table_name='model_class_logs',
                        keys=['holdout_date', 'class_model_id', 'target_variable', 'model_type', 'version'],
                        create_if_missing=True)

        upsert_table_sa(engine=sql_engines['models'],
                        data_frame=combined_prob,
                        table_name='model_class_prob',
                        keys=['period_dt',
                              'class_model_id',
                              'target_variable',
                              'model_type',
                              'version',
                              'identifier',
                              'exchange_country'],
                        create_if_missing=True)

        if run_feature_impact:
             upsert_table_sa(engine=sql_engines['models'],
                                    data_frame=class_model_dict['feature_impact'],
                                    table_name='model_class_feat_imp',
                                    keys=['period_dt',
                                          'class_model_id',
                                          'universe',
                                          'target_variable',
                                          'model_type',
                                          'version',
                                          'feature_name'],
                                    create_if_missing=True)

        # Pivot to make dataset easier to look at the net coefficient effects

         projs_df_coeff = class_model_dict[
             'coefficient_estimates'].pivot(index=['period_dt', 'class_model_id', 'model_type', 'version',
                                                   'universe', 'model_feature'],
                                            columns=['target_variable'],
                                        values=['estimate']).droplevel(level=0, axis=1).reset_index()
        
         upsert_table_sa(engine=sql_engines['models'],
                                data_frame=projs_df_coeff,
                                table_name='model_class_coeff',
                                keys=['period_dt',
                                      'class_model_id',
                                      'model_type',
                                      'version',
                                      'universe',
                                      'model_feature'],
                                create_if_missing=True)

        if len(combined_hyperparams) > 0:
            upsert_table_sa(engine=sql_engines['models'],
                            data_frame=combined_hyperparams,
                            table_name='model_class_hyperparameters',
                            keys=['period_dt',
                                  'model_id',
                                  'target_variable',
                                  'model_type',
                                  'version',
                                  'hyperparameter'],
                            create_if_missing=True)

        if len(combined_trials) > 0:
            upsert_table_sa(engine=sql_engines['research'],
                            data_frame=combined_trials,
                            table_name='class_signal_study_trials',
                            keys=['period_dt',
                                  'model_id',
                                  'model_type',
                                  'version',
                                  'metric',
                                  'trial',
                                  'target_variable'],
                            create_if_missing=True)

        # Output models in serialized format via pickle so they can be referenced later

        period_label = str(period.strftime('%Y_%m_%d'))
        model_pickle_file = open(class_path + '\\models\\' + model_id + '_models_' '_' +
                                 period_label + '.pickle', "wb")
        pickle.dump(combined_models, model_pickle_file)
        model_pickle_file.close()

workflow_set_state(engine=sql_engines['investment'],
                    context='equity.model.classification.py_estimate',
                    state_id='completed',
                    state='1')
