"""
XGBoost for regression on SNP, ORF, and CNO data from Peter et al. 2018

# Required Inputs
    -X      Path to feature matrix
    -Y      Path to label matrix
    -test   File containing list of test instances
    -trait  Column name of target trait in Y matrix
    -save   Path to save output files
    -prefix Prefix of output file names
    
    # Optional
    -type   Feature types (e.g. SNP (default), ORF, CNO)
    -fold   k folds for Cross-Validation (default is 5)
    -n      Number of CV repetitions (default is 10)
    -feat   File containing features (from X) to include in model
    -plot   Plot feature importances and predictions (default is t)

# Outputs (prefixed with <type>_<trait>)
    _GridSearch.csv     Grid Search CV R-sq scores
    _lm_test.pdf        Regression plot of predicted and actual test labels
    _model.save         XGBoost model
    _imp.csv            Feature importance scores
    _top20.pdf          Plot of top 20 features' importance scores
    _cv_results.csv     Cross-validation results (various metrics)
    _test_results.csv   Evaluation results (various metrics)
    RESULTS_xgboost.txt Aggregated results (various metrics)    
"""
__author__ = "Kenia Segura Abá"

from configparser import ExtendedInterpolation
import sys
import os
import argparse
import time
import random
import pickle
import datatable as dt
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from hyperopt import hp, fmin, tpe, Trials
from hyperopt.pyll.base import scope
from sklearn.metrics import mean_squared_error, r2_score, explained_variance_score
from sklearn.model_selection import KFold, cross_validate
from sklearn.model_selection import cross_val_predict


def hyperopt_objective(params):
    # Written by Thejesh Mallidi
    reg = xgb.XGBRegressor(
        learning_rate=params["learning_rate"],
        max_depth=int(params["max_depth"]),
        subsample=params["subsample"],
        colsample_bytree=params["colsample_bytree"],
        gamma=params["gamma"],
        n_estimators=int(params["n_estimators"]),
        random_state=42
    )

    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    validation_loss = cross_validate(
        reg, X_train, y_train,
        scoring="r2",
        cv=cv,
        n_jobs=-1,
        error_score="raise"
    )

    # Note: Hyperopt minimizes the objective, so to maximize R^2, return the negative mean
    return -np.mean(validation_loss["test_score"])


def param_hyperopt(param_grid, max_evals=100):
    # Written by Thejesh Mallidi
    trials = Trials()
    params_best = fmin(
        fn=hyperopt_objective,
        space=param_grid,
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials
    )

    print("\n\nBest parameters:", params_best)
    return params_best, trials


def xgb_reg(trait, fold, n, data_type, prefix, plot):
    global X_train
    global y_train

    """ Train XGBoost Regression Model """
    print(trait)
    y = Y[trait]

    # Train-test split
    X_train = X.loc[~X.index.isin(test[0])]
    X_test = X.loc[X.index.isin(test[0])]
    y_train = y.loc[~y.index.isin(test[0])]
    y_test = y.loc[y.index.isin(test[0])]

    # Ensure rows are in the same order
    X_train = X_train.loc[y_train.index, :]
    X_test = X_test.loc[y_test.index, :]

    # Hyperparameter tuning
    parameters = {"learning_rate": hp.uniform("learning_rate", 0.01, 0.5),  # learning rate
                  # tree depth
                  "max_depth": scope.int(hp.quniform("max_depth", 3, 10, 1)),
                  # instances per tree
                  "subsample": hp.uniform("subsample", 0.5, 1.0),
                  # features per tree
                  "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1.0),
                  "gamma": hp.uniform("gamma", 0.0, 1.0),  # min_split_loss
                  "n_estimators": scope.int(hp.quniform("n_estimators", 50, 500, 2))
                  }  # sample training instances

    start = time.time()
    best_params, trials = param_hyperopt(parameters, 100)
    run_time = time.time() - start
    print("Total hyperparameter tuning time:", run_time)

    ################## Training with Cross-Validation ##################
    results_cv = []  # hold performance metrics of cv reps
    results_test = []  # hold performance metrics on test set
    feature_imp = pd.DataFrame(columns=[f"rep_{i}" for i in range(0, n)])
    preds = {}

    # Training with Cross-validation
    for j in range(0, n):  # repeat cv 10 times
        print(f"Running {j+1} of {n}")
        # Build model using the best parameters
        best_model = xgb.XGBRegressor(
            eta=best_params["learning_rate"],
            max_depth=int(best_params["max_depth"]),
            subsample=best_params["subsample"],
            colsample_bytree=best_params["colsample_bytree"],
            gamma=best_params["gamma"],
            n_estimators=int(best_params["n_estimators"]),
            random_state=j)

        cv_pred = cross_val_predict(
            best_model, X_train, y_train, cv=fold, n_jobs=-1)  # predictions

        # Performance statistics on validation set
        mse_val = mean_squared_error(y_train, cv_pred)
        rmse_val = np.sqrt(mean_squared_error(y_train, cv_pred))
        evs_val = explained_variance_score(y_train, cv_pred)
        r2_val = r2_score(y_train, cv_pred)
        cor_val = np.corrcoef(np.array(y_train), cv_pred)
        print("Val MSE: %f" % (mse_val))
        print("Val RMSE: %f" % (rmse_val))
        print("Val R-sq: %f" % (r2_val))
        print("Val PCC: %f" % (cor_val[0, 1]))
        result_val = [mse_val, rmse_val, evs_val, r2_val, cor_val[0, 1]]
        results_cv.append(result_val)

        # Evaluate the model on the test set
        best_model.fit(X_train, y_train)
        y_pred = best_model.predict(X_test)

        # Performance on the test set
        mse = mean_squared_error(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        evs = explained_variance_score(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        cor = np.corrcoef(np.array(y_test), y_pred)
        print("Test MSE: %f" % (mse))
        print("Test RMSE: %f" % (rmse))
        print("Test R-sq: %f" % (r2))
        print("Test PCC: %f" % (cor[0, 1]))
        result_test = [mse, rmse, evs, r2, cor[0, 1]]
        results_test.append(result_test)

        # Save the fitted model to a file
        filename = f"{args.save}/{prefix}_model_rep_{j}.pkl"
        pickle.dump(best_model, open(filename, "wb"))

        # Save feature importance scores to file
        feature_imp[f"rep_{j}"] = pd.Series(best_model.feature_importances_)

        # Save predicted labels to file
        preds[f"rep_{j}"] = pd.concat([pd.Series(cv_pred, index=X_train.index),
                                       pd.Series(y_pred, index=X_test.index)], axis=0)

        if plot == "t":
            # Plot linear regression of actual and predicted test values
            sns.regplot(x=y_pred, y=y_test, fit_reg=True,
                        ci=95, seed=123, color="black")
            plt.xlabel("Predicted")
            plt.ylabel("Actual")
            plt.title(f"{trait}")
            plt.savefig(
                f"{args.save}/{prefix}_lm_test_rep_{j}.pdf", format="pdf")
            plt.close()

            # Plot feature importances
            xgb.plot_importance(
                best_model, grid=False, max_num_features=20,
                title=f"{trait} Feature Importances", xlabel="Weight")
            plt.savefig(
                f"{args.save}/{prefix}_top20_rep_{j}.pdf", format="pdf")
            plt.close()

    # Write feature importances across reps to file
    feature_imp.to_csv(f"{args.save}/{prefix}_imp.csv")

    # Write predictions across reps to file
    pd.DataFrame.from_dict(preds).to_csv(f"{args.save}/{prefix}_preds.csv")

    return (results_cv, results_test)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(
        description="XGBoost Regression on SNP and ORF data")

    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", help="path to feature data file", required=True)
    req_group.add_argument(
        "-Y", help="path to label data file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-trait", help="name of trait column in y dataframe", required=True)
    req_group.add_argument(
        "-save", help="path to save output files", required=True)
    req_group.add_argument(
        "-prefix", help="prefix of output file names", required=True)

    # Optional input
    req_group.add_argument(
        "-type", help="data type of X matrix (e.g. SNP, PAV, CNV)", default="")
    req_group.add_argument(
        "-fold", help="k number of cross-validation folds", default=5)
    req_group.add_argument(
        "-n", help="number of cross-validation repetitions", default=10)
    req_group.add_argument(
        "-feat", help="file containing features (from X) to include in model", default="all")
    req_group.add_argument(
        "-plot", help="plot feature importances and predictions (t/f)", default="t")

    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()  # Read arguments

    # Read in data
    X = dt.fread(args.X)
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)
    Y = pd.read_csv(args.Y, index_col=0)
    test = pd.read_csv(args.test, sep="\t", header=None)

    # Filter out features not in feat file given - default: keep all
    if args.feat != "all":
        print("Using subset of features from: %s" % args.feat)
        with open(args.feat) as f:
            features = f.read().strip().splitlines()
        X = X.loc[:, features]
        print(f"New dimensions: {X.shape}")

    # Train the model
    start = time.time()
    results_cv, results_test = xgb_reg(
        args.trait, int(args.fold), int(args.n), args.type, args.prefix, args.plot)
    run_time = time.time() - start
    print("Training Run Time: %f" % (run_time))

    # Save results to file
    results_cv = pd.DataFrame(
        results_cv,
        columns=["MSE_val", "RMSE_val", "EVS_val", "R2_val", "PCC_val"])
    # results_cv.to_csv(f"{args.save}/{args.type}_{args.trait}_cv_results.csv")
    results_test = pd.DataFrame(
        results_test,
        columns=["MSE_test", "RMSE_test", "EVS_test", "R2_test", "PCC_test"])
    # results_test.to_csv(
    #     f"{args.save}/{args.type}_{args.trait}_test_results.csv")

    # Aggregate results and save to file
    if not os.path.isfile(f"{args.save}/RESULTS_xgboost.txt"):
        out = open(f"{args.save}/RESULTS_xgboost.txt", "a")
        out.write("Date\tRunTime\tData\tTrait\tNumInstances\tFeatureNum\
            \tCVfold\tCV_rep\tMSE_val\tMSE_val_sd\
            \tRMSE_val\tRMSE_val_sd\tEVS_val\tEVS_val_sd\tR2_val\
            \tR2_val_sd\tPCC_val\tPCC_val_sd\tMSE_test\tMSE_test_sd\
            \tRMSE_test\tRMSE_test_sd\tEVS_test\tEVS_test_sd\tR2_test\
            \tR2_test_sd\tPCC_test\tPCC_test_sd\n")
        out.close()

    out = open(f"{args.save}/RESULTS_xgboost.txt", "a")
    out.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t\
        {run_time}\t{args.type}\t{args.trait}\t{X.shape[0]-len(test)}\t\
        {X.shape[1]}\t{int(args.fold)}\t{int(args.n)}\t\
        {np.mean(results_cv.MSE_val)}\t{np.std(results_cv.MSE_val)}\t\
        {np.mean(results_cv.RMSE_val)}\t{np.std(results_cv.RMSE_val)}\t\
        {np.mean(results_cv.EVS_val)}\t{np.std(results_cv.EVS_val)}\t\
        {np.mean(results_cv.R2_val)}\t{np.std(results_cv.R2_val)}\t\
        {np.mean(results_cv.PCC_val)}\t{np.std(results_cv.PCC_val)}\t\
        {np.mean(results_test.MSE_test)}\t{np.std(results_test.MSE_test)}\t\
        {np.mean(results_test.RMSE_test)}\t{np.std(results_test.RMSE_test)}\t\
        {np.mean(results_test.EVS_test)}\t{np.std(results_test.EVS_test)}\t\
        {np.mean(results_test.R2_test)}\t{np.std(results_test.R2_test)}\t\
        {np.mean(results_test.PCC_test)}\t{np.std(results_test.PCC_test)}")
