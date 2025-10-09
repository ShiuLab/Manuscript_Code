"""
XGBoost Regression

# Required Inputs
    -X      Path to feature matrix file
    -y_name Column name of label in X matrix
    -test   File containing list of test instances
    -save   Path to save output files
    -prefix Prefix of output file names
    
    # Optional
    -Y      Path to label matrix file, if label not in X matrix
    -tag    Feature types/identifier for output file naming
    -fold   k folds for Cross-Validation (default is 5)
    -n      Number of CV repetitions (default is 10)
    -feat   File containing features (from X) to include in model
    -plot   Plot feature importances and predictions (default is t)

# Outputs for each training repetition (prefixed with <prefix>_)
    _lm_test_rep_*.pdf        Regression plot of predicted and actual test labels
    _model_rep_*.save         XGBoost model
    _top20_rep_*.pdf          Plot of top 20 features' importance scores

# Summary outputs (prefixed with <prefix>_)
    _imp.csv                  Feature importance scores
    _cv_results.csv           Cross-validation results (various metrics)
    _test_results.csv         Evaluation results (various metrics)
    RESULTS_xgboost.txt       Aggregated results (various metrics)
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
from sklearn.model_selection import KFold, cross_validate, StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import KBinsDiscretizer


def hyperopt_objective(params, X_train, y_train):
    """
    Create the hyperparameter grid and run Hyperopt hyperparameter tuning
    with K-fold cross-validation
    Written by Thejesh Mallidi
    Modified by Kenia Segura Abá
    """
    reg = xgb.XGBRegressor(
        learning_rate=params["learning_rate"],
        max_depth=int(params["max_depth"]),
        subsample=params["subsample"],
        colsample_bytree=params["colsample_bytree"],
        gamma=params["gamma"],
        # I excluded alpha when I ran the 20240531 dataset models
        alpha=params["alpha"],
        min_child_weight=params["min_child_weight"],
        n_estimators=int(params["n_estimators"]),
        objective=params["objective"],
        eval_metric=params["eval_metric"],
        random_state=42
    )

    # Bin y_train for stratified K-fold sampling
    discretizer = KBinsDiscretizer(
        n_bins=10, encode="ordinal", strategy="quantile")
    y_train_binned = discretizer.fit_transform(
        y_train.values.reshape(-1, 1)).astype(int).flatten()
    y_train_binned = pd.Series(y_train_binned, index=y_train.index)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    validation_loss = []
    for i, (train_idx, val_idx) in enumerate(cv.split(X_train, y_train_binned)):
        X_train_cv, X_val_cv = X_train.iloc[train_idx], X_train.iloc[val_idx]
        y_train_cv, y_val_cv = y_train.iloc[train_idx], y_train.iloc[val_idx]

        reg.fit(X_train_cv, y_train_cv)
        cv_pred = reg.predict(X_val_cv)
        validation_loss.append(r2_score(y_val_cv, cv_pred))

    # Note: Hyperopt minimizes the objective, so to maximize R^2, return the negative mean
    return -np.mean(validation_loss)


def param_hyperopt(param_grid, X_train, y_train, max_evals=100):
    """
    Obtain the best parameters from Hyperopt
    Written by Thejesh Mallidi
    """
    trials = Trials()
    params_best = fmin(
        fn=lambda params: hyperopt_objective(params, X_train, y_train),
        space=param_grid,
        algo=tpe.suggest,
        max_evals=max_evals,
        trials=trials
    )

    print("\n\nBest parameters:", params_best)
    return params_best, trials


def xgb_reg(trait, X, y, fold, n, prefix, plot):
    """ Train XGBoost Regression Model """
    print(trait)

    # Train-test split
    X_train = X.loc[~X.index.isin(test[0])]
    X_test = X.loc[X.index.isin(test[0])]
    y_train = y.loc[~y.index.isin(test[0])]
    y_test = y.loc[y.index.isin(test[0])]

    # Ensure rows are in the same order
    X_train = X_train.loc[y_train.index, :]
    X_test = X_test.loc[y_test.index, :]

    # Hyperparameter tuning
    parameters = {"learning_rate": hp.uniform("learning_rate", 0.01, 0.4),  # learning rate
                  # tree depth
                  "max_depth": scope.int(hp.quniform("max_depth", 2, 10, 1)),
                  # instances per tree
                  "subsample": hp.uniform("subsample", 0.5, 1.0),
                  # features per tree
                  "colsample_bytree": hp.uniform("colsample_bytree", 0.7, 1.0),
                  "gamma": hp.uniform("gamma", 0.1, 5.0),  # min_split_loss
                  "alpha": hp.uniform("alpha", 0.1, 5.0),  # L1 regularization
                  # minimum sum of instance weight needed in a child
                  "min_child_weight": scope.int(hp.quniform("min_child_weight", 5, 20, 2)),
                  "n_estimators": scope.int(hp.quniform("n_estimators", 5, 500, 5)),
                  "objective": "reg:squarederror", "eval_metric": "rmse"}

    start = time.time()
    best_params, trials = param_hyperopt(parameters, X_train, y_train, 100)
    run_time = time.time() - start
    print("Total hyperparameter tuning time:", run_time)
    print("Best parameters: ", best_params)
    print("Trials", trials)

    ################## Training with Cross-Validation ##################
    results_cv = []  # hold performance metrics of cv reps
    results_test = []  # hold performance metrics on test set
    feature_imp = pd.DataFrame(index=X_train.columns)
    Y_preds = pd.DataFrame(y.copy(deep=True))

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
            alpha=best_params["alpha"],
            min_child_weight=int(best_params["min_child_weight"]),
            n_estimators=int(best_params["n_estimators"]),
            objective="reg:squarederror",
            eval_metric="rmse",
            random_state=j)

        # cv_pred = cross_val_predict(
        #     best_model, X_train, y_train, cv=fold, n_jobs=-1) # predictions

        # Bin y_train for stratified K-fold sampling
        discretizer = KBinsDiscretizer(
            n_bins=10, encode="ordinal", strategy="quantile")
        y_train_binned = discretizer.fit_transform(
            y_train.values.reshape(-1, 1)).astype(int).flatten()
        y_train_binned = pd.Series(y_train_binned, index=y_train.index)

        cv = StratifiedKFold(n_splits=fold, shuffle=True, random_state=42)
        Y_preds[f"cv_preds_{j}"] = np.nan
        cv_pred = []
        for i, (train_idx, val_idx) in enumerate(cv.split(X_train, y_train_binned)):
            X_train_cv, X_val_cv = X_train.iloc[train_idx], X_train.iloc[val_idx]
            y_train_cv, y_val_cv = y_train.iloc[train_idx], y_train.iloc[val_idx]

            best_model.fit(X_train_cv, y_train_cv)
            cv_preds = best_model.predict(X_val_cv)
            Y_preds.loc[y_train.index[val_idx], f"cv_preds_{j}"] = cv_preds
            cv_pred.append(r2_score(y_val_cv, cv_preds))

        # Performance statistics on validation set
        mse_val = mean_squared_error(
            y_train, Y_preds.loc[y_train.index, f"cv_preds_{j}"])
        rmse_val = np.sqrt(mean_squared_error(
            y_train, Y_preds.loc[y_train.index, f"cv_preds_{j}"]))
        evs_val = explained_variance_score(
            y_train, Y_preds.loc[y_train.index, f"cv_preds_{j}"])
        r2_val = r2_score(y_train, Y_preds.loc[y_train.index, f"cv_preds_{j}"])
        cor_val = np.corrcoef(
            np.array(y_train), Y_preds.loc[y_train.index, f"cv_preds_{j}"])
        print("Val MSE: %f" % (mse_val))
        print("Val RMSE: %f" % (rmse_val))
        print("Val R-sq: %f" % (r2_val))
        print("Val PCC: %f" % (cor_val[0, 1]))
        result_val = [mse_val, rmse_val, evs_val, r2_val, cor_val[0, 1]]
        results_cv.append(result_val)

        # Evaluate the model on the test set
        best_model.fit(X_train, y_train)
        y_pred = best_model.predict(X_test)
        Y_preds[f"test_preds_{j}"] = np.nan
        Y_preds.loc[y_test.index, f"test_preds_{j}"] = y_pred

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
        feature_imp = pd.concat([feature_imp, pd.Series(best_model.feature_importances_,
                                                        index=best_model.feature_names_in_, name=f"rep_{j}")],
                                ignore_index=False, axis=1)

        # Save predicted labels to file
        # preds[f"rep_{j}"] = pd.concat([pd.Series(cv_pred, index=X_train.index),
        #     pd.Series(y_pred, index=X_test.index)], axis=0)

        if plot == "t":
            # Plot linear regression of actual and predicted test values
            plt.figure(figsize=(3, 3))
            plt.axline([0, 0], [1, 1], color="black", linestyle="--")
            plt.scatter(x=y_pred, y=y_test, color="blue",
                        alpha=0.7, marker=".")
            plt.xlim(np.min(y_pred)-.1, np.max(y_pred)+.1)
            plt.ylim(np.min(y_pred)-.1, np.max(y_pred)+.1)
            plt.xlabel("Predicted")
            plt.ylabel("Actual")
            plt.title(f"{trait}")
            plt.axis("square")
            plt.tight_layout()
            plt.savefig(
                f"{args.save}/{prefix}_lm_test_rep_{j}.pdf", format="pdf")
            plt.close()

            # Plot feature importances
            xgb.plot_importance(
                best_model, grid=False, max_num_features=20,
                title=f"{trait} Feature Importances", xlabel="Weight")
            plt.tight_layout()
            plt.savefig(
                f"{args.save}/{prefix}_top20_rep_{j}.pdf", format="pdf")
            plt.close()

    # Write feature importances across reps to file
    feature_imp.to_csv(f"{args.save}/{prefix}_imp.csv")

    # Write predictions across reps to file
    # pd.DataFrame.from_dict(preds).to_csv(f"{args.save}/{prefix}_preds.csv")
    Y_preds.to_csv(f"{args.save}/{prefix}_cv_preds.csv")

    return (results_cv, results_test)


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(
        description="XGBoost Regression on SNP and ORF data")

    # Required input
    req_group = parser.add_argument_group(title="Required Input")
    req_group.add_argument(
        "-X", help="path to feature table file", required=True)
    req_group.add_argument(
        "-y_name", help="name of label in X file", required=True)
    req_group.add_argument(
        "-test", help="path to file of test set instances", required=True)
    req_group.add_argument(
        "-save", help="path to save output files", required=True)
    req_group.add_argument(
        "-prefix", help="prefix of output file names", required=True)

    # Optional input
    req_group.add_argument(
        "-Y", help="path to label table file", default="")
    req_group.add_argument(
        "-tag", help="description about run to add to results file", default="")
    req_group.add_argument(
        "-fold", help="k number of cross-validation folds", default=5)
    req_group.add_argument(
        "-n", help="number of cross-validation repetitions", default=10)
    req_group.add_argument(
        "-feat", help="file containing features (from X) to include in model", default="all")
    req_group.add_argument(
        "-feat_list", help="comma-separated list of features (from X) to include in model",
        nargs="+", type=lambda s: [col.strip() for col in s.split(",")], default=[])
    req_group.add_argument(
        "-plot", help="plot feature importances and predictions (t/f)", default="t")

    # Help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()  # Read arguments

    # Read in data
    X = dt.fread(args.X)  # feature table
    X = X.to_pandas()
    X.set_index(X.columns[0], inplace=True)

    if args.Y == "":  # get the label from X or Y files
        y = X.loc[:, args.y_name]
        X.drop(columns=args.y_name, inplace=True)
    else:
        Y = dt.fread(args.Y).to_pandas()
        Y.set_index(Y.columns[0], inplace=True)
        y = Y.loc[:, args.y_name]

    test = pd.read_csv(args.test, header=None)  # test instances

    # Filter out features not in the given feat file - default: keep all
    if args.feat != "all":
        print("Using subset of features from: %s" % args.feat)
        with open(args.feat) as f:
            features = f.read().strip().splitlines()
        X = X.loc[:, features]
        print(f"New dimensions: {X.shape}")

    if len(args.feat_list) > 0:
        print("Using subset of features from list", args.feat_list[0])
        X = X.loc[:, args.feat_list[0]]
        print(f"New dimensions: {X.shape}")

    # Train the model
    start = time.time()
    results_cv, results_test = xgb_reg(
        args.y_name, X, y, int(args.fold), int(args.n), args.prefix, args.plot)
    run_time = time.time() - start
    print("Training Run Time: %f" % (run_time))

    # Save results to file
    results_cv = pd.DataFrame(
        results_cv,
        columns=["MSE_val", "RMSE_val", "EVS_val", "R2_val", "PCC_val"])
    results_test = pd.DataFrame(
        results_test,
        columns=["MSE_test", "RMSE_test", "EVS_test", "R2_test", "PCC_test"])

    # Aggregate results and save to file
    if not os.path.isfile(f"{args.save}/RESULTS_xgboost.txt"):
        out = open(f"{args.save}/RESULTS_xgboost.tsv", "a")
        out.write("Date\tRunTime\tTag\tY\tNumInstances\tNumFeatures")
        out.write("\tCV_fold\tCV_rep\tMSE_val\tMSE_val_sd")
        out.write("\tRMSE_val\tRMSE_val_sd\tEVS_val\tEVS_val_sd\tR2_val")
        out.write("\tR2_val_sd\tPCC_val\tPCC_val_sd\tMSE_test\tMSE_test_sd")
        out.write("\tRMSE_test\tRMSE_test_sd\tEVS_test\tEVS_test_sd\tR2_test")
        out.write("\tR2_test_sd\tPCC_test\tPCC_test_sd")
        out.close()

    out = open(f"{args.save}/RESULTS_xgboost.tsv", "a")
    out.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\t")
    out.write(
        f"{run_time}\t{args.tag}\t{args.y_name}\t{X.shape[0]-len(test)}\t")
    out.write(f"{X.shape[1]}\t{int(args.fold)}\t{int(args.n)}\t")
    out.write(f"{np.mean(results_cv.MSE_val)}\t{np.std(results_cv.MSE_val)}\t")
    out.write(f"{np.mean(results_cv.RMSE_val)}\t{np.std(results_cv.RMSE_val)}\t")
    out.write(f"{np.mean(results_cv.EVS_val)}\t{np.std(results_cv.EVS_val)}\t")
    out.write(f"{np.mean(results_cv.R2_val)}\t{np.std(results_cv.R2_val)}\t")
    out.write(f"{np.mean(results_cv.PCC_val)}\t{np.std(results_cv.PCC_val)}\t")
    out.write(
        f"{np.mean(results_test.MSE_test)}\t{np.std(results_test.MSE_test)}\t")
    out.write(
        f"{np.mean(results_test.RMSE_test)}\t{np.std(results_test.RMSE_test)}\t")
    out.write(
        f"{np.mean(results_test.EVS_test)}\t{np.std(results_test.EVS_test)}\t")
    out.write(
        f"{np.mean(results_test.R2_test)}\t{np.std(results_test.R2_test)}\t")
    out.write(
        f"{np.mean(results_test.PCC_test)}\t{np.std(results_test.PCC_test)}")
