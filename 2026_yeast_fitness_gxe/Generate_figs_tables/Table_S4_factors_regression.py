#!/usr/bin/env python3
# use conda environment: /mnt/home/seguraab/miniconda3/envs/shaplinear
import os
import pickle
import joblib
import shap  # shap v0.48.0
import pandas as pd
import numpy as np  # v1.26.4
import statsmodels.formula.api as smf  # v0.14.2
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import r2_score  # v1.6.1
from glob import glob
from scipy.spatial.distance import euclidean  # v1.15.3
from scipy.cluster.hierarchy import leaves_list

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
# TABLE S4
################################################################################
# linear regression of single-env model performance on fitness-related factors
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
var = pheno.var(axis=0)
var.name = "var"
med = pheno.median(axis=0)
med.name = "median"

# how many environments have small variance?
cv = (pheno.std(axis=0) / pheno.mean(axis=0) *
      100).sort_values()  # coefficient of variation (%)
cv.describe()

# trait heritabilities
h2 = pd.read_csv("Data/Peter_2018/Heritability_h2_H2_sommer.csv")
h2.set_index("Conditions", inplace=True)
X = pd.concat([var, med, h2.h2], ignore_index=False, axis=1)  # fitness factors
X_s = (X-X.mean())/X.std()  # center and scale

for data_type in ["SNPs", "PAVs", "CNVs", "PCs"]:
    if data_type == "PCs":  # model performance results
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
    else:
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
    #
    Y.set_index("Y", inplace=True)
    #
    # Regress model performance on all factors at once
    # df = pd.concat([X_s, Y.r2_test], axis=1)
    # mod = smf.ols(formula=f"r2_test ~ {' * '.join(df.columns[:3])}", data=df)
    # Regress model performance on trait variance alone
    df = pd.concat([X_s, Y.r2_test], axis=1, ignore_index=False)
    mod = smf.ols(formula=f"r2_test ~ var", data=df)
    res = mod.fit()
    yhats = res.predict()
    # double check sm.ols R2
    r2_scores = r2_score(df.r2_test, yhats, multioutput=None)
    print(data_type, res.pvalues)
    #
    # with open(f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_results.txt", "w") as out:
    with open(f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_results_var_alone.txt", "w") as out:
        out.write(res.summary().as_text())
        out.write(f"\nR-sq: {r2_scores}")
    # vars(res) # attributes
    #
    pickle.dump(mod, open(
        # f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_model.pkl", 'wb'))  # save the model
        f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_model_var_alone.pkl", 'wb'))  # save the model
    #
    yhats = pd.Series(yhats)
    yhats.index = df.index
    yhats.name = 'y_pred'
    pd.concat([Y.r2_test, yhats], ignore_index=False, axis=1).\
        to_csv(
            # f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_preds.csv")
            f"Scripts/Data_Vis/Section_2/Table_S3_factors_ols_{data_type}_preds_var_alone.csv")
    #
    del df, Y, mod

# SNPs Intercept    2.402510e-09
# var          2.116161e-02 << p-value of fitness variance
# dtype: float64
# 1837
# 26
# PAVs Intercept    2.074725e-09
# var          4.068378e-02
# dtype: float64
# 1837
# 25
# CNVs Intercept    3.954441e-10
# var          4.873478e-04
# dtype: float64
# 1837
# 26
# PCs Intercept    1.775204e-08
# var          5.051

# Regression model performance on the number of features
feat_nums = {}
for data_type in ["SNPs", "PAVs", "CNVs"]:
    if data_type == "PCs":  # model performance results
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
    else:
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
    #
    Y.set_index("Y", inplace=True)
    feat_nums[data_type] = Y.FeatureNum
    #
    mod = smf.ols(formula=f"r2_test ~ FeatureNum", data=Y)
    res = mod.fit()
    print(data_type, res.pvalues)
    pickle.dump(mod, open(
        f"Scripts/Data_Vis/Section_2/Table_S4_featnum_ols_{data_type}_model.pkl", 'wb'))
    with open(f"Scripts/Data_Vis/Section_2/Table_S4_featnum_ols_{data_type}_results.txt", "w") as out:
        out.write(res.summary().as_text())

# SNPs Intercept     5.247553e-07
# FeatureNum    1.827464e-04
# PAVs Intercept     0.006707
# FeatureNum    0.177538
# CNVs Intercept     0.000002
# FeatureNum    0.947710
# PCs Intercept     3.709167e-08
# FeatureNum    3.709167e-08

feat_nums = pd.DataFrame(feat_nums)
feat_nums.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S4_featnums_matrix.csv")

# Interpret the statsmodels linear models with SHAP values
# first need to build the design matrix and add the interaction terms


def make_design_mat(X):
    X_design = X.copy()
    X_design['var:median'] = X_design['var'] * X_design['median']
    X_design['var:h2'] = X_design['var'] * X_design['h2']
    X_design['median:h2'] = X_design['median'] * X_design['h2']
    X_design['var:median:h2'] = X_design['var'] * \
        X_design['median'] * X_design['h2']
    return X_design


Xs_design = make_design_mat(X_s)  # standardized factor values
X_design = make_design_mat(X)  # original factor values
# X_design.to_csv(
#     "Scripts/Data_Vis/Section_2/Table_S4_factors_design_matrix.csv")

fig, ax = plt.subplots(nrows=7, ncols=4, figsize=(14, 22))
points_to_annotate = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500",
                      "YPDCUSO410MM", "YPDSODIUMMETAARSENITE", "YPDANISO10"]

for j, data_type in enumerate(["PCs", "SNPs", "PAVs", "CNVs"]):
    if data_type == "PCs":  # model performance results
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")
    else:
        Y = pd.read_csv(
            "Scripts/Data_Vis/Section_2/RESULTS_RF_%s_FS.txt" % data_type, sep="\t")
    #
    # Load the fitted model, obtain the parameters, and compute SHAP values
    mod = joblib.load(open(
        f"Scripts/Data_Vis/Section_2/Table_S4_factors_ols_{data_type}_model.pkl", "rb"))
    y_preds = pd.Series(mod.fit().predict(), index=X_s.index, name="y_pred")
    # r2 = r2_score(Y.set_index("Y").loc[y_preds.index, "r2_test"], y_preds)
    # print(f"{data_type} R2: {r2}")
    #
    params = pd.Series(mod.fit().params)
    params = pd.DataFrame({"intercept": params.Intercept,
                          "coef": params.drop("Intercept")})
    explainer = shap.LinearExplainer(
        (params.coef.to_numpy(), params.intercept.to_numpy()), Xs_design)
    shap_values = explainer(Xs_design)
    shap_values.shape  # (35, 7)
    # shap_values = pd.DataFrame(shap_values.values, columns=Xs_design.columns, index=Xs_design.index)
    #
    # Visualise the SHAP values
    # Cluster the SHAP values
    clustering = shap.utils.hclust(
        Xs_design, Y.set_index("Y").loc[Xs_design.index, "r2_test"],
        linkage="average", metric=euclidean)
    feat_order = leaves_list(clustering)
    clustering = shap.utils.hclust(
        Xs_design.T, Y.set_index("Y").loc[Xs_design.index, "r2_test"],
        linkage="average", metric=euclidean)
    instance_order = leaves_list(clustering)
    shap_ordered = shap_values.values[instance_order, :][:, feat_order]
    Xs_design_ordered = Xs_design.iloc[instance_order, feat_order]
    X_design_ordered = X_design.iloc[instance_order, feat_order]
    shap_ordered = pd.DataFrame(shap_ordered, index=Xs_design_ordered.index,
                                columns=Xs_design_ordered.columns)
    shap_ordered.to_csv(
        f"Scripts/Data_Vis/Section_2/Table_S4_factors_shap_{data_type}_matrix.csv")
    #
    # # Correlation of shap values with model performance
    # pd.concat([Y.set_index('Y').loc[shap_ordered.index, 'r2_test'], shap_ordered],
    # 	axis=1, ignore_index=False).corr().to_csv(
    # 	f"Scripts/Data_Vis/Section_2/Table_S3_factors_shap_{data_type}_corr_performance.csv")
    #
    # Correlation of shap values with X_design_ordered and r2_test
    shap_factors = pd.concat(
        [Y.set_index("Y").loc[X_design_ordered.index, "r2_test"],
         X_design_ordered, shap_ordered], axis=1,
        ignore_index=False)
    shap_factors.columns = ["r2_test"] + X_design_ordered.columns.tolist() + \
        [f"shap_{col}" for col in shap_ordered.columns]
    # shap_factors.corr().to_csv(
    # 	f"Scripts/Data_Vis/Section_2/Table_S3_factors_shap_{data_type}_corr_factors.csv")
    #
    # ## Heatmap with SHAP values annotated
    # # ax = shap.plots.heatmap(shap_values, feature_order=feat_order,
    # # 	instance_order=instance_order, show=False)
    # fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8.5, 10), sharey=True)
    # sns.heatmap(shap_ordered, cmap="RdBu_r", center=0, vmin=-.4, vmax=.4,
    # 	xticklabels=shap_ordered.columns, yticklabels=shap_ordered.index,
    # 	cbar_kws={"label": "SHAP value"}, square=True, ax=ax[0]) #, annot=True, fmt=".3f")
    # sns.heatmap(X_design_ordered, cmap="RdBu_r", square=True,
    # 	xticklabels=X_design_ordered.columns.tolist(),
    # 	yticklabels=X_design_ordered.index, cbar_kws={"label": "Feature value"},
    # 	ax=ax[1]) #, annot=True, fmt=".3f")
    # sns.heatmap(Y.set_index("Y").loc[X_design_ordered.index, "r2_test"].to_frame(),
    # 	cmap="RdBu_r", xticklabels=["r2_test"], yticklabels=X_design_ordered.index,
    # 	cbar_kws={"label": "Model performance (test R2)"}, ax=ax[2], square=True)
    # # plt.tight_layout()
    # plt.savefig(f"Scripts/Data_Vis/Section_2/Table_S4_factors_shap_{data_type}_heatmap.pdf")
    # plt.close('all')
    #
    # ## Waterfall plots
    # for i, instance in enumerate(Xs_design.index):
    # 	shap_instance = shap.Explanation(
    # 		values=shap_values.values[i],
    # 		base_values=shap_values.base_values[i][0], # the base value is the same for all features
    # 		data=shap_values.data[i],
    # 		feature_names=shap_values.feature_names)
    # 	shap.plots.waterfall(shap_instance, show=False)
    # 	plt.title(instance, fontsize=7)
    # 	plt.savefig(f"Scripts/Data_Vis/Section_2/Table_S4_factors_shap_{data_type}_{instance}_waterfall.pdf")
    # 	plt.close('all')
    #
    # Scatter plots of SHAP values vs feature values
    # fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(8, 11))
    # points_to_annotate = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500",
    # 	"YPDCUSO410MM", "YPDSODIUMMETAARSENITE", "YPDANISO10"]
    #
    for i, feature in enumerate(X_design_ordered.columns):
        sns.scatterplot(x=shap_factors[feature], y=shap_factors[f"shap_{feature}"],
                        hue=Y.set_index(
                            "Y").loc[shap_factors.index, "r2_test"],
                        ax=ax[i, j], palette="RdBu_r", alpha=0.7, edgecolor="black", legend=True)
        #
        # Draw a red regression line
        sns.regplot(x=shap_factors[feature], y=shap_factors[f"shap_{feature}"],
                    ax=ax[i, j], scatter=False, color="red", line_kws={"linestyle": "--"},
                    seed=i,)
        ax[i, j].set_xlabel(f"{feature} values")
        ax[i, j].set_ylabel(
            f"SHAP values for {feature} from {data_type} models")
        #
        # Annotate a select number of points
        for point in points_to_annotate:
            if point in shap_factors.index:
                ax[i, j].annotate(point, xy=(shap_factors[feature][point],
                                  shap_factors[f"shap_{feature}"][point]), fontsize=6,
                                  xytext=(1, 1), textcoords='offset points', color="black")
    #
    # set aspect ratio to box
    ax[i, j].set_box_aspect(1)
    #
    # plt.tight_layout()
    # plt.savefig(
    #     f"Scripts/Data_Vis/Section_2/Table_S4_factors_shap_{data_type}_scatter.pdf")
    # plt.close("all")
    #
    del Y, mod, explainer, shap_values, clustering, Xs_design_ordered
    del X_design_ordered, shap_factors, shap_ordered

plt.tight_layout()
plt.savefig(
    f"Scripts/Data_Vis/Section_2/Table_S4_factors_shap_scatter.pdf")
plt.close("all")
