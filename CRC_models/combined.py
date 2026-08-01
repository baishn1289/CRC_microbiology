#!/usr/bin/env python
# coding: utf-8

# # import packages


import numpy as np
import pandas as pd
import os

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

from sklearn.ensemble import RandomForestClassifier
from catboost import CatBoostClassifier
from tabpfn import TabPFNClassifier

import matplotlib.pyplot as plt
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import classification_report, roc_auc_score, average_precision_score

from sklearn.preprocessing import OneHotEncoder

from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_curve, roc_auc_score

from sklearn.model_selection import LeaveOneGroupOut


# ### Load data


patient_data = pd.read_csv(
    '../data/0.clinical_data.csv',
)

patient_data


for s in range(len(patient_data)):
    if patient_data.iloc[s,2] == 'male':
        patient_data.iloc[s,2] = 1
    elif patient_data.iloc[s,2] == 'female':        
        patient_data.iloc[s,2] = 0
    else:
        raise ValueError(f"Unexpected value in gender column.")

patient_data


def data_summary(df):

    summary = df.describe(include='all').T

    summary['Missing Values'] = df.isnull().sum()
    summary['Missing Rate (%)'] = (df.isnull().mean() * 100).round(2)

    summary['Data Type'] = df.dtypes

    cols = ['Data Type', 'Missing Values', 'Missing Rate (%)'] + [
        c for c in summary.columns
        if c not in ['Data Type', 'Missing Values', 'Missing Rate (%)']
    ]

    return summary[cols]


summary_df = data_summary(patient_data)
print(summary_df)


missing_bmi_samples = patient_data[patient_data['BMI'].isnull()]

missing_bmi_samples


judgement = pd.read_csv(
    '../data/0.model_cohort_determine.csv'
)

judgement


for i in range(len(judgement)):
    if judgement.iloc[i,5] < 15:
        judgement.iloc[i,7] = 'yes'
        judgement.iloc[i,9] = 'no'

judgement


df0 = pd.read_csv(
    '../data/input_data0.ml_input_data_raw.csv',
)

df0


df0 = df0.T
df0.columns = df0.iloc[0]
df0 = df0[1:]
df0 = df0.reset_index(names="Index")
df0.rename(columns={df0.columns[0]: "Sample_ID"}, inplace=True)
df0


aaa = df0.drop(columns=["Sample_ID"])
global_min = aaa.replace(0, np.inf).to_numpy().min()

print(f"min value (excluding zeros): {global_min}")


df = df0.drop(columns=["Sample_ID"])

df = (df.apply(pd.to_numeric, errors="coerce")
                     .replace([np.inf, -np.inf], np.nan)
                     .fillna(0.0)
                     .astype(np.float64))

row_sums = df.sum(axis=1)
X_closed = df.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

eps = 1e-6
X_pc = X_closed + eps

logX = np.log(X_pc)
clr = logX.sub(logX.mean(axis=1), axis=0)

df = clr

df = pd.DataFrame(df, index=df0.index, columns=df0.columns)
df["Sample_ID"] = df0["Sample_ID"]

df


first_dir_name = "../data/"
os.listdir(first_dir_name)


# ### Meta_Sig_EO


EO = pd.DataFrame()


second_dir_name = "LODO_Meta_Sig_EO/"
full_dir_name = first_dir_name + second_dir_name
MEO_files = os.listdir(full_dir_name)
pre = full_dir_name
cohort = []

for post in MEO_files:
    test_data_name = post[13:-4]
    print(test_data_name)
    cohort.append(test_data_name)

EO["cohort"] = cohort


All_record_tab = {}

for post in MEO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "EO_test"].iloc[0] == 'no':
        All_record_tab[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MEO_feature = pd.read_csv(address)
    
    selected_columns = ["Sample_ID"] + MEO_feature["Feature"].tolist()
    df_MEO_feature = df[selected_columns].copy()
    patient_data_MEO = patient_data[patient_data["Age_class"].eq("EO")]
    target_MEO = patient_data_MEO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MEO.rename(columns={"Group": "labels"}, inplace=True)
    df_MEO_feature["Sample_ID"] = df_MEO_feature["Sample_ID"].astype(str)
    target_MEO["Sample_ID"] = target_MEO["Sample_ID"].astype(str)
    merged_MEO_feature = pd.merge(df_MEO_feature, target_MEO, left_on="Sample_ID", right_on="Sample_ID", how="right")

    X_full = merged_MEO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MEO_feature["labels"]
    Groups = merged_MEO_feature["Cohort"]
    
    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    tc = TabPFNClassifier(model_path="...", # The location of the model weights on your computer
                          device="cuda") 
    tc.fit(X_train, y_train)
    proba = tc.predict_proba(X_test)
    auc = roc_auc_score((y_test == tc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_tab[test_data_name] = auc


for key in All_record_tab:
    print(f"{All_record_tab[key]}")
EO["Tab-v2.5-3"] = All_record_tab.values()


All_record_rf = {}

for post in MEO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "EO_test"].iloc[0] == 'no':
        All_record_rf[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MEO_feature = pd.read_csv(address)
    
    selected_columns = ["Sample_ID"] + MEO_feature["Feature"].tolist()
    df_MEO_feature = df[selected_columns].copy()
    patient_data_MEO = patient_data[patient_data["Age_class"].eq("EO")]
    target_MEO = patient_data_MEO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MEO.rename(columns={"Group": "labels"}, inplace=True)
    df_MEO_feature["Sample_ID"] = df_MEO_feature["Sample_ID"].astype(str)
    target_MEO["Sample_ID"] = target_MEO["Sample_ID"].astype(str)
    merged_MEO_feature = pd.merge(df_MEO_feature, target_MEO, left_on="Sample_ID", right_on="Sample_ID", how="right")
    
    X_full = merged_MEO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MEO_feature["labels"]
    Groups = merged_MEO_feature["Cohort"]
    
    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    rfc = RandomForestClassifier(n_estimators=500, random_state=42)
    rfc.fit(X_train, y_train)
    proba = rfc.predict_proba(X_test)
    auc = roc_auc_score((y_test == rfc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_rf[test_data_name] = auc


for key in All_record_rf:
    print(f"{All_record_rf[key]}")
EO["RF"] = All_record_rf.values()


All_record_cb = {}

for post in MEO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "EO_test"].iloc[0] == 'no':
        All_record_cb[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MEO_feature = pd.read_csv(address)
    
    selected_columns = ["Sample_ID"] + MEO_feature["Feature"].tolist()
    df_MEO_feature = df[selected_columns].copy()
    patient_data_MEO = patient_data[patient_data["Age_class"].eq("EO")]
    target_MEO = patient_data_MEO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MEO.rename(columns={"Group": "labels"}, inplace=True)
    df_MEO_feature["Sample_ID"] = df_MEO_feature["Sample_ID"].astype(str)
    target_MEO["Sample_ID"] = target_MEO["Sample_ID"].astype(str)
    merged_MEO_feature = pd.merge(df_MEO_feature, target_MEO, left_on="Sample_ID", right_on="Sample_ID", how="right")
    
    X_full = merged_MEO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MEO_feature["labels"]
    Groups = merged_MEO_feature["Cohort"]

    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    cbc = CatBoostClassifier(iterations=500, learning_rate=0.05, depth=6, verbose=False)
    cbc.fit(X_train, y_train)
    proba = cbc.predict_proba(X_test)
    auc = roc_auc_score((y_test == cbc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_cb[test_data_name] = auc


for key in All_record_cb:
    print(f"{All_record_cb[key]}")
EO["CB"] = All_record_cb.values()


EO.to_excel('Combined_EO.xlsx', index=False)


EO


# ### Meta_Sig_LO


LO = pd.DataFrame()


second_dir_name = "LODO_Meta_Sig_LO/"
full_dir_name = first_dir_name + second_dir_name
MLO_files = os.listdir(full_dir_name)
pre = full_dir_name
cohort = []

for post in MLO_files:
    test_data_name = post[13:-4]
    print(test_data_name)
    cohort.append(test_data_name)

LO["cohort"] = cohort


All_record_tab = {}

for post in MLO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "LO_test"].iloc[0] == 'no':
        All_record_tab[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MLO_feature = pd.read_csv(address)
    
    selected_columns = ["Sample_ID"] + MLO_feature["Feature"].tolist()
    df_MLO_feature = df[selected_columns].copy()
    patient_data_MLO = patient_data[patient_data["Age_class"].eq("LO")]
    target_MLO = patient_data_MLO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MLO.rename(columns={"Group": "labels"}, inplace=True)
    df_MLO_feature["Sample_ID"] = df_MLO_feature["Sample_ID"].astype(str)
    target_MLO["Sample_ID"] = target_MLO["Sample_ID"].astype(str)
    merged_MLO_feature = pd.merge(df_MLO_feature, target_MLO, left_on="Sample_ID", right_on="Sample_ID", how="right")
    
    X_full = merged_MLO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MLO_feature["labels"]
    Groups = merged_MLO_feature["Cohort"]
    
    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    tc = TabPFNClassifier(model_path="...", # The location of the model weights on your computer
                          device="cuda") 
    tc.fit(X_train, y_train)
    proba = tc.predict_proba(X_test)
    auc = roc_auc_score((y_test == tc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_tab[test_data_name] = auc


for key in All_record_tab:
    print(f"{All_record_tab[key]}")
LO["Tab-v2.5-3"] = All_record_tab.values()


All_record_rf = {}

for post in MLO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "LO_test"].iloc[0] == 'no':
        All_record_rf[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MLO_feature = pd.read_csv(address)
     
    selected_columns = ["Sample_ID"] + MLO_feature["Feature"].tolist()
    df_MLO_feature = df[selected_columns].copy()
    patient_data_MLO = patient_data[patient_data["Age_class"].eq("LO")]
    target_MLO = patient_data_MLO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MLO.rename(columns={"Group": "labels"}, inplace=True)
    df_MLO_feature["Sample_ID"] = df_MLO_feature["Sample_ID"].astype(str)
    target_MLO["Sample_ID"] = target_MLO["Sample_ID"].astype(str)
    merged_MLO_feature = pd.merge(df_MLO_feature, target_MLO, left_on="Sample_ID", right_on="Sample_ID", how="right")
    
    X_full = merged_MLO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MLO_feature["labels"]
    Groups = merged_MLO_feature["Cohort"]
    
    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    rfc = RandomForestClassifier(n_estimators=500, random_state=42)
    rfc.fit(X_train, y_train)
    proba = rfc.predict_proba(X_test)
    auc = roc_auc_score((y_test == rfc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_rf[test_data_name] = auc


for key in All_record_rf:
    print(f"{All_record_rf[key]}")
LO["RF"] = All_record_rf.values()


All_record_cb = {}

for post in MLO_files:
    test_data_name = post[13:-4]
    if judgement.loc[judgement["Cohort"] == test_data_name, "LO_test"].iloc[0] == 'no':
        All_record_cb[test_data_name] = 'NaN'
        continue
    address = pre + post
    print(address)
    MLO_feature = pd.read_csv(address)
     
    selected_columns = ["Sample_ID"] + MLO_feature["Feature"].tolist()
    df_MLO_feature = df[selected_columns].copy()
    patient_data_MLO = patient_data[patient_data["Age_class"].eq("LO")]
    target_MLO = patient_data_MLO[["Sample_ID", "Group", "Cohort", "Sex", "BMI"]].copy()
    target_MLO.rename(columns={"Group": "labels"}, inplace=True)
    df_MLO_feature["Sample_ID"] = df_MLO_feature["Sample_ID"].astype(str)
    target_MLO["Sample_ID"] = target_MLO["Sample_ID"].astype(str)
    merged_MLO_feature = pd.merge(df_MLO_feature, target_MLO, left_on="Sample_ID", right_on="Sample_ID", how="right")

    X_full = merged_MLO_feature.drop(columns=["Sample_ID", "labels", "Cohort"])
    
    y_full = merged_MLO_feature["labels"]
    Groups = merged_MLO_feature["Cohort"]
    
    J = Groups == test_data_name

    X_train = X_full[~J].copy()
    y_train = y_full[~J]
    
    X_test = X_full[J].copy()
    y_test = y_full[J]

    train_bmi_median = X_train["BMI"].median()

    X_train["BMI"] = X_train["BMI"].fillna(train_bmi_median)
    X_test["BMI"] = X_test["BMI"].fillna(train_bmi_median)

    cbc = CatBoostClassifier(iterations=500, learning_rate=0.05, depth=6, verbose=False)
    cbc.fit(X_train, y_train)
    proba = cbc.predict_proba(X_test)
    auc = roc_auc_score((y_test == cbc.classes_[1]).astype(int), proba[:, 1])
    print("ROC-AUC:", auc)
    All_record_cb[test_data_name] = auc


for key in All_record_cb:
    print(f"{All_record_cb[key]}")
LO["CB"] = All_record_cb.values()


LO.to_excel('Combined_LO.xlsx', index=False)


LO

