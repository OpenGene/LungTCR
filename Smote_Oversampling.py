#!/users/miniconda3/envs/plot/bin/python
# coding: utf-8


import pandas as pd
import sys
import os
import random
from imblearn.over_sampling import SMOTENC,BorderlineSMOTE
from sklearn.preprocessing import LabelEncoder


def get_index(input_tbl:pd.DataFrame, transformed_df:pd.DataFrame, index_col:str):
    le = LabelEncoder()
    id_lst = le.fit_transform(input_tbl['Sample_ID'])
    transformed_df['index'] = le.inverse_transform(transformed_df['index'])
    transformed_df.set_index('index')
    return transformed_df
    

input_file = os.path.abspath(sys.argv[1])
features_lst_file = os.path.abspath(sys.argv[2])
output_file = os.path.abspath(sys.argv[3])
method = sys.argv[4]  #smotenc
RAND_SEED = 77
BINARY_FEATURES = ["somking","quit_somking","lung_cancer_history","family_history","COPD","left","right","up","down","middle","SN","pSN","GGN", "glitch"]

input_tbl = pd.read_csv(input_file)
feature_lst = pd.read_csv(features_lst_file)['feature'].unique().tolist() # the order of the input feature list will influence the output data
category_cols = list(set(feature_lst) & set(BINARY_FEATURES))
type_counts = input_tbl['Type'].value_counts().to_dict()
# ratio = 'minority'
NEIGHBORS = 5
ratio = type_counts['Benign'] / type_counts['Malignant']
if ratio > 1:
    ratio = 1 / ratio
ratio = min(ratio*3, max(0.7, min(ratio*2, 1)))
if not category_cols and method == 'smotenc':
    print('No category cols change method to borderlinesmote')
    method = 'borderlinesmote'
if method == 'borderlinesmote':
    # sm = BorderlineSMOTE(random_state=RAND_SEED, k_neighbors=5, sampling_strategy='minority')
    sm = BorderlineSMOTE(random_state=RAND_SEED, k_neighbors=NEIGHBORS, sampling_strategy=ratio)
elif method == 'smotenc':
    # sm = SMOTENC(random_state=RAND_SEED, categorical_features=category_cols, k_neighbors=5, sampling_strategy='minority')
    sm = SMOTENC(random_state=RAND_SEED, categorical_features=category_cols, k_neighbors=NEIGHBORS, sampling_strategy=ratio)

input_tbl.set_index('Sample_ID', inplace=True)
print(f'Before {method}:')
print(input_tbl['Type'].value_counts().to_dict())
X, y = sm.fit_resample(input_tbl[feature_lst], input_tbl['Type'])
output_df = pd.DataFrame(X)
print(f'After {method}:')
print(y.value_counts().to_dict())
output_df['Type'] = y
output_df['Sample_ID'] = X.index
output_df.to_csv(output_file, index=None)
