import pandas as pd
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.model_selection import train_test_split

ebv_file = open(sys.argv[1], 'r')
df_ebv_file = pd.read_csv(ebv_file)
df_ebv_file['Study ID'] = df_ebv_file['Study ID'].astype(str)
#df_ebv_file_int = df_ebv_file.astype('Int64')

human_file = open(sys.argv[2], 'r')
df_human_file = pd.read_csv(human_file)
df_human_file['Study ID'] = df_human_file['Study ID'].astype(str)
#df_human_file_int = df_human_file.astype('Int64')

combined_df = pd.merge(df_ebv_file, df_human_file, on='Study ID')

#print(combined_df)
#only filtering for geography=1 aka African samples + removing Brazilian patients
combined_df = combined_df.drop(combined_df[combined_df['Geography'] == 0].index)
combined_df = combined_df.drop(combined_df[combined_df['Location'] == 'Belo Horizonte'].index)


combined_df = combined_df.drop(columns=combined_df.columns[combined_df.eq(0).mean()>0.95])
#print(combined_df)

#so that we only get columns with 2 or fewer missing values
combined_df.dropna(axis=1, thresh = int(0.98*combined_df.shape[0]), inplace=True)
#combined_df.to_csv('df_ebv_and_human_mutations_231115_240225trial_Africanonly.csv', index=False)

#imputation with median value
combined_df[combined_df.columns] = combined_df[combined_df.columns].apply(pd.to_numeric, errors='ignore')

combined_df.fillna(combined_df.median(numeric_only=True), inplace=True)

na_column_df = combined_df.isna().sum()
na_column_df.to_csv('df_na_column.csv', index=False)

combined_df.drop(['Unnamed: 0_x', 'Unnamed: 0_y'], axis=1, inplace=True)

#combined_df.to_csv('df_ebv_and_human_mutations_geography1_231115_Africaonly_240226.csv', index=False)

output_variable_file = open("/Users/isaacekimjr/Desktop/231115_varnames.txt", 'w+')
for varname in list(combined_df.columns):
    if "Unnamed" not in varname:
        output_variable_file.write(varname.lower()[:34])
        output_variable_file.write('\n')

#print(combined_df)
