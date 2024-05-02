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


# #only keep geography, top EBV variants, and top human genes to show correlation matrix
# correlation_df = combined_df[['Geography', '_(168228,C,G)', '_(168409,T,G)', '_(91461,C,T)', '_(168636,G,T)', '_(134730,T,C)', '_(77563,G,A)', '_(167871,C,T)', '_(136560,C,A)', 'HIST1H2BK', 'P2RY8', 'GNAI2', 'TFAP4', 'IGLV3-25', 'ETS1', 'ID3']]
# print(correlation_df)
# corr = correlation_df.corr()
# plt.figure(figsize=(30,30))
# sns.set(font_scale=3)
# sns.heatmap(corr,annot=True,annot_kws={"fontsize":16})
# plt.savefig("seaborn_correlation_matrix.png")
#


#K-Means Clustering
#initialize kmeans parameters
# kmeans_kwargs = {
# "init": "random",
# "n_init": 10,
# "random_state": 1,
# }
#
# #create list to hold SSE values for each k
# sse = []
# num_clusters = 20
#
# #For calculating sum of squared errors to use elbow method and determine ideal cluster number
# # for k in range(1, num_clusters):
# #     kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
# #     print(kmeans.fit(correlation_df))
# #
# #     sse.append(kmeans.inertia_)
# #
# # #visualize results
# # plt.plot(range(1, num_clusters), sse)
# # plt.xticks(range(1, num_clusters))
# # plt.xlabel("Number of Clusters")
# # plt.ylabel("SSE")
# # plt.savefig("kmeans.png")
#
# #For calculating silhouette scores to determine ideal cluster number
# # sil = []
# # # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
# # for k in range(2, num_clusters+1):
# #   kmeans = KMeans(n_clusters = k).fit(correlation_df)
# #   labels = kmeans.labels_
# #   sil.append(silhouette_score(correlation_df, labels, metric = 'euclidean'))
# #
# # plt.plot(range(1, num_clusters), sil)
# # plt.xticks(range(1, num_clusters))
# # plt.xlabel("Number of Clusters")
# # plt.ylabel("Silhouette Score")
# # plt.savefig("silhouette_scores.png")
# #kmeans.predict(X)
#
# kmeans = KMeans(n_clusters=5, random_state=0).fit(combined_df[['Geography', '_(168476,C,A)', '_(168228,C,G)', '_(108730,C,A)', '_(118745,T,G)', '_(168409,T,G)', '_(72465,G,T)', '_(155389,C,T)', '_(168636,G,T)', '_(134730,T,C)', '_(97411,A,T)', '_(77563,G,A)', '_(167871,C,T)', '_(136560,C,A)', 'IGLV3-25', 'EIF4A1', 'HIST1H2BK', 'SIN3A', 'P2RY8', 'MAP3K9']])
# # Save the labels
# combined_df.loc[:,'cluster_number'] = kmeans.labels_
# print(combined_df)
# combined_df.to_csv('df_ebv_and_human_mutations_cluster_number.csv', index=False)
#survival_df = combined_df[['Study ID', 'Survival', 'Dead Status']]
#print(survival_df)
