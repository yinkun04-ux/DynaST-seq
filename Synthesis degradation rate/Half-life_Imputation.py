import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def cp10k_normalization(expression_matrix):
    total_counts = expression_matrix.sum(axis=1)
    return expression_matrix.div(total_counts, axis=0) * 10000

def compute_pca_knn(expression_data, n_neighbors=10, n_components=30):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(expression_data.T)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(pca_result)
    distances, indices = nbrs.kneighbors(pca_result)
    return indices

def compute_mean_expression(indices, expression_data):
    mean_expression = np.zeros_like(expression_data)
    for i, neighbors in enumerate(indices):
        mean_expression[:, i] = expression_data.iloc[:, neighbors].mean(axis=1)
    return pd.DataFrame(mean_expression, index=expression_data.index, columns=expression_data.columns)

def compute_half_life(new_rna_mean, total_rna_mean, t=1.5):
    R = new_rna_mean / total_rna_mean
    return -t * np.log(2) / np.log(1 - R)

def impute_half_life_data(group_data):
    detection_rate = group_data.notna().mean(axis=1)
    genes_to_impute = detection_rate >= 0.1
    group_data_to_impute = group_data.loc[genes_to_impute]
    group_data_to_impute.replace([np.inf, -np.inf], np.nan, inplace=True)
    valid_cols = group_data_to_impute.columns[~group_data_to_impute.isna().all()]
    return group_data_to_impute[valid_cols]

def compute_imputed_half_life_mean(indices, group_data_to_impute):
    half_life_mean = np.zeros_like(group_data_to_impute)
    for i, neighbors in enumerate(indices):
        half_life_mean[:, i] = group_data_to_impute.iloc[:, neighbors].mean(axis=1)
    return pd.DataFrame(half_life_mean, index=group_data_to_impute.index, columns=group_data_to_impute.columns)

total_expression_raw = pd.read_csv("../../totalRNA_data.csv", index_col=0)
new_expression_raw = pd.read_csv("../../newRNA_data.csv", index_col=0)
clusters = pd.read_csv("../../meta_data.csv", index_col=0)

# Step 1: Normalize total and new expression matrices using CP10k
total_expression = cp10k_normalization(total_expression_raw)
new_expression = cp10k_normalization(new_expression_raw)

cluster_ids = clusters['region'].values
unique_clusters = np.unique(cluster_ids)

for cluster_id in unique_clusters:
    print(f"Processing cluster: {cluster_id}")
    
    cluster_indices = clusters[clusters['region'] == cluster_id].index
    cluster_total_expression = total_expression.loc[:, cluster_indices]
    cluster_new_expression = new_expression.loc[:, cluster_indices]
    
    X = pd.concat([cluster_total_expression, cluster_new_expression], axis=0).fillna(0)
    print(f"Expression matrix shape: {X.shape}")
    X = X.T

    # Step 2: PCA and Nearest Neighbors
    indices = compute_pca_knn(X)

    # Step 3: Compute mean expression across neighbors for new and total RNA
    new_rna_mean = compute_mean_expression(indices, cluster_new_expression)
    total_rna_mean = compute_mean_expression(indices, cluster_total_expression)

    new_rna_mean.to_csv(f"/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/5_半衰期/result/1010res/{cluster_id}_new_rna_mean.csv")
    total_rna_mean.to_csv(f"/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/5_半衰期/result/1010res/{cluster_id}_total_rna_mean.csv")
    
    # Step 4: Compute half-life for each gene
    half_life = compute_half_life(new_rna_mean, total_rna_mean)
    print("Half-life computation complete.")
    
    # Step 5: Impute missing half-life data
    group_data = half_life
    group_data_to_impute = impute_half_life_data(group_data)

    # Step 6: PCA and Nearest Neighbors for imputed data
    indices = compute_pca_knn(X, n_neighbors=5)
    print("Nearest Neighbors completed.")
    
    # Step 7: Compute imputed half-life mean
    half_life_mean = compute_imputed_half_life_mean(indices, group_data_to_impute)
    
    half_life_mean.to_csv(f"/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/5_半衰期/result/1010res/{cluster_id}_half_life_mean.csv")
