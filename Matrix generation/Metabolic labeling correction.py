import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from scipy.stats import binom
from scipy.optimize import minimize, minimize_scalar
import ngs_tools as ngs

CONVERSION_COLUMNS = ['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG']
BASE_COLUMNS = ['A', 'C', 'G', 'T']

# Load data
def load_data(mutations_file, reads_file, group_info_file, convs_file):
    mutations_df = pd.read_csv(mutations_file)
    reads_df = pd.read_csv(reads_file)
    group_meta = pd.read_csv(group_info_file)
    conversions = ngs.utils.read_pickle(convs_file)
    
    return mutations_df, reads_df, group_meta, conversions

# Data Preprocessing
def preprocess_data(mutations_df, reads_df):
    count_df = mutations_df.groupby(['read_id', 'conversion']).size().reset_index(name='count')
    pivot_df = count_df.pivot_table(index='read_id', columns='conversion', values='count', fill_value=0).reset_index()
    
    merged_df = pd.merge(reads_df, pivot_df, on='read_id', how='left').fillna(0)
    merged_df['conversion_sum'] = merged_df[CONVERSION_COLUMNS].sum(axis=1)
    merged_df['has_conversions'] = merged_df[CONVERSION_COLUMNS].sum(axis=1) > 0
    merged_df = merged_df.sort_values(['has_conversions', 'conversion_sum']).drop(columns=['conversion_sum', 'has_conversions'])
    
    return merged_df

# Aggregating conversion data for each combination
def aggregate_conversions(df_sorted, conversions, all_conversions):
    results = []
    for convs in conversions:
        convs = sorted(convs)
        other_convs = list(set(all_conversions) - set(convs))
        bases = list(set(f[0] for f in convs))
        
        df_filtered = df_sorted[(df_sorted[other_convs] == 0).all(axis=1)] if other_convs else df_sorted
        df_combined = df_filtered[['barcode', 'gene_id']].copy()
        df_combined['conversion'] = df_filtered[list(convs)].sum(axis=1).astype(int)
        df_combined['base'] = df_filtered[bases].sum(axis=1)
        df_combined = df_combined.groupby(['barcode', 'gene_id', 'conversion', 'base'], sort=False).size().reset_index(name='count')
        
        results.append(df_combined)
    return results

# Merge group info into the data
def merge_group_info(df_sorted, df_combined_list, group_meta):
    group_meta.columns = ['barcode', 'group']
    df_combined = pd.concat(df_combined_list, ignore_index=True)
    df_sorted = pd.merge(df_sorted, group_meta, on='barcode', how='left')
    df_combined = pd.merge(df_combined, group_meta, on='barcode', how='left')
    
    return df_sorted, df_combined

# Calculate mean conversion rates for each group
def calculate_mean_conversion_rate(df_sorted, group_meta, convs):
    group_pe_values = {}
    for group_name, group_data in df_sorted.groupby('group'):
        group_sum = group_data.sum(numeric_only=True).astype(np.uint32).T
        conversion_columns = [conv for conv in convs if conv not in ['TC', 'TA', 'TG']]
        
        for conversion in conversion_columns:
            group_sum[conversion] /= group_sum[conversion[0]]  # Normalizing
        if conversion_columns:
            overall_pe_value = group_sum[conversion_columns].mean()  # Remove axis=1
            group_pe_values[group_name] = overall_pe_value
    result_df = pd.DataFrame(group_pe_values.items(), columns=['Group', 'Mean Conversion Rate'])
    return result_df


# Sampling data
def sample_data(filtered_df, subset_size):
    sampled_dfs = []
    for group, group_df in filtered_df.groupby('group'):
        if len(group_df) < subset_size:
            print(f"Warning: Only {len(group_df)} samples available for group {group}. Unable to sample {subset_size}.")
            sampled_group_df = group_df
        else:
            sampled_group_df = group_df.sample(n=subset_size, random_state=1)
        sampled_dfs.append(sampled_group_df)
    df_combined_group_data = pd.concat(sampled_dfs)
    return df_combined_group_data

# Log likelihood function for optimization
def LL(params, data):
    t1 = binom.pmf(data[:, 1], data[:, 0], params[0])
    t2 = binom.pmf(data[:, 1], data[:, 0], params[1])
    l = params[2] * t1 + (1 - params[2]) * t2
    l = np.clip(l, 1e-10, 1)
    return -np.sum(np.log(l))
def LL_2(params, data, p, q):
    t1 = binom.pmf(data[:, 1], data[:, 0], p)
    t2 = binom.pmf(data[:, 1], data[:, 0], q)
    l = params * t1 + (1 - params) * t2
    l = np.sum(np.log(l))
    return -l
# Estimate p and q values for each group
def estimate_pq_values(df_combined_group_data, group_pe_values):
    results = []
    for idx, row in group_pe_values.iterrows():
        group = row['Group']  # 获取每行的 group 名称
        q = row['Mean Conversion Rate']  # 获取每行的 conversion rate 值
        pqdata = df_combined_group_data[df_combined_group_data['group'] == group][['base', 'conversion']].values
        if len(pqdata) == 0:
            continue
        
        r_mtx = []
        for _ in range(100):
            init_p = np.random.uniform(0.01, 0.99)
            result = minimize(lambda x: LL([init_p, q, x], pqdata), x0=0.01, bounds=[(0.01, 0.99)])
            par = [init_p, q, result.x]
            value = LL(par, pqdata)
            if 0 <= par[2] <= 1:
                r_mtx.append([par[0], par[2], value])
        
        r_mtx_df = pd.DataFrame(r_mtx)
        min_idx = np.argmin(r_mtx_df.iloc[:, 2])
        p = r_mtx_df.iloc[min_idx, 0]
        results.append({'group': group, 'q_value': q, 'optimal_p': p})
    
    return pd.DataFrame(results)

# Calculate theta for each gene
def calculate_theta_for_gene(results_df, genedata, df_combined):
    all_theta_for_gene = pd.DataFrame()
    for index, row in results_df.iterrows():
        group = row['group']
        q_value = row['q_value']
        optimal_p = row['optimal_p']
        
        genedata_group = genedata[genedata['group'] == group].groupby('gene_id')
        theta_for_gene = []
        gene_ids = []
        for gene_id, group_data in genedata_group:
            if len(group_data) >= 100:
                result = minimize_scalar(lambda x: LL_2(x, group_data[['base', 'conversion']].values, optimal_p, q_value), bounds=(0.01, 0.99), method='bounded')
                theta_for_gene.append(result.x)
                gene_ids.append(gene_id)
        theta_for_gene_df = pd.DataFrame({'gene_id': gene_ids, 'theta': theta_for_gene, 'group': group}).set_index('gene_id')
        if optimal_p < q_value:
            theta_for_gene_df['theta'] = 1 - theta_for_gene_df['theta']
        all_theta_for_gene = pd.concat([all_theta_for_gene, theta_for_gene_df])
    
    return all_theta_for_gene

# Boxplot visualization of Theta values
def plot_boxplot_of_theta(all_theta_for_gene, out_dir):
    plt.figure(figsize=(6, 4))
    sns.boxplot(x='group', y='theta', data=all_theta_for_gene, palette=['#00A087FF', '#4DBBD5FF', '#E64B35FF', '#a65628'])
    plt.title('Boxplot of Theta by Group')
    plt.xlabel(' ')
    plt.ylabel('Theta')
    plt.grid(False)
    plt.savefig(os.path.join(out_dir, "Boxplot of Theta by Group.pdf"), bbox_inches='tight')
    plt.show()

# Main workflow
def main(mutations_file, reads_file, group_info_file, convs_file, out_dir):
    mutations_df, reads_df, group_meta, conversions = load_data(mutations_file, reads_file, group_info_file, convs_file)
    df_sorted = preprocess_data(mutations_df, reads_df)
    all_conversions = sorted(ngs.utils.flatten_iter(conversions))
    df_combined_list = aggregate_conversions(df_sorted, conversions, all_conversions)
    df_sorted, df_combined = merge_group_info(df_sorted, df_combined_list, group_meta)
    
    group_mean_conversion_rate = calculate_mean_conversion_rate(df_sorted, group_meta, CONVERSION_COLUMNS)
    group_mean_conversion_rate.to_csv(os.path.join(out_dir, 'grouped_mean_conversion_rate.csv'), index=False)
    
    filtered_df = df_combined[(df_combined['conversion'] != 0) & (df_combined['base'] != 0)]
    subset_size = 100000
    df_combined_group_data = sample_data(filtered_df, subset_size)
    
    group_pe_values = calculate_mean_conversion_rate(df_sorted, group_meta, CONVERSION_COLUMNS)
    results_df = estimate_pq_values(df_combined_group_data, group_pe_values)
    
    genedata = df_combined[['base', 'conversion', 'gene_id', 'group']]
    group_pe_values.to_csv(os.path.join(out_dir, "group_pe_values.csv"))
    
    all_theta_for_gene = calculate_theta_for_gene(results_df, genedata, df_combined)
    all_theta_for_gene.to_csv(os.path.join(out_dir, "all_theta_for_gene.csv"))
    
    plot_boxplot_of_theta(all_theta_for_gene, out_dir)
    
    print("############# Finished ##################")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process conversion rate data")
    parser.add_argument('--mutations_file', help="Path to mutations CSV file")
    parser.add_argument('--reads_file', help="Path to reads CSV file")
    parser.add_argument('--group_info_file', help="Path to group information CSV file")
    parser.add_argument('--convs_file', help="Path to conversions pickle file")
    parser.add_argument('--out_dir', help="Path to output directory")
    
    args = parser.parse_args()
    
    main(args.mutations_file, args.reads_file, args.group_info_file, args.convs_file, args.out_dir)