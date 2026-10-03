import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import os
import matplotlib.pyplot as plt
import dynamo as dyn
from scipy.sparse import csr_matrix

def process_loom_files(loom_files):
    ldata_list = []
    for i, loom_file in enumerate(loom_files, start=1):
        ldata = scv.read(loom_file, cache=True)
        barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
        barcodes = [bc[:-1] + f'_{i}' for bc in barcodes]
        barcodes = [bc[:8] + bc[9:] if len(bc) > 8 and bc[8] == '_' else bc for bc in barcodes]
        ldata.obs.index = barcodes
        ldata.var_names_make_unique()
        ldata_list.append(ldata)
    
    ldata = ldata_list[0].concatenate(ldata_list[1:])
    ldata.obs.index = [bc[:-2] for bc in ldata.obs.index]
    sc.write('../result/ldata.h5ad', ldata)
    return ldata

# Define dynamo workflow function
def dynamo_workflow(adata, **kwargs):
    preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
    preprocessor.config_monocle_recipe(adata)
    preprocessor.filter_cells_by_outliers_kwargs["keep_filtered"] = True
    preprocessor.preprocess_adata_monocle(adata, **kwargs)
    dyn.tl.dynamics(adata)
    dyn.tl.reduceDimension(adata)
    dyn.tl.cell_velocities(adata, calc_rnd_vel=True, transition_genes=adata.var_names)
    dyn.vf.VectorField(adata, basis='umap')

def load_data():
    totalRNA_data = pd.read_csv("../totalRNA_data.csv", index_col=0)
    newRNA_data = pd.read_csv("../newRNA_data.csv", index_col=0)
    meta_data = pd.read_csv("../meta_data.csv", index_col=0)
    regions = dict(tuple(meta_data.groupby('region')))
    return totalRNA_data, newRNA_data, meta_data, regions

def create_ann_data(regions, totalRNA_data, newRNA_data):
    ann_data_new_dict_05 = {}
    regions_list = ['Cortex', 'Inner medulla', 'Inner stripe', 'Outer stripe']

    for region in regions_list:
        region_index = regions[region].index
        totalRNA_region = totalRNA_data.loc[:, region_index]
        newRNA_region = newRNA_data.loc[:, region_index]
        
        # Filter cells based on specific suffixes
        col_filter = ~totalRNA_region.columns.str.endswith('_2') & \
                    ~totalRNA_region.columns.str.endswith('_4') & \
                    ~totalRNA_region.columns.str.endswith('_6') & \
                    ~totalRNA_region.columns.str.endswith('_8')
        
        filtered_totalRNA_region = totalRNA_region.loc[:, col_filter]
        filtered_newRNA_region = newRNA_region.loc[:, col_filter]
        
        # Create AnnData object using both totalRNA and newRNA
        ann_data_labeling_region = ad.AnnData(
            X=csr_matrix(filtered_totalRNA_region.T.values),
            var=pd.DataFrame(filtered_totalRNA_region.index, columns=['gene_short_name']),
            obs=pd.DataFrame(filtered_totalRNA_region.columns, columns=['cellname']),
            layers={
                'total': csr_matrix(filtered_totalRNA_region.T.values),
                'new': csr_matrix(filtered_newRNA_region.T.values)
            }
        )
        ann_data_new_dict_05[region] = ann_data_labeling_region

    return ann_data_new_dict_05

def assign_cell_times(adata, time_key='time'):
    last_digit = adata.obs.index.str[-1]
    adata.obs[time_key] = '0h'
    adata.obs.loc[last_digit.isin(['1', '2']), time_key] = '0h'
    adata.obs.loc[last_digit.isin(['3', '4']), time_key] = '4h'
    adata.obs.loc[last_digit.isin(['5', '6']), time_key] = '12h'
    adata.obs.loc[last_digit.isin(['7', '8']), time_key] = '24h'
    adata.obs['nUMI'] = adata.X.sum(axis=1)
    adata.obs['nGene'] = adata.X.getnnz(axis=1)
    adata.obs['label_time'] = 1.5  # Labeling time

def generate_plots(data_dict, suffix=""):
    fig, axes = plt.subplots(ncols=4, nrows=1, constrained_layout=True, figsize=(22, 4))
    for i, (region, adata) in enumerate(data_dict.items()):
        ax = axes[i]
        dyn.pl.streamline_plot(adata, color='time', color_key_cmap='viridis', basis='umap', ax=ax, show_legend='right', save_show_or_return='return')
        ax[0].set_title(region)
    plt.show()

def main():
    loom_files = [
        '../result/Aligned_AJOKR.loom', '../result/Aligned_REIBP.loom', '../result/Aligned_C524I.loom', 
        '../result/Aligned_HPT6L.loom', '../result/Aligned_XHKIV.loom', '../result/Aligned_3N4PH.loom', 
        '../result/Aligned_7XEAV.loom', '../result/Aligned_N566X.loom'
    ]
    
    ldata = process_loom_files(loom_files)
    totalRNA_data, newRNA_data, meta_data, regions = load_data()
    ann_data_dict = create_ann_data(regions, totalRNA_data, newRNA_data)
    
    ##### Spliced data (using exon data) #####
    sample_dict_spliced = {region: ldata[ldata.obs.index.isin(ann_data_dict[region].obs.index)] for region in ann_data_dict}
    data_dict = {}
    
    for sample, adata in sample_dict_spliced.items():
        assign_cell_times(adata)
        test_adata_tmp_DE = adata[:, top_genes['top80']]
        dynamo_workflow(test_adata_tmp_DE)
        test_adata_tmp_DE.obs['time'] = test_adata_tmp_DE.obs['time'].fillna(method='ffill')
        data_dict[sample] = test_adata_tmp_DE
    
    generate_plots(data_dict)

    ##### Labeled data (using new RNA data) #####
    sample_dict_labeled = {region: ann_data_dict[region] for region in ann_data_dict}
    data_dict_NEW = {}
    
    for sample, adata in sample_dict_labeled.items():
        assign_cell_times(adata)
        test_adata_tmp_DE = adata[:, top_genes['top80']]
        dynamo_workflow(test_adata_tmp_DE, tkey='label_time', experiment_type='one-shot')
        test_adata_tmp_DE.obs['time'] = test_adata_tmp_DE.obs['time'].fillna(method='ffill')
        data_dict_NEW[sample] = test_adata_tmp_DE

    generate_plots(data_dict_NEW, suffix="_labeled")

if __name__ == "__main__":
    main()
