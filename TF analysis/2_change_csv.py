import os
import numpy as np
import scanpy as sc
import loompy as lp

regions = ['IRI_12h_1','IRI_12h_2','IRI_24h_1','IRI_24h_2']
suffixes = [ "total"]
for region in regions:
    for suffix in suffixes:
        filename = f"{region}_{suffix}.csv"
        filepath = os.path.join("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/data/raw_counts", filename)
        x = sc.read_csv(filepath)
        row_attrs = {"Gene": np.array(x.var_names)}
        col_attrs = {"CellID": np.array(x.obs_names)}
        loom_filename = f"{region}_{suffix}.loom"
        loom_filepath = os.path.join("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/data/raw_counts", loom_filename)
        lp.create(loom_filepath, x.X.transpose(), row_attrs, col_attrs)
