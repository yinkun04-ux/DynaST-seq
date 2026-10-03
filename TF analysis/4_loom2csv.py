import loompy as lp
import pandas as pd
import sys
import os
sample_path="/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/program"
out_path="/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/result/scenic_st"
input = os.path.join(sample_path, sys.argv[1]+ '.sample_SCENIC.loom')
lf = lp.connect(input, mode='r+', validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID) 
lf.close()
output = os.path.join(out_path, sys.argv[1]+ '_SCENIC.csv')
auc_mtx.to_csv(output, sep=',', header=True, index=True)