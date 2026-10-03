dir=/home/disk/chengrui/toolfolder/tf_data/
input_dir=/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/data/raw_counts
out_dir=/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/6_转录因子gy/result/scenic

tfs=$dir/mm_mgi_tfs.txt
feather=$dir/mm9-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.mgi-m0.001-o0.0.tbl 

out_name=$1
input_loom=$input_dir/$out_name.loom
ls $tfs  $feather  $tbl  

pyscenic grn \
    --num_workers 2 \
    --output $out_name.adj.sample.tsv \
    --method grnboost2 \
    $input_loom \
    $tfs

pyscenic ctx \
    $out_name.adj.sample.tsv \
    $feather \
    --annotations_fname $tbl \
    --expression_mtx_fname $input_loom \
    --mode "dask_multiprocessing" \
    --output $out_name.reg.csv \
    --num_workers 64 \
    --mask_dropouts 

pyscenic aucell \
    $input_loom \
    $out_name.reg.csv \
    --output $out_name.sample_SCENIC.loom \
    --num_workers 4
