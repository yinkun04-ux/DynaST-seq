sample_folder=/home/disk/chengrui/workspace/时空转录组/data/

for sample in "XM-ZWY-ST146" "XM-ZWY-ST147" "XM-QYF-ST169" "XM-QYF-ST170" "XM-QYF-ST152" "XM-QYF-ST153" "XM-QYF-ST155" "XM-QYF-ST156";
do
    sample=${sample_folder}/${sample}
    gzip -dc ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv
    awk '{print substr($0, 1, 8) "_" substr($0, 9)}' ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv > ${sample}/outs/filtered_feature_bc_matrix/output.tsv
    /home/disk/chengrui/.conda/envs/d2l-zh/bin/velocyto run \
        -b ${sample}/outs/filtered_feature_bc_matrix/output.tsv \
        -o ../result \
        -m /home/disk/chengrui/toolfolder/时空参考/mouse.gtf ${sample}/Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam /home/disk/chengrui/toolfolder/时空参考/mm10_gencodeM13/gencode.vM13.annotation.gtf
done
wait