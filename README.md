# DynaST-seq
Analysis code for DynaST-seq: transcriptome-wide profiling of spatio-temporal RNA dynamics in tissues.
DynaST-seq
DynaST-seq: a Spatially and Temporally Resolved Transcriptomic Sequencing Approach for Profiling Gene Expression Dynamics in Tissues

Abstract
Understanding the spatio-temporal gene expression is crucial for elucidating the molecular mechanisms that regulate cellular functions. Despite recent advances in spatial transcriptomics, most approaches lack temporal resolution, limiting insights into dynamic gene expression. Here, we present DynaST-seq, an advanced spatio-temporal transcriptomics technology that integrates metabolic RNA labeling with spatial transcriptomics, enabling high-resolution and transcriptome-wide profiling of RNA dynamics across space and time.

Set Up Environment
System: centOS
python 3.8.18
R 4.2.1
samtools 1.6
ngs-tools 1.8.5
pysam 0.20.0
Obtaining the Expression Matrix for DynaST-seq Data
Use the DynamicST tool (https://github.com/DynamicBiosystems/DynamicST) to process spatial transcriptomics data, converting raw image data and sequencing data into valuable gene expression information with a spatial context.
 
 DynamicST count \
    --sample demo \
    --id demo \
    --inputdir demo/input_data \
    --gtf xxx.gtf \
    --transcriptome xxx \
    --image demo/input_data/demo_registedimage.tif \
    --alignment demo/input_data/demo_alignment.json \
    --outputdir demo/output_data

Consensus_sequence_construction.py： Used to group reads with identical UMIs, barcodes, and gene labels in the BAM file, and for each genomic position, select the most frequent base within each UMI group as the consensus base to construct the consensus sequence.

 python Consensus_sequence_construction.py \
    --bam_path demo/output_data/DynamicST_result/Aligned.sortedByCoord.UniqueGene.bam \
    --gtf_file xxx.gtf \
    --ref_fa xxx.fa \
    --temp_out_dir demo/output_data 

Expression_matrix_extraction.py： Extracts gene, barcode, and UMI information from the BAM file; calculates the UMI expression matrices and new-to-old ratios.

 python Expression_matrix_extraction.py \
    --sample_path demo/output_data/DynamicST_result/demo \
    --bam_path demo/output_data/Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam \
    --result_path demo/output_data

Conversion_Proportions.py： Calculates and visualizations the substitution rates, and generates summary reports.

 python Conversion_Proportions.py \
    --sample_path demo/output_data/DynamicST_result/demo \
    --gtf_file xxx.gtf \
    --bam_file demo/output_data/Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam \
    --out_dir demo/output_data

Metabolic_Labeling_Correction.py： A binomial mixture model was used to correct the expression of new RNAs in spots of DynaST-seq data and estimate the substitution rate caused by metabolic RNA labeling.

 python Metabolic_Labeling_Correction.py \
    --mutations_file demo/output_data/all_conversions_list.csv \
    --reads_file demo/output_data/all_alignments_list.csv \
    --group_info_file demo/input_data/demo_group_info.csv \
    --convs_file demo/input_data/convs.pkl.gz \
    --out_dir demo/output_data
