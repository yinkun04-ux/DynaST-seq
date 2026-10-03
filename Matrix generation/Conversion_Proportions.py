import os
import pysam
import ngs_tools as ngs
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

def complement(base):
    """Returns the complementary base of a DNA nucleotide."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_map.get(base, base)

def process_bam_file(bam_file_path, barcodes, n_threads):
    """Processes a BAM file and returns mutations and alignments data."""
    conversions, alignments, conversion_counts = [], [], defaultdict(int)
    with pysam.AlignmentFile(bam_file_path, 'rb', threads=n_threads) as bam_file:
        contigs = bam_file.references
        valid_contigs = {contig for contig in contigs if any(True for _ in bam_file.fetch(contig, until_eof=True))}
        
        for contig in valid_contigs:
            for record in bam_file.fetch(contig):
                barcode, gene_id, umi, strand = record.get_tag("CB").replace("_", ""), record.get_tag("GX"), record.get_tag("UB"), record.get_tag('ST')
                if barcode != "-" and gene_id != "-" and umi != "-" and barcode in barcodes:
                    base_counts, mutations = Counter(record.get_reference_sequence().upper()), {}
                    read_sequence = record.query_sequence.upper()
                    
                    for read_index, ref_index, ref_base in record.get_aligned_pairs(matches_only=True, with_seq=True):
                        read_base = read_sequence[read_index]
                        ref_base = ref_base.upper()

                        if 'N' in (ref_base, read_base):
                            continue

                        conversion = f'{ref_base}{read_base}' if strand == '+' else f'{complement(ref_base)}{complement(read_base)}'
                        conversion_counts[conversion] += 1

                        if ref_base != read_base:
                            mutations[ref_index] = {'read_id': record.query_name, 'contig': contig, 'position': ref_index, 'conversion': conversion}
                    
                    conversions.extend(mutations.values())
                    alignments.append({
                        'read_id': record.query_name, 'barcode': barcode, 'umi': umi, 'gene_id': gene_id,
                        'A': base_counts.get("A", 0), 'C': base_counts.get("C", 0),
                        'G': base_counts.get("G", 0), 'T': base_counts.get("T", 0),
                        'status': 'unassigned', 'assigned': bool(gene_id)
                    })
    
    return conversions, alignments, conversion_counts

def delete_temp_bam_files(temp_out_dir):
    for filename in os.listdir(temp_out_dir):
        if filename.startswith('temp_') and (filename.endswith('.bam') or filename.endswith('.bam.bai')):
            file_path = os.path.join(temp_out_dir, filename)
            if os.path.exists(file_path):
                os.remove(file_path)
            
def main(bam_file, barcodes, n_threads, max_workers, out_dir, gtf_file):
    """Main function to process all BAM files and gather results."""
    # Initialize the result containers
    all_conversions_list = []
    all_alignments_list = []
    conversion_counts = defaultdict(int)

    os.makedirs(out_dir, exist_ok=True)
    gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(gtf_file, use_version=False)

    # Use the first item in bam_file list (bam_file[0])
    split_bam_res = [value[0] for value in ngs.bam.split_bam(bam_file[0], out_dir + "/temp", n=64, n_threads=n_threads, show_progress=True).values()]
    
    for bam_file in split_bam_res:
        pysam.index(bam_file, f'{bam_file}.bai', '-@', str(n_threads))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_bam_file, bam_file_path, barcodes, n_threads) for bam_file_path in split_bam_res]
        for future in futures:
            local_conversions, local_alignments, local_conversion_counts = future.result()
            all_conversions_list.extend(local_conversions)  # Add conversions from this future
            all_alignments_list.extend(local_alignments)  # Add alignments from this future
            for conversion, count in local_conversion_counts.items():
                conversion_counts[conversion] += count  # Accumulate conversion counts
    delete_temp_bam_files(out_dir)
    return all_conversions_list, all_alignments_list, conversion_counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM files and calculate conversion statistics.")
    parser.add_argument("--sample_path", help="Sample path for the directory")
    parser.add_argument("--bam_file", help="Path to the BAM file")
    parser.add_argument("--gtf_file", help="Path to the GTF file")
    parser.add_argument("--out_dir", help="Output directory")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of workers for parallel processing")
    parser.add_argument("--n_threads", type=int, default=4, help="Number of threads for BAM processing")
    
    args = parser.parse_args()

    # Extract arguments
    sample_path = args.sample_path
    bam_file = args.bam_file
    gtf_file = args.gtf_file
    out_dir = args.out_dir
    n_threads = args.n_threads
    max_workers = args.max_workers

    # Read barcodes
    barcodes = pd.read_csv(os.path.join(sample_path, 'outs/filtered_feature_bc_matrix/barcodes.tsv.gz'), header=None)[0].to_list()

    # Run the main function
    all_conversions_list, all_alignments_list, conversion_counts = main([bam_file], barcodes, n_threads, max_workers, out_dir, gtf_file)

    # Output results
    pd.DataFrame(all_conversions_list).to_csv(os.path.join(out_dir, 'all_conversions_list.csv'), index=False)
    pd.DataFrame(all_alignments_list).to_csv(os.path.join(out_dir, 'all_alignments_list.csv'), index=False)

    # Print conversion counts
    for conversion, count in conversion_counts.items():
        print(f"Conversion: {conversion}, Count: {count}")

    # Calculate total counts by start base
    start_base_totals = defaultdict(int)
    for conversion, count in conversion_counts.items():
        start_base_totals[conversion[0]] += count
    
    # Calculate conversion proportions
    conversion_proportions = {conversion: count / start_base_totals[conversion[0]] 
                            for conversion, count in conversion_counts.items() if conversion not in ['CC', 'AA', 'GG', 'TT']}
    
    labels = sorted(conversion_proportions.keys())
    proportions = [conversion_proportions[label] for label in labels]
    colors = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

    plt.figure(figsize=(8, 4))
    plt.bar(labels, proportions, color=colors)
    plt.xlabel('Conversion Type')
    plt.ylabel('Proportion')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "conversion_proportions.pdf"), bbox_inches='tight')
    plt.show()

    pd.DataFrame(list(conversion_proportions.items()), columns=['Conversion_Type', 'Proportion']).to_csv(os.path.join(out_dir, 'conversion_proportions.csv'), index=False)
    pd.DataFrame(list(conversion_counts.items()), columns=['Conversion_Type', 'Count']).to_csv(os.path.join(out_dir, 'conversion_counts.csv'), index=False)