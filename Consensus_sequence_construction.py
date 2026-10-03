import argparse
import os
import subprocess
import pysam
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import array
from typing import Any, Dict, List, Optional

def read_gtf(file_path):
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf_data = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=columns)
    return gtf_data

def delete_temp_bam_files(temp_out_dir, num_batches):
    for i in range(num_batches):
        temp_bam = os.path.join(temp_out_dir, f"batch_{i}.bam")
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
            print(f"Deleted temporary BAM file: {temp_bam}")

def parse_attributes(attr_str):
    attributes = {}
    for item in attr_str.split(';'):
        if item.strip():
            key, value = item.strip().split(' ', 1)
            value = value.strip('"')
            attributes[key] = value
    return attributes

def build_gene_infos(gtf_file):
    gtf_data = read_gtf(gtf_file)
    genes = gtf_data[gtf_data["feature"] == "gene"]
    gene_infos = {}
    for index, row in genes.iterrows():
        attributes = parse_attributes(row["attribute"])
        gene_id = attributes.get("gene_id")
        if gene_id:
            gene_infos[gene_id] = {
                "gene_name": attributes.get("gene_name"),
                "gene_type": attributes.get("gene_type"),
                "chromosome": row["seqname"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "segment": (row["start"], row["end"])
            }
    return gene_infos

def parse_md_tag(md_tag):
    mismatches = []
    position = 0
    number = ''
    
    i = 0
    while i < len(md_tag):
        char = md_tag[i]
        
        if char.isdigit():
            number += char
        else:
            if number:
                position += int(number)
                number = ''
            
            if char == '^':
                i += 1
                insert_seq = ''
                while i < len(md_tag) and md_tag[i] in 'ACGT':
                    insert_seq += md_tag[i]
                    i += 1
                continue
            
            mismatches.append((position, char))
            position += 1
            
        i += 1

    if number:
        position += int(number)
        
    return mismatches

def process_batch(batch, gene_infos, header_dict, temp_out_dir, batch_num):
    temp_out = os.path.join(temp_out_dir, f"batch_{batch_num}.bam")
    with pysam.AlignmentFile(temp_out, 'wb', header=pysam.AlignmentHeader.from_dict(header_dict)) as out:
        for gene, barcode_data in batch:
            gene_info = gene_infos[gene]
            strand = gene_info['strand']
            for barcode, umi_data in barcode_data.items():
                for umi, reads in umi_data.items():
                    if len(reads) == 1:
                        read = reads[0]
                        read.set_tag('ST', strand)
                        out.write(read)
                    else:
                        tags = {}
                        tags['CB'] = barcode
                        tags['UB'] = umi
                        tags['GX'] = gene
                        tags['GN'] = gene_info['gene_name']
                        tags['ST'] = strand

                        Base = ('A', 'T', 'C', 'G')
                        Base_index = {b: i for i, b in enumerate(Base)}
                        left_pos = min([read.reference_start for read in reads])
                        right_pos = max([read.reference_end for read in reads])
                        length = right_pos - left_pos
                        seq = np.zeros((length, 4))
                        qua = np.zeros((length, 4))
                        ref = np.full(length, -1)
                        deletions = 0

                        for read in reads:
                            read_seq = read.query_sequence.upper()
                            read_qualities = read.query_qualities
                            for read_i, genome_i, genome_base in read.get_aligned_pairs(matches_only=False, with_seq=True):
                                if genome_i is not None and genome_base is not None:
                                    genome_base = genome_base.upper()
                                    if genome_base != 'N':
                                        i = genome_i - left_pos
                                        if read_i is None:
                                            if ref[i] == -1:
                                                ref[i] = Base_index[genome_base]
                                                deletions += 1
                                            continue

                                        read_base = read_seq[read_i]
                                        if read_base == 'N':
                                            continue

                                        if ref[i] == -1:
                                            ref[i] = Base_index[genome_base]

                                        seq[i, Base_index[read_base]] += 1

                        consensus_length = (seq > 0).any(axis=1).sum()
                        consensus = np.zeros(consensus_length, dtype=np.uint8)
                        frequence = np.zeros(consensus_length, dtype=np.uint8)

                        cigar = []
                        last_cigar_op = None
                        cigar_n = 0
                        md = []
                        md_n = 0
                        md_zero = True
                        md_del = False
                        nm = 0
                        consensus_i = 0

                        for i in range(length):
                            ref_base = ref[i]
                            cigar_op = 'N'

                            if ref_base >= 0:
                                seq_base = seq[i]

                                if (seq_base == 0).all():
                                    cigar_op = 'D'
                                    if md_n > 0 or md_zero:
                                        md.append(str(md_n))
                                        md_n = 0

                                    if not md_del:
                                        md.append('^')
                                    md.append(Base[ref_base])
                                    md_del = True

                                else:
                                    md_del = False

                                    base_q = qua[i].max()
                                    base_freq = seq_base.max()
                                    bases = (seq_base == base_freq).nonzero()[0]
                                    if len(bases) > 0 and ref_base in bases:
                                        base = ref_base
                                    else:
                                        base = bases[0]

                                    cigar_op = 'M'
                                    if ref_base == base:
                                        md_n += 1
                                        md_zero = False
                                    else:
                                        if md_n > 0 or md_zero:
                                            md.append(str(md_n))
                                            md_n = 0
                                        md.append(Base[ref_base])
                                        md_zero = True
                                        nm += 1

                                    consensus[consensus_i] = base
                                    frequence[consensus_i] = base_freq
                                    consensus_i += 1

                            if cigar_op == last_cigar_op:
                                cigar_n += 1
                            else:
                                if last_cigar_op:
                                    cigar.append(str(cigar_n) + last_cigar_op)
                                last_cigar_op = cigar_op
                                cigar_n = 1

                        md.append(str(md_n))
                        cigar.append(str(cigar_n) + last_cigar_op)

                        header = pysam.AlignmentHeader.from_dict(header_dict)
                        al = pysam.AlignedSegment(header)
                        al.query_name = sha256(''.join(read.query_name for read in reads).encode('utf-8')).hexdigest()
                        al.query_sequence = ''.join(Base[i] for i in consensus)
                        al.reference_name = reads[0].reference_name
                        al.reference_id = reads[0].reference_id
                        al.reference_start = left_pos
                        al.mapping_quality = 255
                        al.cigarstring = ''.join(cigar)

                        tags = tags or {}
                        tags.update({'MD': ''.join(md), 'NM': nm})
                        al.set_tags(list(tags.items()))

                        al.is_unmapped = False
                        al.is_paired = False
                        al.is_duplicate = False
                        al.is_qcfail = False
                        al.is_secondary = False
                        al.is_supplementary = False

                        out.write(al)

def merge_bams(output_file, temp_out_dir, num_batches):
    first_bam = os.path.join(temp_out_dir, f"batch_0.bam")
    with pysam.AlignmentFile(first_bam, 'rb') as template_bam:
        with pysam.AlignmentFile(output_file, 'wb', header=template_bam.header) as merged_bam:
            for i in range(num_batches):
                temp_bam = os.path.join(temp_out_dir, f"batch_{i}.bam")
                with pysam.AlignmentFile(temp_bam, 'rb') as temp:
                    for read in temp:
                        merged_bam.write(read)

def create_bam_index(bam_file):
    try:
        subprocess.check_call(['samtools', 'index', bam_file])
        print(f"Index for {bam_file} created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error creating index for {bam_file}: {e}")

def process_tags(input_bam, output_bam):
    def parse_md_tag(md_tag):
        mismatches = []
        position = 0
        number = ''
        
        i = 0
        while i < len(md_tag):
            char = md_tag[i]
            
            if char.isdigit():
                number += char
            else:
                if number:
                    position += int(number)
                    number = ''
                
                if char == '^':
                    i += 1
                    insert_seq = ''
                    while i < len(md_tag) and md_tag[i] in 'ACGT':
                        insert_seq += md_tag[i]
                        i += 1
                    continue
                
                mismatches.append((position, char))
                position += 1
                
            i += 1

        if number:
            position += int(number)
            
        return mismatches

    with pysam.AlignmentFile(input_bam, 'rb') as infile, pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
        for read in infile.fetch():
            md = read.get_tag('MD')
            strand = read.get_tag('ST')
            gene_name = read.get_tag('GN')

            if md.isdigit():
                gene_name += '--T'
                read.set_tag('GN', gene_name)
                outfile.write(read)
                continue

            if strand == '+' and 'T' not in md:
                gene_name += '--T'
                read.set_tag('GN', gene_name)
                outfile.write(read)
                continue

            if strand == '-' and 'A' not in md:
                gene_name += '--T'
                read.set_tag('GN', gene_name)
                outfile.write(read)
                continue

            if strand == '+' and 'T' in md:
                mismatches = parse_md_tag(md)
                found_match = False
                for position, ref_base in mismatches:
                    base = read.query_sequence[position]
                    if base == 'C' and ref_base == 'T':
                        gene_name += '--C'
                        read.set_tag('GN', gene_name)
                        found_match = True
                        break
                if not found_match:
                    gene_name += '--T'
                    read.set_tag('GN', gene_name)
            
            if strand == '-' and 'A' in md:
                mismatches = parse_md_tag(md)
                found_match = False
                for position, ref_base in mismatches:
                    base = read.query_sequence[position]
                    if base == 'G' and ref_base == 'A':
                        gene_name += '--C'
                        read.set_tag('GN', gene_name)
                        found_match = True
                        break
                if not found_match:
                    gene_name += '--T'
                    read.set_tag('GN', gene_name)

            outfile.write(read)
# def run_samtools_calmd(bam_file, ref_fa):
#     try:
#         subprocess.check_call(['samtools', 'calmd', '-AEur', bam_file, ref_fa, '>', f"{bam_file}_calmd.bam"])
#         print(f"Samtools calmd completed for {bam_file}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error during samtools calmd: {e}")
def run_samtools_calmd(bam_file, ref_fa):
    try:
        with open(f"{bam_file}_calmd.bam", "wb") as out_bam:
            process = subprocess.Popen(
                ['samtools', 'calmd', '-AEur', bam_file, ref_fa],
                stdout=out_bam,  # Redirect the output to the file
                stderr=subprocess.PIPE
            )
            _, stderr = process.communicate()

            if process.returncode != 0:
                print(f"Error during samtools calmd: {stderr.decode()}")
            else:
                print(f"Samtools calmd completed for {bam_file}.")
    except Exception as e:
        print(f"Error: {e}")
                
def main():
    parser = argparse.ArgumentParser(description="Process BAM files and GTF annotations.")
    parser.add_argument('--bam_path', required=True, help='Path to the input BAM file.')
    parser.add_argument('--gtf_file', type=str, default='Homo_sapiens.GRCh38.99.gtf', help='Path to the GTF file')
    parser.add_argument('--temp_out_dir', required=True, help='Temporary output directory for BAM files.')
    parser.add_argument('--ref_fa', required=True, help='Path to the reference genome in FASTA format.')  # Add ref_fa argument

    args = parser.parse_args()

    bam_path = args.bam_path
    ref_fa = args.ref_fa  # Get the reference genome path from the command-line argument
    gtf_file = args.gtf_file
    temp_out_dir = args.temp_out_dir
    # Run samtools calmd
    run_samtools_calmd(bam_path, ref_fa)
    create_bam_index(f"{bam_path}_calmd.bam")
    # After calmd, the BAM file will be updated and you can proceed to work with the calmd BAM file.
    bam_path = f"{bam_path}_calmd.bam"  # Update BAM path to the calmd file

    # The rest of your code follows...
    
    temp_out = os.path.join(temp_out_dir, 'Aligned.sortedByCoord.UniqueGene.consensus.bam')
    sorted_bam = os.path.join(temp_out_dir, 'Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam')
    tagged_bam = os.path.join(temp_out_dir, 'Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam')

    os.makedirs(temp_out_dir, exist_ok=True)

    gene_infos = build_gene_infos(gtf_file)

    with pysam.AlignmentFile(bam_path, 'rb') as f:
        header_dict = f.header.to_dict()

    barcode_umi_groups = {}
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        for read in f.fetch():
            if read.is_unmapped or read.is_secondary:
                continue
            barcode = read.get_tag('CB')
            umi = read.get_tag('UB')
            genes = read.get_tag('GX')
            if barcode == '-' or umi == '-' or genes == '-':
                continue
            barcode_umi_groups.setdefault(genes, {}).setdefault(barcode, {}).setdefault(umi, []).append(read)

    batches = [list(barcode_umi_groups.items())[i:i + 2000] for i in range(0, len(barcode_umi_groups), 2000)]

    with ThreadPoolExecutor(max_workers=50) as executor:
        futures = [executor.submit(process_batch, batch, gene_infos, header_dict, temp_out_dir, i) for i, batch in enumerate(batches)]
        for future in futures:
            future.result()

    merge_bams(temp_out, temp_out_dir, len(batches))
    create_bam_index(temp_out)
    pysam.sort("-o", sorted_bam, temp_out)
    create_bam_index(sorted_bam)

    print(f"Final sorted BAM file: {sorted_bam} with index created.")

    # Process tags
    process_tags(sorted_bam, tagged_bam)
    create_bam_index(tagged_bam)
    print(f"Tagged BAM file: {tagged_bam} created.")
    delete_temp_bam_files(temp_out_dir, len(batches))

if __name__ == "__main__":
    main()