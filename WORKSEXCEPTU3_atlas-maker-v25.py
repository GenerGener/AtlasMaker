#!/usr/bin/env python3

import os
import sys
import datetime
import argparse
from pathlib import Path
import subprocess
import re
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import shutil

class AtlasMaker:
    def __init__(self, input_file, reference_file, model_set_file, 
                 similarity_threshold=0.8, coverage_threshold=0.9):
        self.input_file = input_file
        self.reference_file = reference_file
        self.model_set_file = model_set_file
        self.similarity_threshold = similarity_threshold
        self.coverage_threshold = coverage_threshold
        self.timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        self.working_dir = f"results_hiv-atlas-maker_{self.timestamp}"
        
    def setup_directories(self):
        try:
            os.makedirs(self.working_dir, exist_ok=True)
            print(f"Created working directory: {self.working_dir}")
            return True
        except Exception as e:
            print(f"Error creating directory: {e}")
            return False

    def check_dependencies(self):
        """Verify minimap2 is installed and accessible"""
        try:
            result = subprocess.run(['minimap2', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                print(f"Found minimap2 version: {result.stdout.strip()}")
                return True
            else:
                print("Error: minimap2 not found. Please ensure it's installed and in PATH")
                return False
        except Exception as e:
            print(f"Error checking minimap2: {e}")
            return False

    def check_sequence_quality(self):
        try:
            input_records = list(SeqIO.parse(self.input_file, "fasta"))
            ref_records = list(SeqIO.parse(self.reference_file, "fasta"))
            
            if not input_records or not ref_records:
                print("Error: Empty input or reference file")
                return False
            
            input_seq = str(input_records[0].seq)
            ref_seq = str(ref_records[0].seq)
            
            matches = sum(a == b for a, b in zip(input_seq[:min(len(input_seq), len(ref_seq))], 
                                               ref_seq[:min(len(input_seq), len(ref_seq))]))
            similarity = matches / min(len(input_seq), len(ref_seq))
            coverage = len(input_seq) / len(ref_seq)
            
            metadata_file = os.path.join(self.working_dir, f"qc_metadata_{self.timestamp}.txt")
            with open(metadata_file, "w") as f:
                f.write(f"Sequence similarity: {similarity:.2%}\n")
                f.write(f"Reference coverage: {coverage:.2%}\n")
                f.write(f"Input sequence length: {len(input_seq)}\n")
                f.write(f"Reference sequence length: {len(ref_seq)}\n")
                f.write(f"Number of matches: {matches}\n")
            
            if similarity >= self.similarity_threshold and coverage >= self.coverage_threshold:
                print(f"QC passed: Similarity {similarity:.2%}, Coverage {coverage:.2%}")
                return True
            else:
                print(f"QC failed: Similarity {similarity:.2%}, Coverage {coverage:.2%}")
                return False
                
        except Exception as e:
            print(f"Error in sequence quality check: {e}")
            return False

    def duplicate_input(self):
        try:
            input_1 = os.path.join(self.working_dir, "input_1.fa")
            input_2 = os.path.join(self.working_dir, "input_2.fa")
            
            # Get original ID from input file
            record = next(SeqIO.parse(self.input_file, "fasta"))
            input_id = record.id
            
            # Write duplicates with modified IDs that will be easy to track
            for idx, dest_file in enumerate([input_1, input_2], 1):
                with open(dest_file, 'w') as f:
                    f.write(f">{input_id}_INPUT{idx}\n{str(record.seq)}\n")
            
            print(f"Created duplicate input files in {self.working_dir}")
            return True
        except Exception as e:
            print(f"Error duplicating input: {e}")
            return False

    def create_combined_fasta(self, model_id):
        """Combine duplicated input files with model sequence"""
        try:
            combined_file = os.path.join(self.working_dir, f"combined_{model_id}.fa")
            input_1 = os.path.join(self.working_dir, "input_1.fa")
            input_2 = os.path.join(self.working_dir, "input_2.fa")
            temp_model = os.path.join(self.working_dir, f"temp_model_{model_id}.fa")
            
            # Write model sequence to temp file
            with open(self.model_set_file) as f:
                for record in SeqIO.parse(f, "fasta"):
                    if record.id == model_id:
                        with open(temp_model, 'w') as model_out:
                            SeqIO.write(record, model_out, "fasta")
                        break
            
            # Concatenate files
            with open(combined_file, 'wb') as outfile:
                for infile in [input_1, input_2, temp_model]:
                    with open(infile, 'rb') as readfile:
                        shutil.copyfileobj(readfile, outfile)
            
            os.remove(temp_model)
            print(f"Created combined FASTA for model {model_id}")
            return combined_file
            
        except Exception as e:
            print(f"Error creating combined FASTA: {e}")
            return None

    def run_minimap2_alignment(self, query, reference, output):
        try:
            if not os.path.exists(query):
                print(f"Error: Query file not found: {query}")
                return False
            if not os.path.exists(reference):
                print(f"Error: Reference file not found: {reference}")
                return False
            
            cmd = ["minimap2"]
            cmd.extend([
                "-ax", "splice",     # Enable spliced alignment
                "-u", "n",           # No long gaps between seed chains
                "-G", "500k",        # Max intron size
                "-k", "14",          # K-mer size
                "-w", "20",          # Minimizer window size
                "--forward-only",    # Force forward strand only (replaces -F 0)
                "-A", "2",           # Match score
                "-B", "4",           # Mismatch penalty
                "-O", "4,24",        # Gap open penalties
                "-E", "2,1",         # Gap extension penalties
                reference,
                query
            ])
            
            print(f"Running minimap2 command: {' '.join(cmd)}")
            
            with open(output, 'w') as out_file:
                result = subprocess.run(
                    cmd,
                    stdout=out_file,
                    stderr=subprocess.PIPE,
                    text=True
                )
            
            if result.returncode != 0:
                print("Minimap2 execution failed:")
                print(f"Exit code: {result.returncode}")
                print("Error output:")
                print(result.stderr)
                return False
            
            if not os.path.exists(output):
                print(f"Error: Minimap2 did not generate output file: {output}")
                return False
            
            if os.path.getsize(output) == 0:
                print(f"Error: Minimap2 generated empty output file: {output}")
                print("Minimap2 stderr output:")
                print(result.stderr)
                return False
            
            print(f"Successfully completed alignment: {output}")
            print(f"Alignment file size: {os.path.getsize(output)} bytes")
            return True
            
        except Exception as e:
            print(f"Error running minimap2: {e}")
            return False

    def sort_sam_file(self, input_sam, output_sam=None):
        """Sort SAM file"""
        try:
            if output_sam is None:
                output_sam = os.path.join(
                    self.working_dir,
                    os.path.basename(input_sam).replace('.sam', '_sorted.sam')
                )
            
            print(f"Sorting SAM file:")
            print(f"Input: {input_sam}")
            print(f"Output: {output_sam}")
            
            # Read header lines and alignments separately
            header_lines = []
            alignments = []
            with open(input_sam, 'r') as f:
                # Collect header lines
                for line in f:
                    if line.startswith('@'):
                        header_lines.append(line)
                    else:
                        alignments.append(line)
            
            # Sort alignments based on reference position
            sorted_alignments = sorted(alignments, 
                                     key=lambda x: (x.split('\t')[2], # reference name
                                                  int(x.split('\t')[3]) # position
                                                 ) if len(x.split('\t')) > 3 else ("", 0))
            
            # Write sorted SAM
            with open(output_sam, 'w') as out:
                # Write header
                for line in header_lines:
                    out.write(line)
                # Write sorted alignments
                for aln in sorted_alignments:
                    out.write(aln)
            
            print(f"Successfully sorted SAM file: {output_sam}")
            return output_sam
            
        except Exception as e:
            print(f"Error in sort_sam_file: {str(e)}")
            return None

    def parse_sam_alignment(self, sam_file, model_length):
        """Parse SAM file to extract consensus sequence, masking regions without model coverage"""
        try:
            # Get original input ID to match modified headers
            record = next(SeqIO.parse(self.input_file, "fasta"))
            input_base_id = record.id
            
            # Track alignments by source
            model_regions = set()  # Regions where model aligns
            input_data = defaultdict(lambda: {'bases': [], 'depth': 0, 'model_supported': False})
            
            print("\nProcessing alignments...")
            print("Reading SAM file...")

            # First pass - identify model regions
            with open(sam_file, 'r') as f:
                for line in f:
                    if line.startswith('@'):  # Skip header lines
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) < 11:  # Skip malformed lines
                        continue
                        
                    query_name = fields[0]
                    ref_pos = int(fields[3])
                    cigar = fields[5]
                    
                    if cigar == '*':  # Skip unmapped reads
                        continue

                    # Process model alignments first
                    if query_name == self.model_id:
                        pos = ref_pos
                        cigar_nums = re.split(r'[MIDNSHP=X]', cigar)[:-1]
                        cigar_ops = re.split(r'[0-9]+', cigar)[1:]
                        cigar_pairs = zip(map(int, cigar_nums), cigar_ops)
                        
                        for length, op in cigar_pairs:
                            if op == 'M':  # Match or mismatch
                                for i in range(pos, pos + length):
                                    model_regions.add(i)
                                pos += length
                            elif op in 'DN':  # Deletion or skip
                                pos += length

            # Second pass - process input alignments
            with open(sam_file, 'r') as f:
                for line in f:
                    if line.startswith('@') or len(line.strip().split('\t')) < 11:
                        continue
                        
                    fields = line.strip().split('\t')
                    query_name = fields[0]
                    ref_pos = int(fields[3])
                    cigar = fields[5]
                    seq = fields[9]
                    
                    if cigar == '*':
                        continue

                    # Process input reads
                    if query_name.startswith(f"{input_base_id}_INPUT"):
                        pos = ref_pos
                        seq_pos = 0
                        cigar_nums = re.split(r'[MIDNSHP=X]', cigar)[:-1]
                        cigar_ops = re.split(r'[0-9]+', cigar)[1:]
                        cigar_pairs = zip(map(int, cigar_nums), cigar_ops)
                        
                        for length, op in cigar_pairs:
                            if op == 'M':  # Match or mismatch
                                for i in range(length):
                                    curr_pos = pos + i
                                    if curr_pos in model_regions:  # Only collect where model aligns
                                        input_data[curr_pos]['bases'].append(seq[seq_pos + i])
                                        input_data[curr_pos]['depth'] = len(input_data[curr_pos]['bases'])
                                        input_data[curr_pos]['model_supported'] = True
                                pos += length
                                seq_pos += length
                            elif op in 'DN':  # Deletion or skip
                                pos += length
                            elif op in 'IS':  # Insertion or soft clip
                                seq_pos += length
            
            # Clear data from positions without model support
            for pos in list(input_data.keys()):
                if not input_data[pos]['model_supported']:
                    del input_data[pos]
            
            print(f"Model regions: {len(model_regions)}")
            print(f"Positions with input data (model-supported): {len(input_data)}")
            if len(input_data) > 0:
                print("Sample of input data:")
                for pos in sorted(list(input_data.keys()))[:3]:
                    print(f"Position {pos}: depth={input_data[pos]['depth']}, bases={input_data[pos]['bases']}, model_supported={input_data[pos]['model_supported']}")
            
            return input_data, sorted(model_regions)
            
        except Exception as e:
            print(f"Error parsing SAM file: {e}")
            return None, None

    def call_consensus(self, alignment_data, model_regions):
        """Call consensus requiring 3x coverage (model + both inputs)"""
        consensus_seq = []
        
        for pos in sorted(model_regions):
            if pos in alignment_data:
                bases = alignment_data[pos]['bases']
                depth = alignment_data[pos]['depth']
                
                if depth == 2:  # Must have both input copies
                    base_counts = defaultdict(int)
                    for base in bases:
                        base_counts[base] += 1
                    
                    if len(base_counts) == 1:  # Both inputs agree
                        consensus_seq.append(list(base_counts.keys())[0])
                    elif len(base_counts) > 1:
                        # Take majority if exists
                        max_base = max(base_counts.items(), key=lambda x: x[1])
                        if max_base[1] > 1:
                            consensus_seq.append(max_base[0])
        
        return ''.join(consensus_seq)
    
    def generate_consensus(self, sorted_sam, model_id):
        try:
            self.model_id = model_id
            model_length = None
            for record in SeqIO.parse(self.model_set_file, "fasta"):
                if record.id == model_id:
                    model_length = len(record.seq)
                    break
            
            if model_length is None:
                print(f"Error: Could not find model {model_id}")
                return False
            
            alignment_data, model_regions = self.parse_sam_alignment(sorted_sam, model_length)
            if not alignment_data:
                print("No alignment data found")
                return False
            
            consensus_sequence = self.call_consensus(alignment_data, model_regions)
            if not consensus_sequence:
                print(f"Error: No consensus sequence generated for model {model_id}")
                return False
            
            length_diff = abs(len(consensus_sequence) - model_length)
            length_diff_percent = (length_diff / model_length) * 100
            
            ref_header = None
            with open(self.reference_file, 'r') as ref_file:
                ref_header = ref_file.readline().strip().split()[0][1:]
            
            output_file = os.path.join(
                self.working_dir,
                f"input_{ref_header}_MSTRG.2.{model_id}_consensus.fa"
            )
            with open(output_file, 'w') as out:
                header = f">consensus_{ref_header}_MSTRG.2.{model_id}"
                out.write(f"{header}\n{consensus_sequence}\n")
            
            metadata_file = os.path.join(
                self.working_dir,
                f"consensus_metadata_{model_id}_{self.timestamp}.txt"
            )
            with open(metadata_file, 'w') as meta:
                meta.write(f"Model ID: {model_id}\n")
                meta.write(f"Model Length: {model_length}\n")
                meta.write(f"Consensus Length: {len(consensus_sequence)}\n")
                meta.write(f"Length Difference: {length_diff_percent:.2f}%\n")
                positions_with_data = sum(1 for pos in model_regions if pos in alignment_data)
                meta.write(f"Positions with data: {positions_with_data}\n")
            
            return True
            
        except Exception as e:
            print(f"Error generating consensus for model {model_id}: {e}")
            return False

    def run_pipeline(self):
        print("Starting Atlas Maker pipeline...")
        
        if not self.check_dependencies():
            return False
            
        if not self.setup_directories():
            return False
            
        if not self.check_sequence_quality():
            return False
            
        if not self.duplicate_input():
            return False
            
        for model in SeqIO.parse(self.model_set_file, "fasta"):
            model_id = model.id
            print(f"\nProcessing model: {model_id}")
            
            combined_file = self.create_combined_fasta(model_id)
            if not combined_file:
                continue
            
            alignment_output = os.path.join(self.working_dir, f"alignment_{model_id}.sam")
            if not self.run_minimap2_alignment(
                combined_file,
                self.reference_file,
                alignment_output
            ):
                continue
            
            sorted_sam = self.sort_sam_file(alignment_output)
            if not sorted_sam:
                continue
            
            if not self.generate_consensus(sorted_sam, model_id):
                continue
            
            if os.path.exists(combined_file):
                os.remove(combined_file)
            
        print("Atlas Maker pipeline completed")
        return True

def main():
    parser = argparse.ArgumentParser(description='Atlas Maker v22 - HIV RNA Model Prediction Tool')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--reference', required=True, help='Reference sequence file')
    parser.add_argument('--model-set', required=True, help='Model set FASTA file')
    parser.add_argument('--similarity', type=float, default=0.8,
                      help='Minimum sequence similarity threshold (default: 0.8)')
    parser.add_argument('--coverage', type=float, default=0.9,
                      help='Minimum reference coverage threshold (default: 0.9)')
    
    args = parser.parse_args()
    
    for file_path in [args.input, args.reference, args.model_set]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}")
            sys.exit(1)
    
    atlas_maker = AtlasMaker(
        input_file=args.input,
        reference_file=args.reference,
        model_set_file=args.model_set,
        similarity_threshold=args.similarity,
        coverage_threshold=args.coverage
    )
    
    success = atlas_maker.run_pipeline()
    if success:
        print("Atlas Maker completed successfully")
        sys.exit(0)
    else:
        print("Atlas Maker failed")
        sys.exit(1)

if __name__ == "__main__":
    main()