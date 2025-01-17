# AtlasMaker

The goal of this repo will be to abstract methods from hiv-atlas-maker.py development to create species-agnostic method for reference-guided mRNA modeling.

Related repo: [AtlasMatcher](https://github.com/GenerGener/AtlasMatcher)

# Tutorial

## Overview

Atlas Maker "atlas-maker.py" is a Python tool designed for spliced RNA mapping, particularly useful for analyzing HIV genomic sequences. The tool implements a consensus-based approach using duplicate reads and model sequences to ensure high-quality mapping results.

## Prerequisites

- Python 3.x (tested in Python 3.9.19)
- BioPython[^1]
- minimap2[^2][^3] aligner installed and accessible in PATH
- Input files in FASTA format:
  - Input genome sequence
  - Reference sequence
  - Model set containing multiple sequences

## Installation

1. Download the `atlas-maker.py` script
2. Install required Python packages:
```bash
pip install biopython
```
3. Ensure minimap2 is installed and accessible in your PATH

## Basic Usage

```bash
python atlas-maker.py --input input.fa --reference ref.fa --model-set models.fa
```

### Optional Parameters

- `--similarity`: Minimum sequence similarity threshold (default: 0.8)
- `--coverage`: Minimum reference coverage threshold (default: 0.9)

## Input Files Format

### Input Genome Sequence (input.fa)
Single FASTA sequence representing the genome to be analyzed:
```
>sequence_name
ACTG...
```

### Reference Sequence (ref.fa)
Single FASTA sequence used as the mapping reference:
```
>reference_name
ACTG...
```

### Model Set (models.fa)
Multiple FASTA sequences representing different models:
```
>model1
ACTG...
>model2
ACTG...
```

## Workflow Details

1. **Quality Control**
   - Checks sequence similarity against reference
   - Verifies coverage requirements
   - Creates initial metadata file with QC metrics

2. **Input Processing**
   - Duplicates input sequence
   - Creates modified headers (_INPUT1, _INPUT2)
   - Stores in working directory

3. **Model Processing** (for each model)
   - Combines duplicated inputs with model sequence
   - Creates temporary combined FASTA file
   - Performs spliced alignment using minimap2
   - Sorts alignment results

4. **Consensus Generation**
   - Identifies regions where model aligns
   - Collects sequences from input reads in model regions
   - Requires 3x coverage (both inputs + model)
   - Generates consensus sequence
   - Creates output FASTA and metadata

## Background Figure 1: HIV-1 M NLAD8 direct native RNA
From the first full-length native/direct RNA sequencing of HIV-1. Salmon coloring is indicative of reads or sequences which match the forward orientation of the reference genome used. Visulized in IGV[^4]. See: [Gener, Alejandro R.](https://www.biorxiv.org/content/10.1101/845610v2), Kimata, Jason T. Full-coverage native RNA sequencing of HIV-1 viruses. bioRxiv. 2019/01/01. doi: 10.1101/845610.

<img width="1430" alt="Screenshot 2024-12-24 at 10 14 38 AM" src="https://github.com/user-attachments/assets/47ed3abb-2eda-4234-b83f-4a8ba50d8b40" />

## Background Figure 2: HIV-1 mRNA Atlas v1.0
From the first standardized HIV-1 mRNA Atlas. Visulized in IGV. See: [Gener AR](https://journals.lww.com/aidsonline/citation/2022/01010/anticipating_hiv_drug_resistance_with_appropriate.16.aspx). Anticipating HIV drug resistance with appropriate sequencing methods. AIDS. 2022 Jan 1;36(1):147-148. doi: 10.1097/QAD.0000000000003087. PMID: 34873093. 

<img width="1920" alt="Screenshot 2024-12-24 at 10 15 26 AM" src="https://github.com/user-attachments/assets/b3d1cb85-68c6-4f15-9b26-8d93ba721e6a" />

## Workflow Figures Figure 3: Reverse sandwich mapping performed for each model. 
Color-coding denotes distinct sequences which might coincidentaly share homology. For this process, it is important to ensure appropriate sequence is captured and fed forward into the output input-derived transcript models using HIV-1 M NL4-3 model set (HIV-1 mRNA Atlas v1.0). Neither reference genome sequence's nor transcript model's sequences make it into the output using the current implementation.

<img width="1084" alt="Screenshot 2024-12-24 at 9 40 57 AM" src="https://github.com/user-attachments/assets/d5c6a4d5-12b6-49e3-a31d-9c52daa3d84b" />


## Workflow Figures Figure 4: Example intermediate output (before consensus calling).
Blue coloring denotes features from the GenBank record MZ242719.1. Most are novel open reading frames elucidated from thresholding mRNAs seen in CD4 cells infected by suppernatant from pNL4-3 transfected cell lines. Visulized in IGV version 2.8.13 11/20/20 12:50 AM.

<img width="1920" alt="input vs reference with features" src="https://github.com/user-attachments/assets/b1fe88ae-f133-4b70-b977-ae48775ce4c8" />

## Output Files

The tool creates a timestamped working directory containing:

1. **QC Files**
   - `qc_metadata_{timestamp}.txt`: Quality control metrics

2. **Intermediate Files**
   - `input_1.fa`, `input_2.fa`: Duplicated input sequences
   - `combined_{model_id}.fa`: Combined sequences for each model
   - `alignment_{model_id}.sam`: Raw alignment results
   - `alignment_{model_id}_sorted.sam`: Sorted alignments

3. **Results**
   - `input_{ref_header}_MSTRG.2.{model_id}_consensus.fa`: Consensus sequences
   - `consensus_metadata_{model_id}_{timestamp}.txt`: Consensus metrics

## Example Usage

Using the provided example files:

```bash
python atlas-maker.py \
  --input MZ242719.fasta \
  --reference HIV-1_HXB2.fasta \
  --model-set model_set.fa
```

### Expected Output Structure

```
results_hiv-atlas-maker_{timestamp}/
├── qc_metadata_{timestamp}.txt
├── input_1.fa
├── input_2.fa
├── combined_MSTRG.2.4.fa
├── alignment_MSTRG.2.4.sam
├── alignment_MSTRG.2.4_sorted.sam
├── input_HXB2_MSTRG.2.4_consensus.fa
└── consensus_metadata_MSTRG.2.4_{timestamp}.txt
```

### Example Output Files

#### QC Metadata File (qc_metadata_{timestamp}.txt)
```
Sequence similarity: 91.23%
Reference coverage: 95.67%
Input sequence length: 9719
Reference sequence length: 9720
Number of matches: 8868
```

#### Duplicated Input Files

input_1.fa:
```
>MZ242719.1_INPUT1
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
```

input_2.fa:
```
>MZ242719.1_INPUT2
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
```

#### Combined Model File (combined_MSTRG.2.4.fa)
```
>MZ242719.1_INPUT1
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
>MZ242719.1_INPUT2
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
>MSTRG.2.4
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
```

#### Consensus Output (input_HXB2_MSTRG.2.4_consensus.fa)
```
>consensus_HXB2_MSTRG.2.4
GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGC
CTCAATAAAGCTTGCCTTGAGTGCTCAAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAG
ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACTTGAAAGCGA
```

#### Consensus Metadata (consensus_metadata_MSTRG.2.4_{timestamp}.txt)
```
Model ID: MSTRG.2.4
Model Length: 9719
Consensus Length: 9719
Length Difference: 0.00%
Positions with data: 9719
```

Note 9719 bp is the length of input genome of [HIV-1 M HXB2](https://www.ncbi.nlm.nih.gov/nuccore/1906382), in the form of LTR-HIV-LTR.

## Technical Details

### Minimap2 Parameters

The tool uses optimized parameters for spliced alignment:

```bash
minimap2 -ax splice \      # Enable spliced alignment
  -u n \                   # No long gaps between seed chains
  -G 500k \               # Max intron size
  -k 14 \                 # K-mer size
  -w 20 \                 # Minimizer window size
  --forward-only \        # Force forward strand
  -A 2 \                  # Match score
  -B 4 \                  # Mismatch penalty
  -O 4,24 \              # Gap open penalties
  -E 2,1                  # Gap extension penalties
```

(Much thanks to Heng Li and his work on minimap2.)

### Consensus Calling

The consensus sequence is generated with the following requirements:

- Both input copies must align (2x coverage minimum)
- Model sequence must align to the region
- Only bases with agreement between input copies are included
- Positions lacking model support are masked

### Example Alignment

The alignment output shows sequence comparisons between references and consensus sequences:

```
ref_MZ2427 ---------------------------------gggtctctctggttagaccagatctga
ref_HXB2   cctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctga
out_unsp   ---------------------------------gggtctctctggttagaccagatctga
out_19     ---------------------------------gggtctctctggttagaccagatctga
                                            ***************************

ref_MZ2427 gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
ref_HXB2   gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
out_unsp   gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
out_19     gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
           ************************************************************

ref_MZ2427 tgagtgctcaaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
ref_HXB2   tgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
out_unsp   tgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
out_19     tgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
           ********. **************************************************
```
(Pairwise alignment for tutorial example done with [MAFFT online server](https://mafft.cbrc.jp/alignment/software/) [^5]. Method FFT-NS-i (Standard) [^6]. Command: mafft --reorder --auto input.)

Key features visible in the alignment:
- Asterisks (*) indicate perfect conservation
- Dots (.) indicate positions with variation
- Dashes (-) indicate gaps or missing sequence
- Output sequences (out_unsp, out_19) show consensus with *input (*which happens to be the HIV-1 species-specific reference)

Note that because of U3-R ideosyncracies and how users might supply their input and reference genomes, .... The latest implentation works best if the user supplies a genome that is longer than the mRNA's which were used to inform HIV-1 M NL4-3 mRNA (HIV-1 Atlas V1.0). The unspliced HIV-1 M NL4-3 mRNA model [MZ242719.1](https://www.ncbi.nlm.nih.gov/nuccore/MZ242719.1/) is 9173 bases long. Other types of loci which might be amenable to this approach include repetitive elements, whether endogenous or exogenous (e.g., retroviruses), which leverage direct long terminal repeats in their propagation.

TODO: Double-check whether current script performs gap removal and/or requires linebreak removal.

## Troubleshooting

Common issues and solutions:

1. **Minimap2 Not Found**
   - Ensure minimap2 is installed and in your PATH
   - Test with `minimap2 --version`

2. **Low Coverage Results**
   - Check input sequence quality
   - Verify similarity to reference
   - Adjust coverage threshold if needed

3. **Empty Consensus**
   - Verify model sequence alignment
   - Check for sufficient input coverage
   - Examine alignment files for mapping issues

4. **Failed initiation**
   - Check similarity from QC. May need to adust --similarity if --input is an HIV-1 M non-B. Reminder that the HIV-1 Atlas v1.0 was made from HIV-1 M NL4-3 sequence. However, because of conservation across many primate immunodeficiency viruses, and because of flexibility of segment-based alignment which minimap2 -splice performs, the HIV-1 M NL4-3 model set sequence is sufficient for generating input-specific models.

## Best Practices

1. **Input Quality**
   - Use high-quality input sequences
   - HIV reference sets
      - [Los Alamos National Labs HIV Sequence Database Alignments](https://www.hiv.lanl.gov/content/sequence/NEWALIGN/align.html)
      - [HIV-1 in Nextstrain by LANL Group](https://nextstrain.org/groups/LANL-HIV-DB/HIV/genome)
      - [NCBI HIV Subtyping tool](https://www.ncbi.nlm.nih.gov/projects/genotyping/formpage.cgi)
      - [ViralZone Human Immunodeficiency Virus 1](https://viralzone.expasy.org/7)
   - Verify sequence format and content
   - Check for proper FASTA formatting

2. **Resource Management**
   - Monitor disk space for large datasets
   - Clean up intermediate files if needed
   - Use appropriate computing resources

3. **Output Validation**
   - Check consensus lengths against models
   - Verify coverage metrics
   - Examine alignment quality scores

## Performance Considerations

- Memory usage scales with sequence length
- Processing time depends on number of models
- Disk space needed for intermediate files
- Consider cleanup of temporary files for large datasets

[^1]: Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009 Jun 1;25(11):1422-3. doi: 10.1093/bioinformatics/btp163. Epub 2009 Mar 20. PMID: 19304878; PMCID: PMC2682512.
[^2]: Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996.
[^3]: Li H. New strategies to improve minimap2 alignment accuracy. Bioinformatics. 2021 Dec 7;37(23):4572-4574. doi: 10.1093/bioinformatics/btab705. PMID: 34623391; PMCID: PMC8652018.
[^4]: Thorvaldsdóttir H, Robinson JT, Mesirov JP. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. Brief Bioinform. 2013 Mar;14(2):178-92. doi: 10.1093/bib/bbs017. Epub 2012 Apr 19. PMID: 22517427; PMCID: PMC3603213.
[^5]: Kuraku S, Zmasek CM, Nishimura O, Katoh K. aLeaves facilitates on-demand exploration of metazoan gene family trees on MAFFT sequence alignment server with enhanced interactivity. Nucleic Acids Res. 2013 Jul;41(Web Server issue):W22-8. doi: 10.1093/nar/gkt389. Epub 2013 May 15. PMID: 23677614; PMCID: PMC3692103.
[^6]: Katoh K, Misawa K, Kuma K, Miyata T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res. 2002 Jul 15;30(14):3059-66. doi: 10.1093/nar/gkf436. PMID: 12136088; PMCID: PMC135756.






