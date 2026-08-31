# Telomere Motif Pipeline
A computational workflow for identifying and characterizing telomeric repeat sequences in plant genomes. The pipeline detects candidate motifs, particularly the canonical plant repeat TTTAGGG, and analyzes their abundance, genomic distribution, and presence at chromosome or scaffold ends using Inhouse Python pipeline for extracting chromosome termini and quantifying telomeric repeat motifs.

## Features

- Extract first 10 kb and last 10 kb from every chromosome
- Count overlapping motif occurrences
- Support multiple motifs
- Generate Excel reports
- Generate TSV reports
- Generate motif position files
- Create publication-quality plots
- Calculate motif density
- Summarize telomere enrichment

---
## Requirement
biopython>=1.81
pandas>=2.0
numpy>=1.24
matplotlib>=3.7
openpyxl>=3.1
seaborn >=0.13.2


## Installation

bash
git clone https://github.com/Naveen-kumar-Pandey/Telomare_finder.git

cd Telomare_finder

pip install -r requirements.txt


---

## Input

Genome FASTA:

fasta \
>Chr1\
ATGC...\
>Chr2\
ATGC...


---

## Run

bash \
python telomere_pipeline.py \
-i genome.fa \
-m TTTAGGG,CCCTAAA


---

## Output

### motif_counts.xlsx

| Chromosome | Region | Motif | Count |
|------------|---------|--------|--------|
| Chr1 | Start_10kb | TTTAGGG | 542 |
| Chr1 | End_10kb | TTTAGGG | 681 |

### motif_positions.tsv

Coordinates of every motif occurrence.

### chromosome_starts_10kb.fa

First 10 kb of every chromosome.

### chromosome_ends_10kb.fa

Last 10 kb of every chromosome.

### plots/

Contains:

- motif_counts.png
- motif_heatmap.png
- motif_distribution.png

---

## Example

#bash

python telomere_pipeline.py \
-i Arabidopsis.fa \
-m TTTAGGG,CCCTAAA


---

## Citation

If you use this pipeline in your research, please cite the repository.
