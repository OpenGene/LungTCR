# LungTCR - TCR Repertoire Analysis and Lung Cancer Prediction

![GitHub](https://img.shields.io/badge/license-MIT-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)
![Python](https://img.shields.io/badge/Python-%3E%3D3.8-green)

## Overview

LungTCR (https://www.lungtcr.com/) is a comprehensive website for analyzing T-cell receptor (TCR) repertoire data with specialized functions for cancer risk assessment. This repository contains:

- TCR repertoire feature calculation pipeline
- Cancer-associated TCR enrichment scoring 
- Machine learning models for lung cancer/malignant pulmonary nodule risk prediction

## Key Features

| Feature | Description |
|---------|-------------|
| TCR Repertoire Analysis | Calculates 20+ diversity and clonality metrics |
| Cancer TCR Enrichment | Quantifies tumor-associated TCR signatures |
| Risk Prediction Models | Random Forest/GBM models for lung cancer risk assessment |
| Visualization | Plotting of key TCR features and model result |


## Input File Format

The input file format is VDJtools' table format. 
Run Convert routine by VDJtools (https://vdjtools-doc.readthedocs.io/en/master/input.html#vdjtools-format) to geneate the format.

### Clonotype Table Requirements

| Column | Required | Description | Example |
|--------|----------|-------------|---------|
| count | Yes | Read counts of TCR clones | 161853 |
| freq | Yes | Frequency of TCR clones | 0.009385 |
| cdr3_nt | Yes | CDR3 nucleic acid sequence  | TGTGCCAGTTCGTCGTCTAGCTCCTACAATGAGCAGTTCTTC |
| cdr3_aa | Yes | CDR3 amino acid sequence | CASSSSSSYNEQFF |
| v | Yes | V gene segment | TRBV6-4 |
| d | Yes | D gene segment | . |
| j | Yes | J gene segment | TRBJ2-7 |
| VEnd | No | Position of the V gene end | 7 |
| Dstart | No | Position of the d gene start | . |
| Dend | No | Position of the D gene end | . |
| Jstart | No | Position of the J gene strat | 18 |
| sample_id | No | Sample identifier | Patient01_PBMC |


Example file (tabular format):
```
count	freq	cdr3nt	cdr3aa	v	d	j	VEnd	DStart	DEnd	JStart
161853	0.009385105218213133	TGTGCCAGTTCGTCGTCTAGCTCCTACAATGAGCAGTTCTTC	CASSSSSSYNEQFF	TRBV6-4	.	TRBJ2-1	7	-1	-1	18
128851	0.007471472215355789	TGTGCCAGCTCACCATAGGACAGTGCTTCTCTGGAAACACCATATATTTT	CASSP*DS_FSGNTIYF	TRBV18	TRBD1	TRBJ1-3	14	17	22	28
107730	0.006246763329429179	TGTGCCAGCAGTTACGGTCTAAGAGATACGCAGTATTTT	CASSYGLRDTQYF	TRBV6-5	.	TRBJ2-3	14	-1	-1	23
...
```

## Basic Usage

```bash
python TCRfeatureCal.py -m /extdata/metadata.tsv -o output_features/
python diversity_vdjtools_wrapper.py -m /extdata/metadata.txt -o outout_diversity/ -x 10000000
Rscript LungTCR_ModelPrediction.R -i /extdata/test.csv -o output_model/ -m /model/models_list.rds -f /model/SelectedFeatures.csv 
```

Output includes:
- Diversity indices (Shannon, Simpson, etc.)
- Clonality metrics
- V/J gene usage profiles
- CDR3 length distributions
- CDR3 amino acid compositions
- TCR clones frequency distributions
- TCR convergence index
- Lung cancer enrichment score
- Lung cancer prediction result





## Citation

The code files here are linked to the work "Large-Scale TCR Repertoire Profiling Unveils Tumor-Specific Signals for Diagnosing Indeterminate Pulmonary Nodules" by Chen et al .
