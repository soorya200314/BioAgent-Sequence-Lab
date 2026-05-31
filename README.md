# BioAgent-Sequence-Lab

[![Open in Streamlit](https://static.streamlit.io/badge.svg)](https://bioagent-sequence-lab-5yjbyg7rsuvez7zfrypdwt.streamlit.app/)


An interactive genomic workbench built with Streamlit for optimal global sequence alignment, micro-mutational variant scanning, and macro-evolutionary phylogenetic inference.

## 🔍 Overview

This platform provides an end-to-end pipeline for analyzing multi-sequence DNA datasets. By pairing a robust computational biology engine with an intuitive user interface, users can upload raw genomic data to scan for point mutations and reconstruct evolutionary histories in real time.

## 🚀 Features

* **Global Sequence Alignment:** Uses the Needleman-Wunsch dynamic programming algorithm to establish optimal pairwise sequence alignments.
* **Micro-Mutational Analysis Dashboard:** Automatically scans query sequences against a selectable wild-type reference to track base-level variations, mapping mutations to exact sequence alignment positions, codon numbers, and triplet reading frames.
* **Macro-Evolutionary Inference:** Computes a comprehensive pairwise evolutionary distance matrix ($1.0 - \text{Sequence Identity}$) and reconstructs a phylogenetic tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) hierarchical clustering method.
* **Interactive Visualizations:** Dynamically plots publication-ready dendrograms and provides instant CSV export functionality for mutational reports.

## 🛠️ Tech Stack & Architecture

* **Frontend/Deployment:** Streamlit Cloud
* **Data Core:** Pandas, NumPy
* **Scientific/Visual Computing:** SciPy, Matplotlib

## 📦 Local Setup Instructions

If you wish to run this workbench locally on your machine, execute the following commands in your terminal:

```bash
# 1. Clone the repository
git clone [https://github.com/soorya200314/BioAgent-Sequence-Lab.git](https://github.com/soorya200314/BioAgent-Sequence-Lab.git)

# 2. Navigate to the project directory
cd BioAgent-Sequence-Lab

# 3. Install required mathematical and deployment libraries
pip install -r requirements.txt

## 📄 License & Usage

This project is open-source and available under the **MIT License**. Feel free to fork this repository, experiment with the alignment parameters, or adapt the UPGMA visualization modules for your own academic research.

> ⚠️ **Disclaimer:** *This tool is developed for computational biology research and educational exploration. For clinical or diagnostic pipelines, ensure validation against gold-standard benchmarking suites.*

# 4. Launch the local Streamlit server
streamlit run webapp.py
