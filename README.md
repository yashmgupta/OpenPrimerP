# OpenPrimerP: Targeted PCR Primer Design with GenBank Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)

## Abstract

OpenPrimerP is an intuitive web application designed to streamline PCR primer design by directly integrating with GenBank files. This tool addresses a significant gap in molecular biology research by automating the extraction of genomic features and applying advanced algorithms for primer design. Built on Streamlit and leveraging the Primer3 engine, OpenPrimerP provides a user-friendly interface that simplifies primer design for specific genomic features such as coding sequences, tRNAs, and genes. The application's specificity analysis feature further enhances experimental confidence by detecting potential off-target amplifications. OpenPrimerP significantly reduces the complexity of primer design, making sophisticated molecular techniques more accessible to researchers across experience levels.

## Features

- **Direct GenBank Integration**: Upload and parse GenBank files to extract sequence information
- **Multi-Record Support**: Process files containing multiple sequence records
- **Feature-Based Primer Design**: Target specific genomic features like CDS, tRNA, or genes
- **Customizable Parameters**: Adjust PCR product size range and primer counts based on experimental needs
- **Primer Specificity Analysis**: Detect potential off-target amplifications within the same genome
- **User-Friendly Interface**: Intuitive design with progress tracking and downloadable results
- **Comprehensive Output**: Detailed primer information including Tm values, lengths, and product sizes

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)

### Setup Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/yashmgupta/OpenPrimerP.git
   cd OpenPrimerP
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Start the application:
   ```bash
   streamlit run app.py
   ```

2. Access the application through your web browser (typically at http://localhost:8501)

3. Follow the step-by-step workflow:
   - Upload a GenBank (.gb or .gbk) file
   - Select a record if multiple are present
   - Choose a feature type (CDS, tRNA, gene)
   - Select a specific feature
   - Set PCR product size range and primer count
   - Click "Design Primers"
   - Run specificity analysis to check for potential off-target amplifications

## Technical Details

### Dependencies
- `streamlit`: Web application framework
- `primer3-py`: Python API for Primer3 engine
- `biopython`: Tools for biological computation
- `pandas`: Data manipulation and analysis

### Core Components
1. **GenBank Processing Module**: Extracts records and features from GenBank files
2. **Primer Design Engine**: Interfaces with Primer3 for optimal primer design
3. **Specificity Analysis**: Simulates PCR to detect off-target amplifications
4. **UI Components**: User-friendly interface built with Streamlit

## Research Applications

OpenPrimerP is designed to support various molecular biology applications:

- **Gene Expression Studies**: Design primers for RT-PCR and qPCR
- **Cloning Experiments**: Generate primers for amplifying genes of interest
- **Diagnostic PCR**: Create specific primers for pathogen detection
- **Genotyping**: Design primers for SNP detection and genotyping assays
- **Metagenomics**: Develop primers targeting specific taxonomic groups

## Citing OpenPrimerP

If you use OpenPrimerP in your research, please cite:

```
Gupta, Y.M. (2025). OpenPrimerP: An integrated application for feature-specific 
PCR primer design from GenBank files. [Publication details forthcoming]
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

Yash Munnalal Gupta - [yashmunnalalg@nu.ac.th](mailto:yashmunnalalg@nu.ac.th)

## Acknowledgments

- Primer3 developers for the primer design engine
- BioPython team for sequence analysis tools
- Streamlit team for the web application framework
