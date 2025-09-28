# OpenPrimerP

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://openprimerp.streamlit.app/)

An open-source application for PCR Primer Design with GenBank integration.

![OpenPrimerP Banner](https://github.com/yashmgupta/OpenPrimerP/raw/main/assets/banner.png)

## 🧬 Overview

OpenPrimerP is a user-friendly web application designed to streamline the process of PCR primer design by integrating with GenBank files. The tool allows researchers to efficiently design primers targeting specific genomic features such as coding sequences (CDS), tRNAs, and genes directly from GenBank annotations, addressing a critical gap in the current landscape of molecular biology research tools.

## 🔍 Features

- **GenBank Integration**: Direct upload and processing of GenBank (.gb, .gbk) files
- **Multi-Record Support**: Design primers from any record within multi-sequence GenBank files
- **Feature-Specific Design**: Target specific genomic features (CDS, tRNA, genes)
- **Customizable Parameters**: Control product size range and number of primer pairs
- **Primer Specificity Analysis**: Check primers against other sequences to avoid off-target amplification
- **Interactive Interface**: User-friendly web interface powered by Streamlit
- **Results Export**: Download primer designs in CSV format

## 🚀 Quickstart

Visit the hosted application at: [https://openprimerp.streamlit.app/](https://openprimerp.streamlit.app/)

## 📋 Usage Guide

1. **Upload a GenBank File**
   - Click the upload area to select a GenBank (.gb or .gbk) file
   - The application supports both single and multi-record GenBank files

2. **Select a Record** (for multi-record files)
   - Choose the specific record you want to work with from the dropdown menu

3. **Select a Feature**
   - Choose the feature type (CDS, tRNA, gene)
   - Select the specific feature from the list

4. **Set Parameters**
   - Enter the desired PCR product size range (e.g., 150-500 bp)
   - Specify the minimum number of primer pairs to generate

5. **Design Primers**
   - Click the "Design Primers" button to generate primer pairs
   - View the results in a table showing sequences, melting temperatures, and product sizes
   - Download the results as a CSV file

6. **Check Primer Specificity** (optional)
   - Select a primer pair to analyze
   - Run the specificity analysis to check for potential off-target amplification
   - View detailed results and recommendations

## 🔧 Installation for Local Development

```bash
# Clone the repository
git clone https://github.com/yashmgupta/OpenPrimerP.git

# Navigate to the project directory
cd OpenPrimerP

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run app.py
```

## 📚 Dependencies

- streamlit
- pandas
- primer3-py
- biopython
- numpy

## 🧪 Technical Overview

OpenPrimerP leverages several key technologies and algorithms:

1. **Primer3**: For core primer design calculations including Tm, GC content optimization, and self-complementarity checks
2. **BioPython**: For parsing and manipulating GenBank files and sequence data
3. **Streamlit**: For the web-based interactive user interface
4. **PCR Simulation**: Custom algorithms to predict potential amplification products

The application workflow follows these steps:
1. GenBank file parsing and feature extraction
2. Target sequence identification and extraction
3. Primer design using Primer3 algorithms
4. Specificity checking against other sequences in the file
5. Results presentation and interpretation

## 📝 Research Context

This application was developed as part of a research project focused on improving accessibility and efficiency in molecular biology workflows. The tool addresses challenges in primer design for targeted gene amplification by automating the extraction and analysis of genomic features from GenBank files.

Key innovations include:
- Integration of GenBank parsing with primer design algorithms
- Feature-specific targeting with automated sequence extraction
- Cross-reference specificity checking for multi-sequence files
- Research-oriented output formatting for publication-ready results

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📊 Citation Information

If you use OpenPrimerP in your research, please cite:

```
Gupta, Y.M. (2025). OpenPrimerP: Application for Primer Design for Targeted Gene with GenBank Integration. 
DOI: 10.xxxx/xxxxx
```

## 📞 Contact

Yash Munnalal Gupta - [yashmunnalalg@nu.ac.th](mailto:yashmunnalalg@nu.ac.th)

Project Link: [https://github.com/yashmgupta/OpenPrimerP](https://github.com/yashmgupta/OpenPrimerP)

Live Demo: [https://openprimerp.streamlit.app/](https://openprimerp.streamlit.app/)
