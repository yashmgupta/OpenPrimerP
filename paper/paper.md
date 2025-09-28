---
title: 'OpenPrimerP: An Intuitive Application for PCR Primer Design from GenBank Sequences'
tags:
  - Python
  - PCR primers
  - GenBank
  - primer design
  - biodiversity
  - molecular biology
  - bioinformatics
  - Streamlit
authors:
  - name: Yash Munnalal Gupta
    orcid: 0000-0003-3306-832X
    corresponding: true
    affiliation: 1
affiliations:
 - name: Department of Biology, Faculty of Science, Naresuan University, Thailand
   index: 1
date: 28 September 2024
bibliography: paper.bib
---

# Summary

Molecular identification and phylogenetic analysis of biodiversity rely heavily on PCR amplification of target genes, yet primer design for taxonomic studies remains technically challenging for many researchers. `OpenPrimerP` addresses this barrier by providing an intuitive web application that automatically designs PCR primers for specific genomic features directly from GenBank sequences. This tool integrates Primer3 algorithms [@Untergasser2012] with automated GenBank file processing to enable rapid primer development for coding sequences, tRNAs, and genes commonly used in biodiversity studies such as COI, 16S rRNA, and ITS regions. The application significantly reduces the technical expertise required for primer design while maintaining high specificity and efficiency standards essential for reliable taxonomic identification. Validation demonstrates `OpenPrimerP`'s ability to generate robust primers suitable for species identification, population genetics, and phylogenetic reconstruction across diverse taxa.

# Statement of need

The design of effective PCR primers remains a significant technical barrier for many biodiversity researchers, particularly field biologists and taxonomists who may lack extensive bioinformatics training. While sophisticated primer design tools exist, including Primer3 [@Untergasser2012], BatchPrimer3 [@You2008], and Primer-BLAST [@Ye2012], these often require considerable bioinformatics expertise and do not provide seamless integration with GenBank file processing.

A significant research gap persists in providing an accessible and user-friendly solution tailored specifically for designing primers within GenBank files. The process of identifying and extracting specific genomic features from GenBank files, such as coding sequences, tRNAs, and genes, can be cumbersome and time-consuming, hindering the efficiency of primer design workflows.

`OpenPrimerP` fills this gap by democratizing molecular techniques in biodiversity research. The tool is specifically designed with biodiversity applications in mind, incorporating features essential for taxonomic studies and enabling researchers to rapidly develop primers for species identification, population studies, and phylogenetic analysis directly from GenBank sequences without requiring extensive bioinformatics knowledge.

# Methods

The development methodology for `OpenPrimerP` followed a systematic software engineering approach specifically designed to address the needs of biodiversity researchers requiring efficient primer design capabilities. The foundational algorithm selection process prioritized the widely-recognized Primer3 software [@Untergasser2012] as the core computational engine due to its comprehensive thermodynamic modeling capabilities and proven reliability in generating specific and efficient primers essential for taxonomic and phylogenetic studies.

The technical architecture was developed through a modular approach that integrated the Primer3 library directly into a custom Python-based application framework. A specialized GenBank file processing module was engineered to automatically parse and extract genomic features commonly utilized in biodiversity research, including coding sequences (CDS), transfer RNAs (tRNAs), and target genes frequently used in molecular taxonomy such as COI, 16S rRNA, and ITS regions. This processing capability was implemented using the Biopython library [@Chapman2000] to ensure robust handling of diverse GenBank file formats and annotation standards encountered in taxonomic databases.

The user interface was constructed using the Streamlit framework to provide an accessible web-based environment that eliminates the need for specialized bioinformatics training. The integrated primer design workflow combines automated feature extraction with parameter customization options that allow researchers to specify PCR product size ranges appropriate for their experimental protocols and specify the minimum number of primer pairs required for robust amplification strategies.

Validation procedures involved extensive testing across GenBank files representing diverse taxonomic groups including insects (COI genes), fungi (ITS regions), and bacteria (16S rRNA genes) to ensure broad applicability and reliability for biodiversity research applications.

# Results

`OpenPrimerP` provides an intuitive interface where users can upload GenBank files containing single or multiple records, select desired feature types (CDS, tRNA, or gene), and choose specific features from available options. Users can specify desired PCR product size ranges and the minimum number of primer pairs to generate. The application processes the input through the following workflow:

1. **GenBank file parsing**: Automatic extraction of genomic features using Biopython
2. **Feature selection**: User-friendly interface for selecting target sequences
3. **Primer design**: Integration with Primer3 for robust primer generation
4. **Results presentation**: Interactive tables with comprehensive primer information
5. **Data export**: CSV format for integration with laboratory workflows

The application generates comprehensive output data including primer sequences, melting temperatures, primer lengths, and predicted PCR product sizes, presented in both interactive tabular format and exportable CSV files for integration with laboratory information management systems.

Validation testing demonstrated successful primer design across diverse taxonomic groups, confirming the tool's broad applicability for biodiversity research. The web application is freely accessible at https://OpenPrimerP.streamlit.app/, ensuring widespread availability to the research community.

# Availability

The complete software implementation, including installation instructions, dependencies, and detailed usage protocols, is freely available through the GitHub repository at https://github.com/yashmgupta/OpenPrimerP. The repository provides both access to the hosted web application and local installation options for institutional deployment. The open-source nature under the MIT License ensures reproducibility and enables community-driven improvements to support evolving needs in biodiversity genomics research.

# Acknowledgements

This research was supported by the Department of Biology, Faculty of Science, Naresuan University, Thailand (Grant No. R2565B018) and partially supported by the Global and Frontier Research University Fund, Naresuan University (Grant No. R2566C051). The authors gratefully acknowledge the graduate and undergraduate students of the Department of Biology, Faculty of Science, Naresuan University, for their invaluable contributions as early adopters and beta testers of `OpenPrimerP`. Their enthusiastic participation in the initial testing phases and constructive feedback were instrumental in identifying usability issues, refining the user interface, and enhancing the overall functionality of the application.

# References
