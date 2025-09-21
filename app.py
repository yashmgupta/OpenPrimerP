import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
import os
from io import StringIO


# Ensure the 'temp' directory exists for saving temporary files
temp_dir = "temp"
os.makedirs(temp_dir, exist_ok=True)

# Streamlit UI setup
st.set_page_config(page_title="PCR Primer Design", page_icon="🧬", layout="wide")

# User Documentation
st.sidebar.header("User Guide")
st.sidebar.info(
    """
    This application allows you to design PCR primers for specific features within a GenBank file. 
    Follow these steps:
    1. Upload a GenBank file.
    2. Select a feature type and then a specific feature.
    3. Enter the desired PCR product size range and the minimum number of primer pairs.
    4. Click 'Design Primers' to generate your primers.
    """
)

# File uploader with additional help
uploaded_file = st.file_uploader(
    "Upload a GenBank file", 
    type=['gb', 'gbk'],
    help="Upload a GenBank (.gb or .gbk) file containing the DNA sequence from which to design primers."
)

# Custom CSS for styling
st.markdown("""
<style>
h1, h2, h3 {
    color: #1f77b4;
    font-weight: bold;
}
.stProgress {
    margin-top: 20px;
}
</style>
""", unsafe_allow_html=True)

def extract_features_from_genbank(genbank_content, feature_types=['CDS', 'tRNA', 'gene']):
    """Extracts specified features from GenBank content."""
    text_stream = StringIO(genbank_content.decode("utf-8")) if isinstance(genbank_content, bytes) else genbank_content
    record = SeqIO.read(text_stream, "genbank")
    features = {ftype: [] for ftype in feature_types}
    for feature in record.features:
        if feature.type in feature_types:
            features[feature.type].append(feature)
    return features, record

def design_primers_for_region(sequence, product_size_range, num_to_return=10):
    """Design primers for a specific sequence."""
    size_min, size_max = map(int, product_size_range.split('-'))
    return primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': str(sequence),
            'PRIMER_PRODUCT_SIZE_RANGE': [[size_min, size_max]]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 23,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_NUM_RETURN': num_to_return,
        }
    )

if uploaded_file is not None:
    genbank_content = StringIO(uploaded_file.getvalue().decode("utf-8"))
    features, record = extract_features_from_genbank(genbank_content)
    
    st.write("## Feature Selection")
    feature_type = st.selectbox(
        'Select feature type:', 
        ['CDS', 'tRNA', 'gene'],
        help="Choose the type of genomic feature for which you want to design primers."
    )

    if features[feature_type]:
        feature_options = [f"{feature.qualifiers.get('gene', [''])[0]} ({feature.location})" for feature in features[feature_type]]
        selected_index = st.selectbox(
            f'Select a {feature_type}:', 
            options=range(len(feature_options)), 
            format_func=lambda x: feature_options[x],
            help="Select a specific feature based on its gene name and location."
        )
        selected_feature = features[feature_type][selected_index]

        feature_sequence = selected_feature.extract(record.seq)
        st.write(f"Selected {feature_type} sequence (length: {len(feature_sequence)} bp):")
        st.code(str(feature_sequence), language="text")  # Display sequence in code format

        st.write("## Primer Design Parameters")
        product_size_range = st.text_input(
            "Enter desired PCR product size range (e.g., 150-500):", 
            value="150-500",
            help="Specify the range of the desired PCR product size in base pairs (e.g., 150-500)."
        )
        min_num_primers = st.number_input(
            "Enter minimum number of primer pairs to return:", 
            min_value=5, value=5, step=1,
            help="Determine the minimum number of primer pairs to generate."
        )

        if st.button(f'Design Primers for selected {feature_type}'):
            with st.spinner('Designing primers...'):  # Show a spinner while primers are being designed
                primers = design_primers_for_region(feature_sequence, product_size_range, num_to_return=min_num_primers)

            primer_data = []
            for i in range(min_num_primers):
                left_sequence = primers.get(f'PRIMER_LEFT_{i}_SEQUENCE', 'N/A')
                right_sequence = primers.get(f'PRIMER_RIGHT_{i}_SEQUENCE', 'N/A')
                if left_sequence != 'N/A' and right_sequence != 'N/A':
                    primer_info = {
                        'Primer Pair': i + 1,
                        'Left Sequence': left_sequence,
                        'Right Sequence': right_sequence,
                        'Left TM (°C)': primers.get(f'PRIMER_LEFT_{i}_TM', 'N/A'),
                        'Right TM (°C)': primers.get(f'PRIMER_RIGHT_{i}_TM', 'N/A'),
                        'Left Length': len(left_sequence),
                        'Right Length': len(right_sequence),
                        'PCR Product Size (bp)': primers.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 'N/A')
                    }
                    primer_data.append(primer_info)

            if primer_data:
                st.subheader('Designed Primers')
                primer_df = pd.DataFrame(primer_data)
                st.table(primer_df)  # Use st.table to display the primer data

                csv = primer_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    "Download Primers as CSV",
                    csv,
                    "primers.csv",
                    "text/csv",
                    key='download-csv'
                )
            else:
                st.error('No primers were found. Please adjust your parameters and try again.')


# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2024 Yash Munnalal Gupta. All rights reserved.

For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
