import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
import os
import re
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
    1. Upload a GenBank file (single or multiple records).
    2. Select a record, then a feature type and specific feature.
    3. Enter the desired PCR product size range and the minimum number of primer pairs.
    4. Click 'Design Primers' to generate your primers.
    5. Use the 'Primer Specificity Check' to see if primers might amplify other regions.
    """
)

# File uploader with additional help
uploaded_file = st.file_uploader(
    "Upload a GenBank file", 
    type=['gb', 'gbk'],
    help="Upload a GenBank (.gb or .gbk) file containing one or more DNA sequences."
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

def extract_records_from_genbank(file_content):
    """Extract all records from a GenBank file."""
    # Reset the file pointer if it exists
    if hasattr(file_content, 'seek'):
        file_content.seek(0)
    
    # Convert to StringIO if bytes or file object
    if isinstance(file_content, bytes):
        text_stream = StringIO(file_content.decode("utf-8"))
    elif hasattr(file_content, 'read'):
        text_stream = StringIO(file_content.getvalue().decode("utf-8"))
    else:
        text_stream = file_content
        
    # Parse all records
    records = list(SeqIO.parse(text_stream, "genbank"))
    if not records:
        st.error("No valid GenBank records found in the file.")
    
    return records

def extract_features_from_record(record, feature_types=['CDS', 'tRNA', 'gene']):
    """Extracts specified features from a single GenBank record."""
    features = {ftype: [] for ftype in feature_types}
    for feature in record.features:
        if feature.type in feature_types:
            features[feature.type].append(feature)
    return features

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

def simulate_pcr(forward_primer, reverse_primer, sequence, min_size=50, max_size=5000):
    """
    Simple in silico PCR simulation to find potential amplification sites.
    Returns a list of (start, end, size) for potential amplification products.
    """
    # Convert primers to uppercase
    forward_primer = forward_primer.upper()
    reverse_primer = reverse_primer.upper()
    sequence_str = str(sequence).upper()
    
    # Find all occurrences of forward primer (allowing for some mismatches - simplified)
    # In a real implementation, you might want to use a more sophisticated algorithm
    # that allows for a specific number of mismatches
    forward_sites = [m.start() for m in re.finditer(forward_primer, sequence_str, re.IGNORECASE)]
    
    # Find reverse complement of reverse primer
    reverse_complement = str(reverse_primer).translate(str.maketrans('ATCG', 'TAGC'))[::-1]
    reverse_sites = [m.start() for m in re.finditer(reverse_complement, sequence_str, re.IGNORECASE)]
    
    # Find potential amplification products
    products = []
    for f_site in forward_sites:
        for r_site in reverse_sites:
            if f_site < r_site:  # Forward primer comes before reverse primer
                product_size = r_site + len(reverse_primer) - f_site
                if min_size <= product_size <= max_size:
                    products.append((f_site, r_site + len(reverse_primer), product_size))
    
    return products

if uploaded_file is not None:
    try:
        # Extract all records from the GenBank file
        records = extract_records_from_genbank(uploaded_file)
        
        if records:
            # Record selection if there are multiple records
            if len(records) > 1:
                record_names = [f"{record.id} - {record.description[:50]}..." for record in records]
                selected_record_idx = st.selectbox(
                    "Select a record:",
                    range(len(records)),
                    format_func=lambda x: record_names[x],
                    help="Choose which record to design primers for."
                )
                current_record = records[selected_record_idx]
                st.success(f"Selected record: {current_record.id} ({len(current_record.seq)} bp)")
            else:
                current_record = records[0]
                st.success(f"Loaded record: {current_record.id} ({len(current_record.seq)} bp)")
            
            # Extract features for the selected record
            features = extract_features_from_record(current_record)
            
            st.write("## Feature Selection")
            # Filter to only show feature types that exist in this record
            available_feature_types = [ft for ft in ['CDS', 'tRNA', 'gene'] if features[ft]]
            
            if available_feature_types:
                feature_type = st.selectbox(
                    'Select feature type:', 
                    available_feature_types,
                    help="Choose the type of genomic feature for which you want to design primers."
                )
                
                if features[feature_type]:
                    feature_options = []
                    for feature in features[feature_type]:
                        # Try to get gene name or locus_tag or product name
                        name = (
                            feature.qualifiers.get('gene', [''])[0] or 
                            feature.qualifiers.get('locus_tag', [''])[0] or
                            feature.qualifiers.get('product', [''])[0] or
                            'Unknown'
                        )
                        feature_options.append(f"{name} ({feature.location})")
                    
                    selected_index = st.selectbox(
                        f'Select a {feature_type}:', 
                        options=range(len(feature_options)), 
                        format_func=lambda x: feature_options[x],
                        help="Select a specific feature based on its name and location."
                    )
                    selected_feature = features[feature_type][selected_index]
                    
                    try:
                        feature_sequence = selected_feature.extract(current_record.seq)
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
                            min_value=1, max_value=20, value=5, step=1,
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
                                st.session_state['primer_data'] = primer_data
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
                                
                                # Add primer specificity check section
                                st.markdown("---")
                                st.header("Primer Specificity Check")
                                st.write("Check if your primers could amplify unintended regions in other sequences.")
                                
                                # Create options for primer pair selection
                                primer_options = [f"Pair {data['Primer Pair']}: {data['Left Sequence']} / {data['Right Sequence']}" 
                                                for data in primer_data]
                                
                                selected_pair_idx = st.selectbox(
                                    "Select a primer pair to check:",
                                    range(len(primer_options)),
                                    format_func=lambda x: primer_options[x]
                                )
                                
                                # Get the selected primer pair
                                selected_pair = primer_data[selected_pair_idx]
                                forward_primer = selected_pair['Left Sequence']
                                reverse_primer = selected_pair['Right Sequence']
                                
                                if st.button("Check Primer Specificity"):
                                    with st.spinner("Analyzing primer specificity..."):
                                        # Check primers against all records in the file
                                        specificity_results = []
                                        
                                        for idx, record in enumerate(records):
                                            record_id = record.id
                                            record_desc = record.description[:50]
                                            potential_products = simulate_pcr(
                                                forward_primer, 
                                                reverse_primer, 
                                                record.seq,
                                                min_size=50,
                                                max_size=5000
                                            )
                                            
                                            # Add all potential products to results
                                            for start, end, size in potential_products:
                                                specificity_results.append({
                                                    'Record': f"{record_id} - {record_desc}",
                                                    'Start': start + 1,  # 1-based indexing for display
                                                    'End': end,
                                                    'Product Size (bp)': size,
                                                    'Is Target': idx == selected_record_idx and \
                                                                start >= selected_feature.location.start and \
                                                                end <= selected_feature.location.end
                                                })
                                        
                                        # Display results
                                        if specificity_results:
                                            st.success(f"Found {len(specificity_results)} potential amplification products")
                                            
                                            # Convert results to DataFrame and display
                                            results_df = pd.DataFrame(specificity_results)
                                            # Highlight the target amplification
                                            st.dataframe(
                                                results_df.style.apply(
                                                    lambda row: ['background-color: lightgreen' if row['Is Target'] else '' 
                                                                for _ in row],
                                                    axis=1
                                                )
                                            )
                                            
                                            # Create downloadable CSV
                                            csv = results_df.to_csv(index=False).encode('utf-8')
                                            st.download_button(
                                                "Download Specificity Results as CSV",
                                                csv,
                                                "primer_specificity_results.csv",
                                                "text/csv",
                                                key='download-specificity-csv'
                                            )
                                        else:
                                            st.warning("No potential amplification products found for these primers.")
                            else:
                                st.error('No primers were found. Please adjust your parameters and try again.')
                    except Exception as e:
                        st.error(f"Error processing feature sequence: {str(e)}")
                else:
                    st.warning(f"No {feature_type} features found in this record.")
            else:
                st.warning("No recognizable features (CDS, tRNA, gene) found in this record.")
    except Exception as e:
        st.error(f"Error processing GenBank file: {str(e)}")

# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2024 Yash Munnalal Gupta. All rights reserved.

For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
