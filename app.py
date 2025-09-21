import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import re
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
    1. Upload a GenBank file (single or multiple records).
    2. Select a record (if multiple are present).
    3. Select a feature type and then a specific feature.
    4. Enter the desired PCR product size range and the minimum number of primer pairs.
    5. Click 'Design Primers' to generate your primers.
    6. For multi-record files, you can also check primer specificity against other sequences.
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
.highlight {
    background-color: #ffff99;
    padding: 2px;
}
</style>
""", unsafe_allow_html=True)

def extract_records_from_genbank(genbank_file):
    """Extract all records from a GenBank file."""
    # Reset the file pointer to the beginning
    genbank_file.seek(0)
    # Read the file content as text
    content = genbank_file.read().decode("utf-8")
    # Create a StringIO object
    text_stream = StringIO(content)
    # Parse all records
    records = list(SeqIO.parse(text_stream, "genbank"))
    return records

def extract_features_from_record(record, feature_types=['CDS', 'tRNA', 'gene']):
    """Extracts specified features from a GenBank record."""
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

def simulate_pcr(forward_primer, reverse_primer, sequence, max_product_size=3000):
    """
    Simulate PCR to find potential amplification products.
    Fixed logic: Forward primer binds as-is, reverse primer binds as reverse complement.
    """
    results = []
    
    # Convert to string if sequence is a Seq object
    if isinstance(sequence, Seq):
        sequence = str(sequence).upper()
    else:
        sequence = str(sequence).upper()
    
    # Convert primers to uppercase for consistent matching
    forward_primer = forward_primer.upper()
    reverse_primer = reverse_primer.upper()
    
    # Find all forward primer binding sites (exact match)
    forward_matches = []
    for i in range(len(sequence) - len(forward_primer) + 1):
        if sequence[i:i+len(forward_primer)] == forward_primer:
            forward_matches.append(i)
    
    # Find all reverse primer binding sites (reverse complement)
    reverse_primer_rc = str(Seq(reverse_primer).reverse_complement()).upper()
    reverse_matches = []
    for i in range(len(sequence) - len(reverse_primer_rc) + 1):
        if sequence[i:i+len(reverse_primer_rc)] == reverse_primer_rc:
            reverse_matches.append(i)
    
    # Check all possible combinations for valid PCR products
    for f_pos in forward_matches:
        for r_pos in reverse_matches:
            # Check if reverse primer is downstream and within max size
            if r_pos > f_pos and r_pos - f_pos <= max_product_size:
                product_size = r_pos - f_pos + len(reverse_primer_rc)
                results.append({
                    'start': f_pos,
                    'end': r_pos + len(reverse_primer_rc),
                    'size': product_size
                })
    
    return results

def check_primer_specificity(forward_primer, reverse_primer, records, target_record_id, target_feature=None, max_product_size=3000):
    """
    Check primer specificity across all records.
    Fixed logic: Always check target feature first, then check all other features.
    """
    specificity_results = []
    
    for record in records:
        record_id = record.id
        is_target_record = (record_id == target_record_id)
        
        # First, if this is the target record and we have a target feature, check it specifically
        if is_target_record and target_feature is not None:
            try:
                # Extract target feature sequence
                region_seq = target_feature.extract(record.seq)
                pcr_products = simulate_pcr(forward_primer, reverse_primer, region_seq, max_product_size)
                
                # Add target results
                for product in pcr_products:
                    result = {
                        'Record': f"{record_id} (target)",
                        'Location': f"{target_feature.location} (target feature)",
                        'Product Size': product['size'],
                        'Is Target': True
                    }
                    specificity_results.append(result)
            except Exception as e:
                st.warning(f"Could not extract target feature sequence: {e}")
        
        # Then check all features in all records (including target record for off-target effects)
        features = extract_features_from_record(record)
        for feature_type in features:
            for feature in features[feature_type]:
                # Skip the exact target feature to avoid duplication
                if is_target_record and target_feature is not None and feature == target_feature:
                    continue
                
                try:
                    # Get feature name if available
                    feature_name = feature.qualifiers.get('gene', ['Unknown'])[0]
                    # Extract feature sequence
                    region_seq = feature.extract(record.seq)
                    # Simulate PCR
                    pcr_products = simulate_pcr(forward_primer, reverse_primer, region_seq, max_product_size)
                    
                    for product in pcr_products:
                        result = {
                            'Record': record_id,
                            'Location': f"{feature_type}:{feature_name} ({feature.location})",
                            'Product Size': product['size'],
                            'Is Target': False  # These are all non-target results
                        }
                        specificity_results.append(result)
                except Exception as e:
                    # Skip features that can't be processed
                    continue
    
    return specificity_results

if uploaded_file is not None:
    try:
        # Extract all records from the GenBank file
        records = extract_records_from_genbank(uploaded_file)
        
        if not records:
            st.error("No valid GenBank records found in the file.")
        else:
            # Display record count
            st.write(f"## File contains {len(records)} record(s)")
            
            # If multiple records, let user select which one to use for primer design
            if len(records) > 1:
                record_options = [f"{record.id}: {record.description[:50]}..." for record in records]
                selected_record_index = st.selectbox(
                    'Select a record for primer design:',
                    options=range(len(record_options)),
                    format_func=lambda x: record_options[x]
                )
                record = records[selected_record_index]
                st.write(f"Selected record: **{record.id}**")
            else:
                record = records[0]
                st.write(f"Record: **{record.id}**")
            
            # Extract features from the selected record
            features = extract_features_from_record(record)
            
            st.write("## Feature Selection")
            feature_type = st.selectbox(
                'Select feature type:', 
                ['CDS', 'tRNA', 'gene'],
                help="Choose the type of genomic feature for which you want to design primers."
            )
            
            if features[feature_type]:
                feature_options = []
                for feature in features[feature_type]:
                    # Get the gene name if available, otherwise use a placeholder
                    gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
                    feature_options.append(f"{gene_name} ({feature.location})")
                
                selected_index = st.selectbox(
                    f'Select a {feature_type}:', 
                    options=range(len(feature_options)), 
                    format_func=lambda x: feature_options[x],
                    help="Select a specific feature based on its gene name and location."
                )
                selected_feature = features[feature_type][selected_index]
                feature_sequence = selected_feature.extract(record.seq)
                
                st.write(f"Selected {feature_type} sequence (length: {len(feature_sequence)} bp):")
                st.code(str(feature_sequence), language="text")
                
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
                    with st.spinner('Designing primers...'):
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
                        st.session_state['selected_record'] = record.id
                        st.session_state['records'] = records
                        st.session_state['selected_feature'] = selected_feature
                        
                        st.subheader('Designed Primers')
                        primer_df = pd.DataFrame(primer_data)
                        st.table(primer_df)
                        
                        csv = primer_df.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            "Download Primers as CSV",
                            csv,
                            "primers.csv",
                            "text/csv",
                            key='download-csv'
                        )
                        
                        # Primer specificity check option (available for single or multiple records)
                        if 'primer_data' in st.session_state:
                            st.write("## Primer Specificity Check")
                            if len(records) > 1:
                                st.write("Check if your primers could amplify unintended regions in other sequences.")
                            else:
                                st.write("Check primer specificity within the current sequence.")
                            
                            primer_options = [f"Pair {i+1}: {p['Left Sequence']} / {p['Right Sequence']}" 
                                             for i, p in enumerate(st.session_state['primer_data'])]
                            
                            selected_pair = st.selectbox(
                                "Select a primer pair to check:", 
                                options=range(len(primer_options)),
                                format_func=lambda x: primer_options[x]
                            )
                            
                            if st.button("Check Primer Specificity"):
                                selected_primer = st.session_state['primer_data'][selected_pair]
                                forward_primer = selected_primer['Left Sequence']
                                reverse_primer = selected_primer['Right Sequence']
                                
                                with st.spinner("Checking primer specificity..."):
                                    specificity_results = check_primer_specificity(
                                        forward_primer,
                                        reverse_primer,
                                        st.session_state['records'],
                                        st.session_state['selected_record'],
                                        st.session_state['selected_feature']
                                    )
                                
                                if specificity_results:
                                    st.write(f"Found {len(specificity_results)} potential amplification products:")
                                    
                                    # Convert to DataFrame for better display
                                    spec_df = pd.DataFrame(specificity_results)
                                    
                                    # Style the target row
                                    def highlight_target(row):
                                        if row['Is Target']:
                                            return ['background-color: #c6efce'] * len(row)  # Green for target
                                        else:
                                            return ['background-color: #ffc7ce'] * len(row)  # Red for off-target
                                    
                                    # Display styled dataframe
                                    st.dataframe(spec_df.style.apply(highlight_target, axis=1))
                                    
                                    # Count target and non-target amplifications
                                    target_count = sum(1 for r in specificity_results if r['Is Target'])
                                    non_target_count = sum(1 for r in specificity_results if not r['Is Target'])
                                    
                                    st.write(f"**Results Summary:**")
                                    st.write(f"- Target amplifications: {target_count}")
                                    st.write(f"- Non-specific amplifications: {non_target_count}")
                                    
                                    if non_target_count > 0:
                                        st.warning(f"⚠️ Found {non_target_count} potential non-specific amplification products!")
                                    else:
                                        st.success("✅ No non-specific amplification products found!")
                                    
                                    if target_count == 0:
                                        st.error("⚠️ No target amplification found! This suggests the primers may not work as expected.")
                                else:
                                    st.warning("⚠️ No amplification products found. This could indicate:")
                                    st.write("- Primers may not bind to any sequences")
                                    st.write("- PCR product size may be outside the search range")
                                    st.write("- Sequences may be too short for the primers to bind properly")
                    else:
                        st.error('No primers were found. Please adjust your parameters and try again.')
            else:
                st.warning(f"No {feature_type} features found in the selected record.")
    
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Please check that your GenBank file is valid and try again.")
        # Add debug information
        import traceback
        st.error("Debug information:")
        st.code(traceback.format_exc())

# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2024 Yash Munnalal Gupta. All rights reserved.
For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
