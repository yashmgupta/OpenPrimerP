import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
import os
from io import StringIO
import re
from Bio.Seq import Seq

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
    1. Upload a GenBank file (can contain multiple records).
    2. Select a record from the file.
    3. Select a feature type and then a specific feature.
    4. Enter the desired PCR product size range and the minimum number of primer pairs.
    5. Click 'Design Primers' to generate your primers.
    6. Check primer specificity against all sequences in your file.
    """
)

# File uploader with additional help
uploaded_file = st.file_uploader(
    "Upload a GenBank file", 
    type=['gb', 'gbk'],
    help="Upload a GenBank (.gb or .gbk) file containing the DNA sequence(s) from which to design primers. Multiple records are supported."
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
    background-color: yellow;
    padding: 2px;
}
.primer-match {
    background-color: lightgreen;
    padding: 2px;
}
.primer-mismatch {
    background-color: lightcoral;
    padding: 2px;
}
</style>
""", unsafe_allow_html=True)

def extract_records_from_genbank(genbank_content):
    """Extracts all records from GenBank content."""
    text_stream = StringIO(genbank_content.decode("utf-8")) if isinstance(genbank_content, bytes) else genbank_content
    records = list(SeqIO.parse(text_stream, "genbank"))
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

def check_primer_specificity(forward_primer, reverse_primer, records, max_product_size=5000):
    """
    Check if the primer pair can potentially amplify regions in other sequences.
    Returns a list of potential amplification products.
    """
    results = []
    
    # Convert primers to uppercase for case-insensitive matching
    forward_primer_upper = forward_primer.upper()
    reverse_primer_upper = reverse_primer.upper()
    reverse_complement = str(Seq(reverse_primer).reverse_complement())
    reverse_complement_upper = reverse_complement.upper()
    
    for i, record in enumerate(records):
        seq_str = str(record.seq).upper()
        
        # Find all occurrences of forward primer
        forward_matches = [m.start() for m in re.finditer(forward_primer_upper, seq_str)]
        
        # Find all occurrences of reverse primer (reverse complement)
        reverse_matches = [m.start() for m in re.finditer(reverse_complement_upper, seq_str)]
        
        # Check potential amplification products
        for f_pos in forward_matches:
            for r_pos in reverse_matches:
                # Check if reverse primer is downstream of forward primer
                if r_pos > f_pos:
                    product_size = r_pos - f_pos + len(reverse_primer)
                    if product_size <= max_product_size:
                        results.append({
                            "record_id": record.id,
                            "record_index": i,
                            "forward_position": f_pos,
                            "reverse_position": r_pos,
                            "product_size": product_size
                        })
    
    return results

def visualize_primer_binding(sequence, forward_primer, reverse_primer, forward_pos, reverse_pos):
    """
    Create a visualization of where primers bind to the sequence.
    Returns HTML with highlighted regions.
    """
    sequence_str = str(sequence)
    result = ""
    
    # Add sequence before forward primer
    result += sequence_str[:forward_pos]
    
    # Add forward primer (highlighted)
    result += f'<span class="primer-match">{sequence_str[forward_pos:forward_pos+len(forward_primer)]}</span>'
    
    # Add sequence between primers
    middle_start = forward_pos + len(forward_primer)
    middle_end = reverse_pos
    result += sequence_str[middle_start:middle_end]
    
    # Add reverse primer (highlighted)
    result += f'<span class="primer-match">{sequence_str[reverse_pos:reverse_pos+len(reverse_primer)]}</span>'
    
    # Add sequence after reverse primer
    result += sequence_str[reverse_pos+len(reverse_primer):]
    
    return result

if uploaded_file is not None:
    # Parse all records from the GenBank file
    records = extract_records_from_genbank(uploaded_file)
    
    if not records:
        st.error("No valid records found in the GenBank file.")
    else:
        st.write(f"## Found {len(records)} records in the GenBank file")
        
        # Let user select which record to work with
        record_options = [f"{record.id} ({len(record.seq)} bp)" for record in records]
        selected_record_index = st.selectbox(
            'Select a record to design primers for:', 
            options=range(len(record_options)), 
            format_func=lambda x: record_options[x],
            help="Choose which sequence record to use for primer design."
        )
        
        selected_record = records[selected_record_index]
        st.write(f"Working with record: **{selected_record.id}** ({len(selected_record.seq)} bp)")
        
        # Extract features from the selected record
        features = extract_features_from_record(selected_record)
        
        st.write("## Feature Selection")
        feature_type = st.selectbox(
            'Select feature type:', 
            ['CDS', 'tRNA', 'gene'],
            help="Choose the type of genomic feature for which you want to design primers."
        )

        if features[feature_type]:
            feature_options = [f"{feature.qualifiers.get('gene', [''])[0] or feature.qualifiers.get('locus_tag', ['Unknown'])[0]} ({feature.location})" for feature in features[feature_type]]
            selected_index = st.selectbox(
                f'Select a {feature_type}:', 
                options=range(len(feature_options)), 
                format_func=lambda x: feature_options[x],
                help="Select a specific feature based on its gene name and location."
            )
            selected_feature = features[feature_type][selected_index]

            feature_sequence = selected_feature.extract(selected_record.seq)
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
                    st.dataframe(primer_df)  # Use st.dataframe instead of st.table for better interactivity
                    
                    # Allow user to select a primer pair for specificity check
                    st.subheader("Primer Specificity Check")
                    st.info("Select a primer pair to check its specificity against all sequences in your file.")
                    
                    selected_primer_idx = st.selectbox(
                        "Select primer pair to check:", 
                        options=range(len(primer_data)),
                        format_func=lambda x: f"Primer Pair {primer_data[x]['Primer Pair']}"
                    )
                    
                    if st.button("Check Primer Specificity"):
                        selected_primer = primer_data[selected_primer_idx]
                        forward_primer = selected_primer['Left Sequence']
                        reverse_primer = selected_primer['Right Sequence']
                        
                        with st.spinner("Analyzing primer specificity..."):
                            specificity_results = check_primer_specificity(
                                forward_primer, reverse_primer, records
                            )
                        
                        st.write(f"### Specificity Results for Primer Pair {selected_primer['Primer Pair']}")
                        st.write(f"Forward primer: `{forward_primer}`")
                        st.write(f"Reverse primer: `{reverse_primer}`")
                        
                        if specificity_results:
                            st.write(f"Found {len(specificity_results)} potential amplification products:")
                            
                            # Create a DataFrame for all potential products
                            specificity_df = pd.DataFrame(specificity_results)
                            st.dataframe(specificity_df)
                            
                            # Allow user to visualize a specific binding site
                            if st.checkbox("Show primer binding visualization"):
                                binding_idx = st.selectbox(
                                    "Select binding site to visualize:",
                                    options=range(len(specificity_results)),
                                    format_func=lambda x: f"Record {specificity_results[x]['record_id']} - Product size: {specificity_results[x]['product_size']} bp"
                                )
                                
                                binding_site = specificity_results[binding_idx]
                                record_idx = binding_site["record_index"]
                                forward_pos = binding_site["forward_position"]
                                reverse_pos = binding_site["reverse_position"]
                                
                                # Extract sequence region with some context
                                context_size = 50  # Amount of sequence to show before and after
                                start_pos = max(0, forward_pos - context_size)
                                end_pos = min(len(records[record_idx].seq), reverse_pos + len(reverse_primer) + context_size)
                                
                                sequence_region = records[record_idx].seq[start_pos:end_pos]
                                
                                # Adjust positions for the extracted region
                                adj_forward_pos = forward_pos - start_pos
                                adj_reverse_pos = reverse_pos - start_pos
                                
                                # Create visualization
                                html_vis = visualize_primer_binding(
                                    sequence_region, 
                                    forward_primer, 
                                    reverse_primer, 
                                    adj_forward_pos, 
                                    adj_reverse_pos
                                )
                                
                                st.write("Primer binding visualization (primers highlighted in green):")
                                st.markdown(f"<div style='font-family:monospace; overflow-wrap:break-word;'>{html_vis}</div>", unsafe_allow_html=True)
                        else:
                            st.success("No non-specific amplification products found. These primers appear to be specific to your target region.")
                    
                    # Download option
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
        else:
            st.warning(f"No {feature_type} features found in the selected record.")


# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2024 Yash Munnalal Gupta. All rights reserved.

For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
