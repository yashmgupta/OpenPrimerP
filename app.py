import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import os
from io import StringIO
import re

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
    
    ### Basic Primer Design:
    1. Upload a GenBank file (single or multiple records)
    2. Select a record (if multiple are present)
    3. Select a feature type and then a specific feature
    4. Enter the desired PCR product size range and the minimum number of primer pairs
    5. Click 'Design Primers' to generate your primers
    
    ### Primer Specificity Check:
    1. Design primers or select existing ones
    2. Run specificity check against all uploaded sequences
    3. View potential off-target amplifications
    """
)

# File uploader with additional help
uploaded_file = st.file_uploader(
    "Upload a GenBank file", 
    type=['gb', 'gbk'],
    help="Upload a GenBank (.gb or .gbk) file containing the DNA sequence(s) from which to design primers."
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
.specificity-hit {
    background-color: #ffcccc;
}
.specificity-clean {
    background-color: #ccffcc;
}
.primer-table {
    margin-bottom: 20px;
}
</style>
""", unsafe_allow_html=True)

# Session state initialization
if 'records' not in st.session_state:
    st.session_state.records = []
if 'designed_primers' not in st.session_state:
    st.session_state.designed_primers = []

def extract_features_from_genbank(genbank_content, feature_types=['CDS', 'tRNA', 'gene']):
    """Extracts specified features from GenBank content, handling multiple records."""
    text_stream = StringIO(genbank_content.decode("utf-8")) if isinstance(genbank_content, bytes) else genbank_content
    records = list(SeqIO.parse(text_stream, "genbank"))
    
    if not records:
        st.error("No valid GenBank records found in the uploaded file.")
        return None, None
    
    st.session_state.records = records
    return records

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

def check_primer_binding(primer_seq, template_seq, max_mismatches=3):
    """Check if a primer can bind to a template with up to max_mismatches."""
    primer_len = len(primer_seq)
    template_len = len(template_seq)
    
    # Convert to strings if they're not already
    if isinstance(primer_seq, Seq):
        primer_seq = str(primer_seq)
    if isinstance(template_seq, Seq):
        template_seq = str(template_seq)
    
    # Check forward strand
    binding_sites = []
    for i in range(template_len - primer_len + 1):
        template_segment = template_seq[i:i+primer_len]
        mismatches = sum(a != b for a, b in zip(primer_seq, template_segment))
        if mismatches <= max_mismatches:
            binding_sites.append((i, i+primer_len-1, mismatches, '+'))
    
    # Check reverse strand
    rev_comp_template = str(Seq(template_seq).reverse_complement())
    for i in range(template_len - primer_len + 1):
        template_segment = rev_comp_template[i:i+primer_len]
        mismatches = sum(a != b for a, b in zip(primer_seq, template_segment))
        if mismatches <= max_mismatches:
            binding_sites.append((template_len-i-primer_len, template_len-i-1, mismatches, '-'))
    
    return binding_sites

def simulate_pcr(forward_primer, reverse_primer, template_seq, max_product_size=5000, max_mismatches=3):
    """Simulate PCR to predict amplification products."""
    forward_sites = check_primer_binding(forward_primer, template_seq, max_mismatches)
    reverse_sites = check_primer_binding(reverse_primer, template_seq, max_mismatches)
    
    products = []
    
    for f_start, f_end, f_mm, f_strand in forward_sites:
        for r_start, r_end, r_mm, r_strand in reverse_sites:
            # Ensure correct orientation (forward primer binds to sense strand, reverse primer to antisense)
            if f_strand == '+' and r_strand == '-':
                if f_start < r_end:  # Correct orientation
                    product_size = r_end - f_start + 1
                    if product_size <= max_product_size:
                        products.append({
                            'f_start': f_start,
                            'f_end': f_end,
                            'r_start': r_start,
                            'r_end': r_end,
                            'f_mismatches': f_mm,
                            'r_mismatches': r_mm,
                            'product_size': product_size,
                            'total_mismatches': f_mm + r_mm
                        })
                        
    return sorted(products, key=lambda x: x['product_size'])

# Create tabs for different functionalities
tab1, tab2 = st.tabs(["Primer Design", "Specificity Check"])

with tab1:
    st.header("Primer Design")
    
    if uploaded_file is not None:
        records = extract_features_from_genbank(uploaded_file.getvalue())
        
        if records and len(records) > 0:
            # If multiple records, let user select which one to work with
            if len(records) > 1:
                record_names = [rec.id for rec in records]
                selected_record_idx = st.selectbox(
                    "Select a record to design primers for:",
                    range(len(record_names)),
                    format_func=lambda i: f"{record_names[i]} ({len(records[i].seq)} bp)"
                )
                record = records[selected_record_idx]
                st.info(f"Working with record: {record.id} ({len(record.seq)} bp)")
            else:
                record = records[0]
                st.info(f"Loaded record: {record.id} ({len(record.seq)} bp)")
            
            # Extract features for the selected record
            feature_types = ['CDS', 'tRNA', 'gene', 'misc_feature']
            features = {ftype: [] for ftype in feature_types}
            
            for feature in record.features:
                if feature.type in feature_types:
                    features[feature.type].append(feature)
            
            st.write("## Feature Selection")
            available_types = [ft for ft in feature_types if features[ft]]
            
            if not available_types:
                st.warning("No common features (CDS, tRNA, gene) found in this record.")
            else:
                feature_type = st.selectbox(
                    'Select feature type:', 
                    available_types,
                    help="Choose the type of genomic feature for which you want to design primers."
                )
    
                if features[feature_type]:
                    # Create more descriptive labels for features
                    feature_options = []
                    for feature in features[feature_type]:
                        label = ""
                        if 'gene' in feature.qualifiers:
                            label += f"Gene: {feature.qualifiers['gene'][0]}, "
                        elif 'product' in feature.qualifiers:
                            label += f"Product: {feature.qualifiers['product'][0]}, "
                        elif 'note' in feature.qualifiers:
                            note = feature.qualifiers['note'][0]
                            label += f"Note: {note[:20]}..., " if len(note) > 20 else f"Note: {note}, "
                        label += f"Location: {feature.location}"
                        feature_options.append(label)
                    
                    selected_index = st.selectbox(
                        f'Select a {feature_type}:', 
                        options=range(len(feature_options)), 
                        format_func=lambda x: feature_options[x],
                        help="Select a specific feature based on its gene name and location."
                    )
                    selected_feature = features[feature_type][selected_index]
    
                    feature_sequence = selected_feature.extract(record.seq)
                    st.write(f"Selected {feature_type} sequence (length: {len(feature_sequence)} bp):")
                    
                    # Show a snippet of the sequence if it's very long
                    if len(feature_sequence) > 300:
                        st.code(f"{str(feature_sequence[:100])}...{str(feature_sequence[-100:])}", language="text")
                        st.info(f"Showing first and last 100 bp of a {len(feature_sequence)} bp sequence")
                    else:
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
                                    'Left TM (°C)': round(primers.get(f'PRIMER_LEFT_{i}_TM', 0), 2),
                                    'Right TM (°C)': round(primers.get(f'PRIMER_RIGHT_{i}_TM', 0), 2),
                                    'Left Length': len(left_sequence),
                                    'Right Length': len(right_sequence),
                                    'PCR Product Size (bp)': primers.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 'N/A')
                                }
                                primer_data.append(primer_info)
    
                        if primer_data:
                            st.session_state.designed_primers = primer_data
                            
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
                            
                            st.success("Primers designed successfully! Go to the 'Specificity Check' tab to analyze primer specificity.")
                        else:
                            st.error('No primers were found. Please adjust your parameters and try again.')

with tab2:
    st.header("Primer Specificity Check")
    
    if not st.session_state.records:
        st.warning("Please upload a GenBank file and design primers first.")
    elif not st.session_state.designed_primers:
        st.warning("Please design primers in the 'Primer Design' tab first.")
    else:
        st.write("## Select Primers for Specificity Analysis")
        
        # Let user select which primer pair to analyze
        primer_pairs = [
            f"Pair {p['Primer Pair']}: {p['Left Sequence']} + {p['Right Sequence']} (Product: {p['PCR Product Size (bp)']} bp)" 
            for p in st.session_state.designed_primers
        ]
        
        selected_pair_idx = st.selectbox(
            "Select primer pair to analyze:",
            range(len(primer_pairs)),
            format_func=lambda i: primer_pairs[i]
        )
        
        selected_pair = st.session_state.designed_primers[selected_pair_idx]
        forward_primer = selected_pair['Left Sequence']
        reverse_primer = selected_pair['Right Sequence']
        
        st.write("## Specificity Parameters")
        max_product_size = st.number_input("Maximum product size to consider (bp):", min_value=100, value=5000, step=100)
        max_mismatches = st.number_input("Maximum allowed mismatches per primer:", min_value=0, value=3, step=1)
        
        if st.button("Analyze Primer Specificity"):
            with st.spinner("Analyzing primer specificity across all sequences..."):
                results = []
                
                for idx, record in enumerate(st.session_state.records):
                    template_seq = str(record.seq)
                    
                    # Simulate PCR
                    products = simulate_pcr(
                        forward_primer, reverse_primer, template_seq, 
                        max_product_size, max_mismatches
                    )
                    
                    if products:
                        for product in products:
                            results.append({
                                'Record': record.id,
                                'Record Index': idx,
                                'Forward Start': product['f_start'] + 1,  # Convert to 1-based indexing
                                'Reverse End': product['r_end'] + 1,      # Convert to 1-based indexing
                                'Product Size': product['product_size'],
                                'Forward Mismatches': product['f_mismatches'],
                                'Reverse Mismatches': product['r_mismatches'],
                                'Total Mismatches': product['total_mismatches']
                            })
                
                if results:
                    st.subheader("Potential Amplification Products")
                    
                    # Create DataFrame and sort by mismatches and product size
                    results_df = pd.DataFrame(results).sort_values(
                        by=['Total Mismatches', 'Product Size']
                    )
                    
                    # Highlight the perfect matches (0 mismatches)
                    def highlight_matches(val):
                        if val == 0:
                            return 'background-color: #ccffcc'  # Light green
                        elif val <= max_mismatches:
                            return 'background-color: #ffffcc'  # Light yellow
                        else:
                            return 'background-color: #ffcccc'  # Light red
                    
                    # Apply styling
                    styled_df = results_df.style.applymap(
                        highlight_matches, 
                        subset=['Forward Mismatches', 'Reverse Mismatches', 'Total Mismatches']
                    )
                    
                    st.dataframe(styled_df)
                    
                    # Generate summary
                    st.subheader("Specificity Summary")
                    perfect_matches = results_df[results_df['Total Mismatches'] == 0].shape[0]
                    partial_matches = results_df[(results_df['Total Mismatches'] > 0) & 
                                                (results_df['Total Mismatches'] <= max_mismatches)].shape[0]
                    
                    st.write(f"**Perfect matches (0 mismatches):** {perfect_matches}")
                    st.write(f"**Partial matches (1-{max_mismatches} mismatches):** {partial_matches}")
                    
                    if perfect_matches > 1:
                        st.warning(f"⚠️ These primers may amplify {perfect_matches} different regions across all sequences.")
                    elif perfect_matches == 1:
                        st.success("✅ These primers are highly specific (only one perfect match found).")
                    else:
                        st.error("❌ No perfect matches found. Check your primer sequences.")
                        
                    # Offer download option
                    csv = results_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        "Download Specificity Results as CSV",
                        csv,
                        "primer_specificity_results.csv",
                        "text/csv",
                        key='download-specificity-csv'
                    )
                else:
                    st.success("✅ No potential amplification products found. These primers appear to be highly specific.")

# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2025 Yash Munnalal Gupta. All rights reserved.

For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
