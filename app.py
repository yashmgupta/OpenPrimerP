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

def find_primer_matches(sequence, primer, allow_mismatches=0):
    """Find all positions where primer matches in sequence."""
    matches = []
    sequence = str(sequence).upper()
    primer = str(primer).upper()
    
    for i in range(len(sequence) - len(primer) + 1):
        subseq = sequence[i:i+len(primer)]
        mismatches = sum(1 for a, b in zip(primer, subseq) if a != b)
        if mismatches <= allow_mismatches:
            matches.append(i)
    
    return matches

def simulate_pcr(forward_primer, reverse_primer, sequence, max_product_size=5000):
    """Simulate PCR to find potential amplification products."""
    results = []
    
    # Convert to string and uppercase
    sequence = str(sequence).upper()
    forward_primer = str(forward_primer).upper()
    reverse_primer = str(reverse_primer).upper()
    
    # Find forward primer binding sites
    forward_matches = find_primer_matches(sequence, forward_primer)
    
    # Find reverse primer binding sites (look for reverse complement)
    reverse_primer_rc = str(Seq(reverse_primer).reverse_complement()).upper()
    reverse_matches = find_primer_matches(sequence, reverse_primer_rc)
    
    # Find valid PCR products
    for f_pos in forward_matches:
        for r_pos in reverse_matches:
            # Reverse primer should be downstream of forward primer
            if r_pos > f_pos:
                product_size = r_pos + len(reverse_primer_rc) - f_pos
                if product_size <= max_product_size:
                    results.append({
                        'forward_start': f_pos,
                        'reverse_start': r_pos,
                        'product_size': product_size,
                        'amplicon_start': f_pos,
                        'amplicon_end': r_pos + len(reverse_primer_rc)
                    })
    
    return results

def check_primer_specificity(forward_primer, reverse_primer, records, target_record_id, target_feature=None, target_feature_name="Unknown"):
    """Check primer specificity across all records and features."""
    specificity_results = []
    target_found = False
    
    st.write(f"**Checking specificity for primers:**")
    st.write(f"- Forward: {forward_primer}")
    st.write(f"- Reverse: {reverse_primer}")
    st.write(f"- Target: {target_record_id} - {target_feature_name}")
    
    # Check each record
    for record in records:
        record_id = record.id
        st.write(f"Analyzing record: {record_id}...")
        
        # If this is the target record and we have a target feature, check it first
        if record_id == target_record_id and target_feature is not None:
            try:
                target_seq = target_feature.extract(record.seq)
                pcr_products = simulate_pcr(forward_primer, reverse_primer, target_seq)
                
                if pcr_products:
                    target_found = True
                    for i, product in enumerate(pcr_products):
                        result = {
                            'Record ID': record_id,
                            'Feature Type': 'Target Feature',
                            'Feature Name': target_feature_name,
                            'Location': str(target_feature.location),
                            'Product Size (bp)': product['product_size'],
                            'Amplicon Position': f"{product['amplicon_start']}-{product['amplicon_end']}",
                            'Is Target': 'YES',
                            'Status': '✅ Target'
                        }
                        specificity_results.append(result)
            except Exception as e:
                st.warning(f"Error checking target feature: {e}")
        
        # Check all features in the record
        features = extract_features_from_record(record)
        for feature_type in features:
            for feature in features[feature_type]:
                # Skip the target feature to avoid duplication
                if (record_id == target_record_id and target_feature is not None and 
                    feature.location == target_feature.location):
                    continue
                
                try:
                    # Get feature information
                    feature_name = feature.qualifiers.get('gene', 
                                    feature.qualifiers.get('product', 
                                    feature.qualifiers.get('locus_tag', ['Unknown'])))[0]
                    
                    # Extract feature sequence
                    feature_seq = feature.extract(record.seq)
                    
                    # Skip very short sequences
                    if len(feature_seq) < len(forward_primer) + len(reverse_primer):
                        continue
                    
                    # Simulate PCR
                    pcr_products = simulate_pcr(forward_primer, reverse_primer, feature_seq)
                    
                    # Add results
                    for product in pcr_products:
                        result = {
                            'Record ID': record_id,
                            'Feature Type': feature_type,
                            'Feature Name': feature_name,
                            'Location': str(feature.location),
                            'Product Size (bp)': product['product_size'],
                            'Amplicon Position': f"{product['amplicon_start']}-{product['amplicon_end']}",
                            'Is Target': 'NO',
                            'Status': '⚠️ Off-target'
                        }
                        specificity_results.append(result)
                        
                except Exception as e:
                    # Skip problematic features silently
                    continue
        
        # Also check the entire record sequence for comprehensive analysis
        try:
            whole_seq_products = simulate_pcr(forward_primer, reverse_primer, record.seq)
            for product in whole_seq_products:
                # Check if this product overlaps with any feature we already found
                overlaps_existing = False
                for existing in specificity_results:
                    if (existing['Record ID'] == record_id and 
                        abs(product['product_size'] - existing['Product Size (bp)']) < 10):
                        overlaps_existing = True
                        break
                
                if not overlaps_existing:
                    is_target_match = (record_id == target_record_id and not target_found)
                    result = {
                        'Record ID': record_id,
                        'Feature Type': 'Whole sequence',
                        'Feature Name': 'Non-annotated region',
                        'Location': f"{product['amplicon_start']}-{product['amplicon_end']}",
                        'Product Size (bp)': product['product_size'],
                        'Amplicon Position': f"{product['amplicon_start']}-{product['amplicon_end']}",
                        'Is Target': 'MAYBE' if is_target_match else 'NO',
                        'Status': '🔍 Unannotated' if is_target_match else '⚠️ Off-target'
                    }
                    specificity_results.append(result)
                    if is_target_match:
                        target_found = True
                        
        except Exception as e:
            st.warning(f"Could not analyze whole sequence for {record_id}: {e}")
    
    return specificity_results, target_found

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
                    gene_name = feature.qualifiers.get('gene', 
                                feature.qualifiers.get('product', 
                                feature.qualifiers.get('locus_tag', ['Unknown'])))[0]
                    feature_options.append(f"{gene_name} ({feature.location})")
                
                selected_index = st.selectbox(
                    f'Select a {feature_type}:', 
                    options=range(len(feature_options)), 
                    format_func=lambda x: feature_options[x],
                    help="Select a specific feature based on its gene name and location."
                )
                selected_feature = features[feature_type][selected_index]
                selected_feature_name = feature_options[selected_index].split(' (')[0]
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
                                'Left TM (°C)': round(primers.get(f'PRIMER_LEFT_{i}_TM', 0), 2),
                                'Right TM (°C)': round(primers.get(f'PRIMER_RIGHT_{i}_TM', 0), 2),
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
                        st.session_state['selected_feature_name'] = selected_feature_name
                        
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
                        
                        # Primer specificity check option
                        if 'primer_data' in st.session_state:
                            st.write("## Primer Specificity Check")
                            st.write("Check if your primers could amplify unintended regions.")
                            
                            primer_options = [f"Pair {i+1}: {p['Left Sequence'][:15]}... / {p['Right Sequence'][:15]}..." 
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
                                
                                with st.spinner("Analyzing primer specificity across all sequences..."):
                                    specificity_results, target_found = check_primer_specificity(
                                        forward_primer,
                                        reverse_primer,
                                        st.session_state['records'],
                                        st.session_state['selected_record'],
                                        st.session_state['selected_feature'],
                                        st.session_state['selected_feature_name']
                                    )
                                
                                # Display results
                                st.write("---")
                                if specificity_results:
                                    st.write(f"## Specificity Analysis Results")
                                    st.write(f"Found **{len(specificity_results)}** potential amplification products:")
                                    
                                    # Convert to DataFrame for display
                                    spec_df = pd.DataFrame(specificity_results)
                                    
                                    # Color coding function
                                    def highlight_specificity(row):
                                        if row['Is Target'] == 'YES':
                                            return ['background-color: #d4edda'] * len(row)  # Green
                                        elif row['Is Target'] == 'MAYBE':
                                            return ['background-color: #fff3cd'] * len(row)  # Yellow
                                        else:
                                            return ['background-color: #f8d7da'] * len(row)  # Red
                                    
                                    # Display the results table
                                    st.dataframe(spec_df.style.apply(highlight_specificity, axis=1), use_container_width=True)
                                    
                                    # Summary statistics
                                    target_count = sum(1 for r in specificity_results if r['Is Target'] in ['YES', 'MAYBE'])
                                    off_target_count = sum(1 for r in specificity_results if r['Is Target'] == 'NO')
                                    
                                    st.write("### Summary:")
                                    col1, col2, col3 = st.columns(3)
                                    with col1:
                                        st.metric("Target Amplifications", target_count)
                                    with col2:
                                        st.metric("Off-target Amplifications", off_target_count)
                                    with col3:
                                        if off_target_count == 0:
                                            st.success("✅ SPECIFIC")
                                        else:
                                            st.warning("⚠️ NON-SPECIFIC")
                                    
                                    # Detailed interpretation
                                    st.write("### Interpretation:")
                                    if target_count > 0 and off_target_count == 0:
                                        st.success("🎯 **Excellent specificity!** These primers should amplify only your target sequence.")
                                    elif target_count > 0 and off_target_count > 0:
                                        st.warning(f"⚠️ **Moderate specificity.** Primers will amplify the target but also {off_target_count} off-target sequence(s). Consider redesigning for better specificity.")
                                    elif target_count == 0 and off_target_count > 0:
                                        st.error("❌ **Poor specificity.** Primers do not amplify the intended target but will amplify off-target sequences.")
                                    elif target_count == 0 and off_target_count == 0:
                                        st.warning("❓ **No amplification predicted.** Primers may not work under standard PCR conditions.")
                                
                                else:
                                    st.success("🎯 **Perfect Specificity!**")
                                    st.write("No amplification products were found in any sequence other than your target.")
                                    st.write("These primers appear to be highly specific to your selected sequence.")
                    else:
                        st.error('No primers were found. Please adjust your parameters and try again.')
            else:
                st.warning(f"No {feature_type} features found in the selected record.")
    
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Please check that your GenBank file is valid and try again.")
        # Add debug information
        import traceback
        with st.expander("Debug Information"):
            st.code(traceback.format_exc())

# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2024 Yash Munnalal Gupta. All rights reserved.
For inquiries or permissions, contact: [yash.610@live.com](mailto:yash.610@live.com)
""", unsafe_allow_html=True)
