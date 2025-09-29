import streamlit as st
import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os
from io import StringIO

# Make sure we have a temp folder for saving little intermediate files
temp_dir = "temp"
os.makedirs(temp_dir, exist_ok=True)

# Streamlit UI setup
st.set_page_config(page_title="PCR Primer Design", page_icon="🧬", layout="wide")

# Sidebar user guide – quick “how-to” for anyone opening the app
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

# File uploader with additional help (we only care about GenBank format here)
uploaded_file = st.file_uploader(
    "Upload a GenBank file", 
    type=['gb', 'gbk'],
    help="Upload a GenBank (.gb or .gbk) file containing the DNA sequence from which to design primers."
)

# Some CSS tweaks to make things prettier
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

def simple_pcr_check(forward_primer, reverse_primer, sequence):
    """Simple PCR amplification check - just look for primer binding sites."""
    try:
        # Convert everything to uppercase strings
        seq_str = str(sequence).upper()
        forward = str(forward_primer).upper()
        reverse = str(reverse_primer).upper()
        reverse_rc = str(Seq(reverse).reverse_complement()).upper()
        
        # Find primer positions
        forward_pos = seq_str.find(forward)
        reverse_pos = seq_str.find(reverse_rc)
        
        # Check if both primers bind and in correct orientation
        if forward_pos != -1 and reverse_pos != -1 and reverse_pos > forward_pos:
            product_size = reverse_pos - forward_pos + len(reverse)
            return [{
                'start': forward_pos,
                'end': reverse_pos + len(reverse_rc),
                'size': product_size
            }]
        return []
    except Exception as e:
        return []

def count_total_features(records):
    """Count total features across all records for progress tracking."""
    total = 0
    for record in records:
        features = extract_features_from_record(record)
        for feature_type in ['CDS', 'tRNA', 'gene']:
            if feature_type in features:
                total += len(features[feature_type])
    return total

def check_primer_specificity_with_progress(forward_primer, reverse_primer, records, target_record_id, target_feature, target_feature_name):
    """Primer specificity check with progress bar."""
    
    # Initialize results
    all_results = []
    
    try:
        # Step 1: Check target feature
        st.write("🎯 **Checking target sequence...**")
        
        if target_feature is not None:
            target_record = None
            for record in records:
                if record.id == target_record_id:
                    target_record = record
                    break
            
            if target_record:
                target_seq = target_feature.extract(target_record.seq)
                pcr_results = simple_pcr_check(forward_primer, reverse_primer, target_seq)
                
                if pcr_results:
                    for result in pcr_results:
                        all_results.append({
                            'Record': target_record_id,
                            'Feature': f"{target_feature_name} (TARGET)",
                            'Location': str(target_feature.location),
                            'Product Size': result['size'],
                            'Status': '🎯 TARGET',
                            'Is_Target': True
                        })
                    st.success(f"✅ Target sequence will be amplified (Product size: {pcr_results[0]['size']} bp)")
                else:
                    st.warning("⚠️ No amplification detected in target sequence")
        
        # Step 2: Check all other features with progress bar
        st.write("🔍 **Checking for off-target amplification...**")
        
        # Count total features for progress tracking
        total_features = count_total_features(records)
        
        # Create progress bar
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        features_checked = 0
        
        for record in records:
            # Get all features
            features = extract_features_from_record(record)
            
            for feature_type in ['CDS', 'tRNA', 'gene']:
                if feature_type in features:
                    for feature in features[feature_type]:
                        # Update progress
                        features_checked += 1
                        progress = features_checked / total_features
                        progress_bar.progress(progress)
                        status_text.text(f"Analyzing {record.id}: {feature_type} features... ({features_checked}/{total_features})")
                        
                        # Skip target feature to avoid duplication
                        if (record.id == target_record_id and target_feature is not None and 
                            str(feature.location) == str(target_feature.location)):
                            continue
                        
                        try:
                            # Get feature name
                            gene_name = (feature.qualifiers.get('gene', [None])[0] or 
                                        feature.qualifiers.get('product', [None])[0] or 
                                        feature.qualifiers.get('locus_tag', ['Unknown'])[0])
                            
                            # Extract sequence
                            feature_seq = feature.extract(record.seq)
                            
                            # Skip very short sequences
                            if len(feature_seq) < len(forward_primer) + len(reverse_primer) + 10:
                                continue
                            
                            # Check for PCR amplification
                            pcr_results = simple_pcr_check(forward_primer, reverse_primer, feature_seq)
                            
                            if pcr_results:
                                for result in pcr_results:
                                    all_results.append({
                                        'Record': record.id,
                                        'Feature': f"{gene_name} ({feature_type})",
                                        'Location': str(feature.location),
                                        'Product Size': result['size'],
                                        'Status': '⚠️ OFF-TARGET',
                                        'Is_Target': False
                                    })
                            
                        except Exception as e:
                            continue
        
        # Complete progress bar
        progress_bar.progress(1.0)
        status_text.text(f"✅ Analysis complete! Checked {total_features} features across {len(records)} record(s)")
        
        return all_results
        
    except Exception as e:
        st.error(f"Error during specificity analysis: {str(e)}")
        return []

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
                    gene_name = (feature.qualifiers.get('gene', [None])[0] or 
                                feature.qualifiers.get('product', [None])[0] or 
                                feature.qualifiers.get('locus_tag', ['Unknown'])[0])
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
                st.code(str(feature_sequence)[:200] + "..." if len(str(feature_sequence)) > 200 else str(feature_sequence), language="text")
                
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
                        # Store in session state
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
                    else:
                        st.error('No primers were found. Please adjust your parameters and try again.')
                
                # Primer specificity check section - show always after primers are designed
                if 'primer_data' in st.session_state and st.session_state['primer_data']:
                    st.write("---")
                    st.write("## 🎯 Primer Specificity Check")
                    st.info("This tool will check if your primers might amplify unintended sequences.")
                    
                    # Show primer options
                    selected_pair = st.selectbox(
                        "Select a primer pair to analyze:", 
                        options=range(len(st.session_state['primer_data'])),
                        format_func=lambda x: f"Pair {x+1}",
                        key="primer_pair_selection"
                    )
                    
                    # Show selected primer details
                    selected_primer = st.session_state['primer_data'][selected_pair]
                    st.write(f"**Selected Primer Pair {selected_pair + 1}:**")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write(f"Forward: `{selected_primer['Left Sequence']}`")
                    with col2:
                        st.write(f"Reverse: `{selected_primer['Right Sequence']}`")
                    
                    # Specificity check button
                    if st.button("🔍 Run Specificity Analysis", key="specificity_button"):
                        st.write("---")
                        
                        # Run the analysis
                        forward_primer = selected_primer['Left Sequence']
                        reverse_primer = selected_primer['Right Sequence']
                        
                        try:
                            # Run specificity check with progress bar
                            results = check_primer_specificity_with_progress(
                                forward_primer,
                                reverse_primer,
                                st.session_state['records'],
                                st.session_state['selected_record'],
                                st.session_state['selected_feature'],
                                st.session_state['selected_feature_name']
                            )
                            
                            st.write("---")
                            
                            # Display results
                            if results:
                                st.write(f"## 📊 Specificity Results ({len(results)} amplifications found)")
                                
                                # Create results dataframe
                                df = pd.DataFrame(results)
                                
                                # Display with color coding
                                def color_rows(row):
                                    if row['Is_Target']:
                                        return ['background-color: #d4edda'] * len(row)  # Light green
                                    else:
                                        return ['background-color: #f8d7da'] * len(row)  # Light red
                                
                                st.dataframe(df.style.apply(color_rows, axis=1), use_container_width=True)
                                
                                # Summary
                                target_count = sum(1 for r in results if r['Is_Target'])
                                offtarget_count = len(results) - target_count
                                
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Target Amplifications", target_count)
                                with col2:
                                    st.metric("Off-target Amplifications", offtarget_count)
                                with col3:
                                    if offtarget_count == 0 and target_count > 0:
                                        st.success("✅ SPECIFIC")
                                    elif offtarget_count > 0:
                                        st.warning("⚠️ NON-SPECIFIC")
                                    else:
                                        st.error("❌ NO TARGET")
                                
                                # Interpretation
                                st.write("### 📝 Interpretation:")
                                if target_count > 0 and offtarget_count == 0:
                                    st.success("🎯 **Excellent!** These primers are specific to your target sequence.")
                                elif target_count > 0 and offtarget_count > 0:
                                    st.warning(f"⚠️ **Caution:** Primers will amplify the target but also {offtarget_count} off-target sequence(s).")
                                elif target_count == 0 and offtarget_count > 0:
                                    st.error("❌ **Problem:** Primers don't amplify the target but will amplify off-target sequences.")
                                else:
                                    st.info("❓ **No amplification detected.** Primers may not work under standard conditions.")
                            
                            else:
                                st.success("🎯 **Perfect Specificity!**")
                                st.write("No off-target amplification products detected.")
                                st.write("These primers appear to be highly specific to your target sequence.")
                        
                        except Exception as e:
                            st.error(f"Error during specificity analysis: {str(e)}")
                            import traceback
                            with st.expander("Full Error Details"):
                                st.code(traceback.format_exc())
                
            else:
                st.warning(f"No {feature_type} features found in the selected record.")
    
    except Exception as e:
        st.error(f"An error occurred while processing the file: {str(e)}")
        import traceback
        with st.expander("Debug Information"):
            st.code(traceback.format_exc())

else:
    st.info("👆 Please upload a GenBank file to get started.")

# Add copyright information section at the end of the main page
st.markdown("""
---
**Copyright Notice**: © 2025 Yash Munnalal Gupta. All rights reserved.
For inquiries or permissions, contact: [yashmunnalalg@nu.ac.th ](mailto:yashmunnalalg@nu.ac.th)
""", unsafe_allow_html=True)
