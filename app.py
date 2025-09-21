import streamlit as st
import os
import io
import logging
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.stats import chi2
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import plotly.graph_objects as go

st.set_page_config(page_title="GeneBank Genie", layout="wide")
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ---------------------------
# Core Classes and Functions
# ---------------------------
class GenBankParser:
    def __init__(self, text_io):
        self.text_io = text_io
        self.records = []

    def load_records(self):
        # Always seek to start (for re-reads)
        self.text_io.seek(0)
        self.records = list(SeqIO.parse(self.text_io, "genbank"))
        return self.records

class AnalysisEngine:
    @staticmethod
    def perform_pca(features, n_components=2):
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(features_scaled)
        return pca_result, pca.explained_variance_ratio_

    @staticmethod
    def duplicate_feature(X):
        return np.hstack([X, X])

    @staticmethod
    def compute_mahalanobis(points):
        mean_vec = np.mean(points, axis=0)
        cov_matrix = np.cov(points, rowvar=False)
        inv_cov_matrix = np.linalg.inv(cov_matrix)
        diff = points - mean_vec
        md_squared = np.sum(diff.dot(inv_cov_matrix) * diff, axis=1)
        return np.sqrt(md_squared)

def extract_gene_sequences(records, selected_gene):
    gene_seqs = []
    for rec in records:
        for feature in rec.features:
            if feature.type.lower() in ["gene", "cds"]:
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    if gene_name.upper() == selected_gene.upper():
                        seq = feature.extract(rec.seq)
                        new_rec = SeqRecord(seq, id=rec.id, description=rec.annotations.get("organism", ""))
                        gene_seqs.append(new_rec)
                        break
    return gene_seqs

def extract_full_sequences(records):
    full_seqs = []
    for rec in records:
        new_rec = SeqRecord(rec.seq, id=rec.id, description=rec.description)
        full_seqs.append(new_rec)
    return full_seqs

def save_fasta_download(sequences, file_label):
    fasta_io = io.StringIO()
    SeqIO.write(sequences, fasta_io, "fasta")
    fasta_io.seek(0)
    st.download_button(
        label=f"Download {file_label} as FASTA",
        data=fasta_io.getvalue(),
        file_name=f"{file_label}.fasta",
        mime="text/plain"
    )

# ---------------------------
# Streamlit UI
# ---------------------------

st.title("GeneBank Genie (Streamlit Version)")
st.markdown("""
A versatile tool for the analysis of GenBank records.  
**Version 1.0 – © Dr. Yash Munnalal Gupta**  
""")

with st.expander("ℹ️ About GeneBank Genie"):
    st.markdown("""
    **GeneBank Genie** is a comprehensive tool for analyzing GenBank files.  
    It includes modules for general analysis, gene analysis,  
    taxonomic visualization, additional visualizations, dendrogram analysis, and sequence extraction.
    """)

# Sidebar: GenBank File Upload
st.sidebar.header("Step 1. Upload GenBank File")
gb_file = st.sidebar.file_uploader("Upload a GenBank file (.gb, .gbk)", type=["gb", "gbk"])
if gb_file is not None:
    # Always wrap the file as text for Biopython
    gb_file.seek(0)
    gb_bytes = gb_file.read()
    if isinstance(gb_bytes, bytes):
        gb_text = gb_bytes.decode("utf-8")
    else:
        gb_text = gb_bytes
    gb_text_io = io.StringIO(gb_text)
    parser = GenBankParser(gb_text_io)
    try:
        records = parser.load_records()
    except Exception as e:
        st.error(f"Error parsing GenBank file: {e}")
        records = []
    else:
        st.sidebar.success(f"Loaded {len(records)} GenBank records.")
else:
    records = []
    st.warning("Upload a GenBank file to get started.")

if records:
    # Sidebar: Taxonomy Level Selection
    st.sidebar.header("Step 2. Set Analysis Options")
    tax_levels = []
    for rec in records:
        taxonomy = rec.annotations.get("taxonomy", [])
        if len(taxonomy) > len(tax_levels):
            tax_levels = taxonomy
    if tax_levels:
        max_level = len(tax_levels) - 1
        tax_level_index = st.sidebar.number_input("Taxonomy Level Index", min_value=0, max_value=max_level, value=max_level)
    else:
        tax_level_index = 0
    color_palette = st.sidebar.selectbox("Color Palette", ["tab20", "viridis", "plasma", "inferno", "magma", "cividis"])

    # Main Tabs
    tabs = st.tabs([
        "General Analysis", "Gene Analysis", "Sankey Diagram",
        "Additional Visualizations", "Dendrogram Analysis", "Sequence Extraction"
    ])

    # --- General Analysis ---
    with tabs[0]:
        st.header("General Analysis")
        # Populate tax groups for selection
        tax_groups = set()
        for rec in records:
            taxonomy = rec.annotations.get("taxonomy", [])
            if len(taxonomy) > tax_level_index:
                tax_groups.add(taxonomy[tax_level_index])
            elif taxonomy:
                tax_groups.add(taxonomy[-1])
        tax_groups = sorted(list(tax_groups))
        selected_tax_group = st.selectbox("Selected Taxonomic Group", tax_groups)
        default_marker = st.text_input("Default Marker", "o")
        outlier_marker = st.text_input("Outlier Marker", "D")
        run_general = st.button("Run General Analysis")
        summary_general = st.empty()
        if run_general:
            # Compute features
            data = []
            for rec in records:
                seq = str(rec.seq).upper()
                total_length = len(seq)
                if total_length == 0:
                    continue
                countA = seq.count("A")
                countC = seq.count("C")
                countG = seq.count("G")
                countT = seq.count("T")
                propA = countA / total_length
                propC = countC / total_length
                propG = countG / total_length
                propT = countT / total_length
                gc_content = (countG + countC) / total_length
                gene_count = sum(1 for f in rec.features if f.type.lower() == 'gene')
                taxonomy = rec.annotations.get("taxonomy", [])
                if len(taxonomy) > tax_level_index:
                    chosen_tax = taxonomy[tax_level_index]
                elif taxonomy:
                    chosen_tax = taxonomy[-1]
                else:
                    chosen_tax = "Unknown"
                organism = rec.annotations.get("organism", "Unknown")
                data.append([propA, propC, propG, propT, total_length, gc_content, gene_count,
                             organism, " | ".join(taxonomy), chosen_tax])
            df = pd.DataFrame(data, columns=['A','C','G','T','SeqLength','GC','GeneCount',
                                             'Organism','Full_Taxonomy','TaxLevel'])

            # PCA
            features_nuc = df[['A','C','G','T']].values
            features_add = df[['SeqLength','GC','GeneCount']].values
            features_combined = df[['A','C','G','T','SeqLength','GC','GeneCount']].values
            features_seq = df[['SeqLength']].values
            features_gc = df[['GC']].values
            features_gene = df[['GeneCount']].values

            pca_nuc, var_nuc = AnalysisEngine.perform_pca(features_nuc)
            pca_add, var_add = AnalysisEngine.perform_pca(features_add)
            pca_comb, var_comb = AnalysisEngine.perform_pca(features_combined)
            pca_seq, var_seq = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_seq))
            pca_gc, var_gc = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_gc))
            pca_gene, var_gene = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_gene))

            df['PC1_nuc'] = pca_nuc[:, 0]; df['PC2_nuc'] = pca_nuc[:, 1]
            df['PC1_add'] = pca_add[:, 0]; df['PC2_add'] = pca_add[:, 1]
            df['PC1_comb'] = pca_comb[:, 0]; df['PC2_comb'] = pca_comb[:, 1]
            df['PC1_seq'] = pca_seq[:, 0]; df['PC2_seq'] = pca_seq[:, 1]
            df['PC1_gc'] = pca_gc[:, 0]; df['PC2_gc'] = pca_gc[:, 1]
            df['PC1_gene'] = pca_gene[:, 0]; df['PC2_gene'] = pca_gene[:, 1]

            # Outlier detection
            group_mask = df['TaxLevel'] == selected_tax_group
            df.loc[:, 'Outlier'] = False
            group_df = df[group_mask]
            outlier_txt = ""
            if not group_df.empty:
                group_points = group_df[['PC1_comb', 'PC2_comb']].values
                distances = AnalysisEngine.compute_mahalanobis(group_points)
                threshold = np.sqrt(chi2.ppf(0.95, df=2))
                outlier_flags = distances > threshold
                df.loc[group_mask, 'Outlier'] = outlier_flags
                outlier_species = df[(df['TaxLevel'] == selected_tax_group) & (df['Outlier'] == True)]['Organism']
                if not outlier_species.empty:
                    unique_species = outlier_species.unique()
                    outlier_txt = f"Outlier species in {selected_tax_group}: {', '.join(unique_species)}"
                else:
                    outlier_txt = f"No outlier species detected in {selected_tax_group}."
            else:
                outlier_txt = f"No records found for {selected_tax_group}."
            summary_general.info(outlier_txt)

            st.subheader("PCA Plots")
            fig, axes = plt.subplots(2, 3, figsize=(20, 10))
            axes = axes.flatten()
            unique_tax = df['TaxLevel'].unique()
            cmap = plt.cm.get_cmap(color_palette, len(unique_tax))
            color_dict = {tax: cmap(i) for i, tax in enumerate(unique_tax)}
            for i, (col1, col2, title, var, pca_data) in enumerate([
                ('PC1_nuc', 'PC2_nuc', f"Nucleotide PCA (var: {var_nuc[0]:.2f}, {var_nuc[1]:.2f})", var_nuc, pca_nuc),
                ('PC1_seq', 'PC2_seq', f"SeqLength PCA (var: {var_seq[0]:.2f}, {var_seq[1]:.2f})", var_seq, pca_seq),
                ('PC1_gc', 'PC2_gc', f"GC PCA (var: {var_gc[0]:.2f}, {var_gc[1]:.2f})", var_gc, pca_gc),
                ('PC1_gene', 'PC2_gene', f"Gene Count PCA (var: {var_gene[0]:.2f}, {var_gene[1]:.2f})", var_gene, pca_gene),
                ('PC1_add', 'PC2_add', f"Additional PCA (var: {var_add[0]:.2f}, {var_add[1]:.2f})", var_add, pca_add),
                ('PC1_comb', 'PC2_comb', f"Combined PCA (var: {var_comb[0]:.2f}, {var_comb[1]:.2f})", var_comb, pca_comb)
            ]):
                for tax in unique_tax:
                    subset = df[df['TaxLevel'] == tax]
                    marker = outlier_marker if (col1 == 'PC1_comb' and tax == selected_tax_group and subset['Outlier'].any()) else default_marker
                    axes[i].scatter(subset[col1], subset[col2], label=tax, alpha=0.7, marker=marker)
                axes[i].set_title(title)
                axes[i].set_xlabel("PC1")
                axes[i].set_ylabel("PC2")
                axes[i].legend(fontsize=8)
            st.pyplot(fig)
            st.dataframe(df.head(), use_container_width=True)
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button("Download PCA Features CSV", data=csv, file_name="pca_features.csv", mime='text/csv')

    # --- Gene Analysis ---
    with tabs[1]:
        st.header("Gene Analysis")
        # Get all gene names
        gene_set = set()
        for rec in records:
            for feature in rec.features:
                if feature.type.lower() in ["gene", "cds"]:
                    if "gene" in feature.qualifiers:
                        gene_set.add(feature.qualifiers["gene"][0].upper())
        gene_list = sorted(list(gene_set))
        selected_gene = st.selectbox("Select Gene", gene_list)
        selected_tax_group_gene = st.text_input("Selected Taxonomic Group (for outlier detection)", tax_groups[0] if tax_groups else "")
        run_gene = st.button("Run Gene Analysis")
        summary_gene = st.empty()
        if run_gene:
            gene_data = []
            for rec in records:
                gene_seq = None
                for feature in rec.features:
                    if feature.type.lower() in ["gene", "cds"]:
                        if "gene" in feature.qualifiers:
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name.upper() == selected_gene.upper():
                                gene_seq = str(feature.extract(rec.seq)).upper()
                                break
                if gene_seq is None or len(gene_seq) == 0:
                    continue
                gene_length = len(gene_seq)
                countA = gene_seq.count("A")
                countC = gene_seq.count("C")
                countG = gene_seq.count("G")
                countT = gene_seq.count("T")
                propA = countA / gene_length
                propC = countC / gene_length
                propG = countG / gene_length
                propT = countT / gene_length
                gc_content = (countG + countC) / gene_length
                taxonomy = rec.annotations.get("taxonomy", [])
                if len(taxonomy) > tax_level_index:
                    chosen_tax = taxonomy[tax_level_index]
                elif taxonomy:
                    chosen_tax = taxonomy[-1]
                else:
                    chosen_tax = "Unknown"
                organism = rec.annotations.get("organism", "Unknown")
                gene_data.append([propA, propC, propG, propT, gene_length, gc_content,
                                  organism, " | ".join(taxonomy), chosen_tax])
            df_gene = pd.DataFrame(gene_data, columns=["A", "C", "G", "T", "GeneLength", "GC",
                                                       "Organism", "Full_Taxonomy", "TaxLevel"])
            features = df_gene[["A", "C", "G", "T", "GeneLength", "GC"]].values
            pca_result, var_gene = AnalysisEngine.perform_pca(features)
            df_gene["PC1"] = pca_result[:, 0]
            df_gene["PC2"] = pca_result[:, 1]
            group_mask = df_gene["TaxLevel"] == selected_tax_group_gene
            df_gene["Outlier"] = False
            group_df = df_gene[group_mask]
            outlier_txt = ""
            if not group_df.empty:
                group_points = group_df[["PC1", "PC2"]].values
                distances = AnalysisEngine.compute_mahalanobis(group_points)
                threshold = np.sqrt(chi2.ppf(0.95, df=2))
                outlier_flags = distances > threshold
                df_gene.loc[group_mask, "Outlier"] = outlier_flags
                outlier_species = df_gene[(df_gene["TaxLevel"] == selected_tax_group_gene) & (df_gene["Outlier"] == True)]["Organism"]
                if not outlier_species.empty:
                    unique_species = outlier_species.unique()
                    outlier_txt = f"Outlier species in {selected_tax_group_gene}: {', '.join(unique_species)}"
                else:
                    outlier_txt = f"No outlier species detected in {selected_tax_group_gene}."
            else:
                outlier_txt = f"No records found for {selected_tax_group_gene}."
            summary_gene.info(outlier_txt)
            st.subheader("Gene PCA Scatter Plot")
            fig, ax = plt.subplots(figsize=(8, 6))
            unique_tax = df_gene["TaxLevel"].unique()
            cmap = plt.cm.get_cmap(color_palette, len(unique_tax))
            color_dict = {tax: cmap(i) for i, tax in enumerate(unique_tax)}
            for tax in unique_tax:
                subset = df_gene[df_gene["TaxLevel"] == tax]
                marker = outlier_marker if (tax == selected_tax_group_gene and subset["Outlier"].any()) else default_marker
                ax.scatter(subset["PC1"], subset["PC2"], label=tax, alpha=0.7, marker=marker)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.legend(fontsize=8)
            st.pyplot(fig)
            st.dataframe(df_gene.head(), use_container_width=True)
            csv_gene = df_gene.to_csv(index=False).encode('utf-8')
            st.download_button("Download Gene PCA Features CSV", data=csv_gene, file_name="selected_gene_pca_features.csv", mime='text/csv')

    # --- Sankey Diagram ---
    with tabs[2]:
        st.header("Sankey Diagram")
        start_level = st.number_input("Start Level (for taxonomy)", min_value=0, value=0)
        run_sankey = st.button("Run Sankey Diagram")
        if run_sankey:
            import collections
            link_counts = collections.defaultdict(int)
            nodes_set = set()
            for rec in records:
                taxonomy = rec.annotations.get("taxonomy", [])
                if len(taxonomy) <= start_level:
                    continue
                for i in range(start_level, len(taxonomy) - 1):
                    source_node = f"L{i}-{taxonomy[i]}"
                    target_node = f"L{i+1}-{taxonomy[i+1]}"
                    nodes_set.add(source_node)
                    nodes_set.add(target_node)
                    link_counts[(source_node, target_node)] += 1
            def sort_key(label):
                try:
                    level = int(label.split("-")[0][1:])
                except:
                    level = 9999
                return (level, label)
            nodes_list = sorted(list(nodes_set), key=sort_key)
            node_to_index = {node: i for i, node in enumerate(nodes_list)}
            sources = []
            targets = []
            values = []
            for (src, tgt), count in link_counts.items():
                sources.append(node_to_index[src])
                targets.append(node_to_index[tgt])
                values.append(count)
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=nodes_list,
                    color="blue"
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values
                ))])
            fig.update_layout(title_text=f"Sankey Diagram of Taxonomy Flow (Starting at Level {start_level})", font_size=10)
            st.plotly_chart(fig, use_container_width=True)

    # --- Additional Visualizations ---
    with tabs[3]:
        st.header("Additional Visualizations")
        if st.button("Show Correlation Heatmap (A, C, G, T, SeqLength, GC, GeneCount)"):
            data = []
            for rec in records:
                seq = str(rec.seq).upper()
                total_length = len(seq)
                if total_length == 0:
                    continue
                countA = seq.count("A")
                countC = seq.count("C")
                countG = seq.count("G")
                countT = seq.count("T")
                propA = countA / total_length
                propC = countC / total_length
                propG = countG / total_length
                propT = countT / total_length
                gc_content = (countG + countC) / total_length
                gene_count = sum(1 for f in rec.features if f.type.lower() == 'gene')
                taxonomy = rec.annotations.get("taxonomy", [])
                organism = rec.annotations.get("organism", "Unknown")
                data.append([propA, propC, propG, propT, total_length, gc_content, gene_count, organism])
            df = pd.DataFrame(data, columns=['A','C','G','T','SeqLength','GC','GeneCount','Organism'])
            plt.figure(figsize=(8, 6))
            corr_matrix = df[['A','C','G','T','SeqLength','GC','GeneCount']].corr()
            sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
            st.pyplot(plt.gcf())
        n_clusters = st.number_input("Number of Clusters for KMeans", min_value=2, value=2)
        if st.button("Run KMeans Clustering"):
            data = []
            for rec in records:
                seq = str(rec.seq).upper()
                total_length = len(seq)
                if total_length == 0:
                    continue
                countA = seq.count("A")
                countC = seq.count("C")
                countG = seq.count("G")
                countT = seq.count("T")
                propA = countA / total_length
                propC = countC / total_length
                propG = countG / total_length
                propT = countT / total_length
                gc_content = (countG + countC) / total_length
                gene_count = sum(1 for f in rec.features if f.type.lower() == 'gene')
                taxonomy = rec.annotations.get("taxonomy", [])
                organism = rec.annotations.get("organism", "Unknown")
                data.append([propA, propC, propG, propT, total_length, gc_content, gene_count, organism])
            df = pd.DataFrame(data, columns=['A','C','G','T','SeqLength','GC','GeneCount','Organism'])
            features_combined = df[['A','C','G','T','SeqLength','GC','GeneCount']].values
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            clusters = kmeans.fit_predict(features_combined)
            df['KMeans_Cluster'] = clusters
            sil_score = silhouette_score(features_combined, clusters)
            st.info(f"Silhouette score: {sil_score:.2f}")
            st.write("KMeans Cluster Assignment (first 10 rows):")
            st.dataframe(df[['Organism','KMeans_Cluster']].head(10))
            # PCA
            pca_comb, _ = AnalysisEngine.perform_pca(features_combined)
            df['PC1_comb'] = pca_comb[:, 0]
            df['PC2_comb'] = pca_comb[:, 1]
            plt.figure(figsize=(8, 6))
            sns.scatterplot(x='PC1_comb', y='PC2_comb', hue='KMeans_Cluster', data=df, palette='tab10', s=80)
            plt.title("K-Means Clustering (Combined PCA Projection)")
            st.pyplot(plt.gcf())
        if st.button("Show Pair Plot (A, C, G, T)"):
            data = []
            for rec in records:
                seq = str(rec.seq).upper()
                total_length = len(seq)
                if total_length == 0:
                    continue
                countA = seq.count("A")
                countC = seq.count("C")
                countG = seq.count("G")
                countT = seq.count("T")
                propA = countA / total_length
                propC = countC / total_length
                propG = countG / total_length
                propT = countT / total_length
                data.append([propA, propC, propG, propT])
            df = pd.DataFrame(data, columns=['A','C','G','T'])
            sns.pairplot(df, diag_kind='kde')
            st.pyplot(plt.gcf())

    # --- Dendrogram Analysis ---
    with tabs[4]:
        st.header("Dendrogram Analysis")
        feature_option = st.selectbox(
            "Select Feature Set for Dendrogram",
            ["Nucleotide Composition", "Additional Features", "Combined Features"]
        )
        if st.button("Run Dendrogram"):
            data = []
            for rec in records:
                seq = str(rec.seq).upper()
                total_length = len(seq)
                if total_length == 0:
                    continue
                countA = seq.count("A")
                countC = seq.count("C")
                countG = seq.count("G")
                countT = seq.count("T")
                propA = countA / total_length
                propC = countC / total_length
                propG = countG / total_length
                propT = countT / total_length
                gc_content = (countG + countC) / total_length
                gene_count = sum(1 for f in rec.features if f.type.lower() == 'gene')
                taxonomy = rec.annotations.get("taxonomy", [])
                if len(taxonomy) > tax_level_index:
                    chosen_tax = taxonomy[tax_level_index]
                elif taxonomy:
                    chosen_tax = taxonomy[-1]
                else:
                    chosen_tax = "Unknown"
                data.append([propA, propC, propG, propT, total_length, gc_content, gene_count, chosen_tax])
            df = pd.DataFrame(data, columns=['A','C','G','T','SeqLength','GC','GeneCount','TaxLevel'])
            if feature_option == "Nucleotide Composition":
                features = df[['A', 'C', 'G', 'T']].values
            elif feature_option == "Additional Features":
                features = df[['SeqLength', 'GC', 'GeneCount']].values
            else:
                features = df[['A', 'C', 'G', 'T', 'SeqLength', 'GC', 'GeneCount']].values
            scaler = StandardScaler()
            features_scaled = scaler.fit_transform(features)
            linkage_matrix = sch.linkage(features_scaled, method='ward')
            plt.figure(figsize=(12, 8))
            sch.dendrogram(linkage_matrix, labels=df['TaxLevel'].values, leaf_rotation=90)
            plt.title("Dendrogram Analysis")
            plt.xlabel("Taxonomic Level")
            plt.ylabel("Euclidean Distance")
            st.pyplot(plt.gcf())

    # --- Sequence Extraction ---
    with tabs[5]:
        st.header("Sequence Extraction")
        st.subheader("Full Sequences")
        if st.button("Download All Full Sequences as FASTA"):
            full_seqs = extract_full_sequences(records)
            save_fasta_download(full_seqs, "full_sequences")
        st.subheader("Gene Sequences")
        gene_list = sorted(list({feature.qualifiers["gene"][0].upper()
                                 for rec in records for feature in rec.features
                                 if feature.type.lower() in ["gene", "cds"] and "gene" in feature.qualifiers}))
        selected_gene_extract = st.selectbox("Select Gene for Extraction", gene_list)
        if st.button("Download Selected Gene Sequences as FASTA"):
            gene_seqs = extract_gene_sequences(records, selected_gene_extract)
            if gene_seqs:
                save_fasta_download(gene_seqs, f"{selected_gene_extract}_extracted")
            else:
                st.warning(f"No sequences found for gene {selected_gene_extract}.")
