import streamlit as st
import pandas as pd
import numpy as np
import base64
from Bio import Entrez, SeqIO
import io
import json
import time

# Set your email address for NCBI Entrez API (mandatory)
Entrez.email = "kaletejal05@mail.com"  # Replace with your actual email address

# Page configuration
st.set_page_config(
    page_title="Annotrax",
    page_icon="🧬",
    layout="wide"
)

# App styling
st.markdown("""
<style>
    body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    .main-header {
        font-size: 2.8rem;
        color: #1E6091;
        text-align: center;
        margin-bottom: 0px;
        font-weight: 600;
    }
    .tagline {
        font-size: 1.2rem;
        color: #168AAD;
        text-align: center;
        margin-top: 5px;
        margin-bottom: 25px;
    }
    .sub-header { /* General sub-header for tabs */
        font-size: 1.8rem; /* Main section title within a tab */
        color: #1E6091;
        font-weight: 500;
        margin-top: 10px;
        margin-bottom: 15px;
        border-bottom: 2px solid #e0e0e0;
        padding-bottom: 5px;
    }
    .about-section-header { /* Specific for sub-sections within "About" tab */
        font-size: 1.4rem;
        color: #1A5276;
        font-weight: 500;
        margin-top: 25px;
        margin-bottom: 10px;
    }
    .result-box {
        background-color: #f4f8fb;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
        border: 1px solid #d0e0eb;
    }
    .info-text {
        color: #31708f;
        font-size: 0.9rem;
    }
    .acknowledgement-highlight {
        color: #C58940;
        font-weight: bold;
    }
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 5px;
        border-bottom: 2px solid #1E6091;
    }
    .stTabs [data-baseweb="tab"] {
        height: auto;
        min-height: 45px;
        white-space: normal;
        background-color: #e9f1f7;
        border-radius: 6px 6px 0px 0px;
        padding: 10px 15px;
        font-weight: 500;
        color: #1E6091;
        border: 1px solid #d0e0eb;
        border-bottom: none;
    }
    .stTabs [aria-selected="true"] {
        background-color: #BDE0FE;
        color: #0d3b5c;
        border-color: #1E6091;
    }

    /* Developer's Desk specific styling */
    .developer-desk-container {
        padding: 20px;
        background-color: #fdfdff;
        border-radius: 8px;
    }
    .developer-profile-image {
        border-radius: 15px !important;
        border: 3px solid #BDE0FE !important;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1) !important;
        max-width: 100%; 
        height: auto;   
        display: block; 
    }
    .developer-intro-text {
        font-size: 1.05rem; /* Slightly adjusted for readability */
        line-height: 1.75; /* Slightly increased for better spacing */
        color: #2c3e50;
        text-align: justify;
        margin-top: 15px; /* Add space above the intro text */
    }
    .developer-name-caption { /* For the name under the image */
        text-align: center; 
        font-weight: 500; 
        color: #1E6091; 
        margin-top: 10px; /* Space between image and name */
        font-size: 1.1rem;
    }
    .developer-role {
        font-size: 1.1rem;
        color: #168AAD;
        margin-bottom: 10px; /* Space below role */
        text-align: left; /* Ensure role is left-aligned */
    }
     .developer-signature {
        text-align: right;
        font-style: italic;
        color: #555;
        margin-top: 20px;
        font-size: 0.95rem;
    }


    /* Content styling for About tab sections */
    .about-content-wrapper {
        font-size: 1.05rem;
        line-height: 1.7;
        color: #333;
        text-align: justify;
    }
    .about-content-wrapper ul, .about-content-wrapper ol {
        padding-left: 25px;
        list-style-position: outside;
    }
    .about-content-wrapper li {
        margin-bottom: 8px;
    }
    .about-content-wrapper p {
        margin-bottom: 12px;
    }
    .about-content-wrapper a {
        color: #1E6091;
        text-decoration: none;
    }
    .about-content-wrapper a:hover {
        text-decoration: underline;
    }
    .stExpander {
        border: 1px solid #dcecf8;
        border-radius: 6px;
        background-color: #f9fcff;
        margin-top: 15px;
    }
    .stExpander header {
        font-size: 1.2rem !important;
        font-weight: 500 !important;
        color: #1A5276 !important;
        background-color: #eef6fc;
        border-bottom: 1px solid #dcecf8;
        border-radius: 6px 6px 0 0;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 class='main-header'>🧬 Annotrax</h1>", unsafe_allow_html=True)
st.markdown("<p class='tagline'>Annotating Genes with Computational Excellence</p>", unsafe_allow_html=True)

# --- NCBI Data Fetching Functions (Keep as is) ---
@st.cache_data(ttl=3600)
def fetch_gene_annotation(gene_name, organism="Homo sapiens"):
    try:
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return None
        results = []
        for gene_id in record["IdList"]:
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            docsum = summary_record['DocumentSummarySet']['DocumentSummary'][0]
            nuccore_ids = []
            try:
                link_handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id)
                link_results = Entrez.read(link_handle)
                link_handle.close()
                if link_results[0]["LinkSetDb"]:
                    for link_info in link_results[0]["LinkSetDb"]:
                        if link_info["DbTo"] == "nuccore":
                            for link in link_info["Link"]: nuccore_ids.append(link["Id"])
                            break
            except Exception: pass
            genomic_info = docsum.get("GenomicInfo", [{}])
            chromosome = genomic_info[0].get("ChrLoc", "N/A") if genomic_info else "N/A"
            chr_start = genomic_info[0].get("ChrStart", "N/A") if genomic_info else "N/A"
            chr_stop = genomic_info[0].get("ChrStop", "N/A") if genomic_info else "N/A"
            gene_data = {
                "Gene ID": gene_id, "Gene Symbol": docsum.get("Name", gene_name),
                "Official Name": docsum.get("Description", "N/A"),
                "Aliases": ", ".join(docsum.get("OtherAliases", "").split(", ")[:5]),
                "Function": docsum.get("Summary", "No function description available"),
                "Chromosome": chromosome,
                "Position": f"{chr_start}-{chr_stop}" if chr_start != "N/A" and chr_stop != "N/A" else "N/A",
                "Exons": docsum.get("GenomicInfo", [{}])[0].get("ExonCount", "N/A") if genomic_info and docsum.get("GenomicInfo")[0] else "N/A",
                "Source Organism": organism,
                "Annotation Status": "Complete" if len(docsum.get("Summary", "")) > 50 else "Partial",
                "Nucleotide IDs": ", ".join(nuccore_ids[:3]) if nuccore_ids else "N/A"
            }
            results.append(gene_data)
        return results
    except Exception as e:
        st.error(f"Error fetching from NCBI: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_gene_sequence(nuccore_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        features = []
        for feature in record.features:
            if feature.type in ["CDS", "gene", "exon", "mRNA"]:
                feature_info = {
                    "type": feature.type, "location": str(feature.location),
                    "qualifiers": {k: str(v[0]) if isinstance(v, list) and v else str(v) for k,v in feature.qualifiers.items()}
                }
                features.append(feature_info)
        sequence_data = {
            "id": record.id, "name": record.name, "description": record.description,
            "length": len(record.seq), "features": features[:10],
            "sequence": str(record.seq)[:1000] + ("..." if len(record.seq) > 1000 else "")
        }
        return sequence_data
    except Exception as e:
        st.error(f"Error fetching GenBank sequence: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_fasta_sequence_data(nuccore_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="fasta", retmode="text")
        fasta_sequence = handle.read()
        handle.close()
        return fasta_sequence
    except Exception as e:
        st.error(f"Error fetching FASTA sequence: {e}")
        return None

def extract_feature_data(gene_data, feature_type):
    result = {}
    if feature_type == "function": result["Function"] = gene_data.get("Function", "")
    elif feature_type == "exons":
        result["Exons"] = gene_data.get("Exons", "")
        result["Position"] = gene_data.get("Position", "")
    elif feature_type == "organism": result["Source Organism"] = gene_data.get("Source Organism", "")
    return result
# --- END NCBI Data Fetching Functions ---

# Sidebar Search Options
st.sidebar.markdown("<h2 class='sub-header' style='font-size:1.8rem; margin-top:0px;'>Search Options</h2>", unsafe_allow_html=True)
database_option = st.sidebar.radio(
    "Filter by Annotation Status:",
    ["All Genes", "Fully Annotated Only", "Partially Annotated"]
)
organism_options = ["Homo sapiens", "Mus musculus", "Drosophila melanogaster", 
                   "Caenorhabditis elegans", "Saccharomyces cerevisiae"]
selected_organism = st.sidebar.selectbox("Select Organism:", organism_options)
search_mode = st.sidebar.radio(
    "Select Search Mode:",
    ["Single Gene", "Multiple Genes (Batch)"]
)
# Sidebar About Annotrax Quick Info
with st.sidebar.expander("ℹ️ About Annotrax (Quick Info)", expanded=True):
    st.markdown("""
    **Annotrax** simplifies gene annotation exploration.
    - Search single or multiple genes.
    - View function, location, exons.
    - Filter by organism & annotation status.
    - Download annotations (CSV, JSON, TXT) & FASTA sequences.
    
    Full details in the "About Annotrax" tab.
    """)

# Tabs
tab_home, tab_about, tab_dev_desk = st.tabs([
    "🏠 Home", 
    "ℹ️ About Annotrax", 
    "💡 Developer's Desk"
])

with tab_home:
    # --- HOME TAB CONTENT (Keep as is) ---
    if search_mode == "Single Gene":
        gene_name_input = st.text_input("🔍 Enter gene name (e.g., INS, BRCA1, TP53):", key="single_gene_input_home") 
        if gene_name_input:
            with st.spinner(f"Searching for {gene_name_input} in {selected_organism}..."):
                gene_annotations = fetch_gene_annotation(gene_name_input, selected_organism)
            if gene_annotations:
                if database_option == "Fully Annotated Only":
                    gene_annotations = [g for g in gene_annotations if g["Annotation Status"] == "Complete"]
                elif database_option == "Partially Annotated":
                    gene_annotations = [g for g in gene_annotations if g["Annotation Status"] == "Partial"]
                if not gene_annotations:
                    st.warning(f"⚠️ No genes matching '{database_option}' criteria were found for '{gene_name_input}'.")
                else:
                    st.success(f"✅ Found {len(gene_annotations)} results for {gene_name_input}")
                    if len(gene_annotations) > 1:
                        tabs_results = st.tabs([f"{g['Gene Symbol']} ({g['Gene ID']})" for g in gene_annotations])
                        for i, tab_result in enumerate(tabs_results):
                            with tab_result:
                                gene_data_tab = gene_annotations[i]
                                st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                                col1, col2 = st.columns([3, 1])
                                with col1:
                                    st.markdown(f"#### {gene_data_tab['Gene Symbol']} - {gene_data_tab['Official Name']}")
                                    st.markdown(f"**Aliases:** {gene_data_tab['Aliases']}")
                                    st.markdown(f"**Location:** Chromosome {gene_data_tab['Chromosome']}, Position {gene_data_tab['Position']}")
                                    st.markdown(f"**Exons:** {gene_data_tab['Exons']}")
                                    st.markdown(f"**Function:**")
                                    st.caption(f"{gene_data_tab['Function']}")
                                with col2:
                                    st.markdown(f"**Organism:** {gene_data_tab['Source Organism']}")
                                    st.markdown(f"**Gene ID:** {gene_data_tab['Gene ID']}")
                                    st.markdown(f"**Annotation:** {gene_data_tab['Annotation Status']}")
                                    if gene_data_tab['Nucleotide IDs'] != "N/A":
                                        first_nuccore_id = gene_data_tab['Nucleotide IDs'].split(", ")[0]
                                        gene_id_key_part = gene_data_tab['Gene ID']
                                        if st.button(f"View Sequence", key=f"seq_view_tab_{gene_id_key_part}_{i}", help="View GenBank details and sequence"):
                                            with st.spinner("Fetching GenBank details..."): seq_data = fetch_gene_sequence(first_nuccore_id)
                                            if seq_data:
                                                st.markdown("##### Sequence Information (GenBank)")
                                                st.markdown(f"**ID:** {seq_data['id']} | **Length:** {seq_data['length']} bp")
                                                with st.expander("View Sequence (First 1000 bp)"): st.code(seq_data['sequence'], language=None)
                                                with st.expander("View Features (first 10)"):
                                                    for feat in seq_data['features']:
                                                        st.markdown(f"**{feat['type']}** at {feat['location']}")
                                                        for k, v_list in feat['qualifiers'].items():
                                                            if k in ['product', 'gene', 'note', 'protein_id']: st.markdown(f"- {k}: {v_list}")
                                            else: st.warning("Could not fetch GenBank sequence data.")
                                        if st.button(f"Download FASTA", key=f"fasta_dl_tab_{gene_id_key_part}_{i}"):
                                            with st.spinner(f"Fetching FASTA for {first_nuccore_id}..."): fasta_content = fetch_fasta_sequence_data(first_nuccore_id)
                                            if fasta_content:
                                                st.download_button(label="⬇️ FASTA Ready", data=fasta_content, file_name=f"{gene_data_tab['Gene Symbol']}_{first_nuccore_id}.fasta", mime="text/fasta", key=f"actual_fasta_dl_tab_{gene_id_key_part}_{i}")
                                            else: st.warning("Could not fetch FASTA sequence.")
                                st.markdown("</div>", unsafe_allow_html=True)
                    else: 
                        gene_data_single = gene_annotations[0]
                        st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                        col1, col2 = st.columns([3, 1])
                        with col1:
                            st.markdown(f"#### {gene_data_single['Gene Symbol']} - {gene_data_single['Official Name']}")
                            st.markdown(f"**Aliases:** {gene_data_single['Aliases']}")
                            st.markdown(f"**Location:** Chromosome {gene_data_single['Chromosome']}, Position {gene_data_single['Position']}")
                            st.markdown(f"**Exons:** {gene_data_single['Exons']}")
                            st.markdown(f"**Function:**")
                            st.caption(f"{gene_data_single['Function']}")
                        with col2:
                            st.markdown(f"**Organism:** {gene_data_single['Source Organism']}")
                            st.markdown(f"**Gene ID:** {gene_data_single['Gene ID']}")
                            st.markdown(f"**Annotation:** {gene_data_single['Annotation Status']}")
                            if gene_data_single['Nucleotide IDs'] != "N/A":
                                first_nuccore_id = gene_data_single['Nucleotide IDs'].split(", ")[0]
                                gene_id_key_part = gene_data_single['Gene ID']
                                if st.button(f"View Sequence", key=f"seq_view_single_{gene_id_key_part}", help="View GenBank details and sequence"):
                                    with st.spinner("Fetching GenBank details..."): seq_data = fetch_gene_sequence(first_nuccore_id)
                                    if seq_data:
                                        st.markdown("##### Sequence Information (GenBank)")
                                        st.markdown(f"**ID:** {seq_data['id']} | **Length:** {seq_data['length']} bp")
                                        with st.expander("View Sequence (First 1000 bp)"): st.code(seq_data['sequence'], language=None)
                                        with st.expander("View Features (first 10)"):
                                            for feat in seq_data['features']:
                                                st.markdown(f"**{feat['type']}** at {feat['location']}")
                                                for k, v_list in feat['qualifiers'].items():
                                                    if k in ['product', 'gene', 'note', 'protein_id']: st.markdown(f"- {k}: {v_list}")
                                    else: st.warning("Could not fetch GenBank sequence data.")
                                if st.button(f"Download FASTA", key=f"fasta_dl_single_{gene_id_key_part}"):
                                    with st.spinner(f"Fetching FASTA for {first_nuccore_id}..."): fasta_content = fetch_fasta_sequence_data(first_nuccore_id)
                                    if fasta_content:
                                        st.download_button(label="⬇️ FASTA Ready",data=fasta_content,file_name=f"{gene_data_single['Gene Symbol']}_{first_nuccore_id}.fasta",mime="text/fasta",key=f"actual_fasta_dl_single_{gene_id_key_part}")
                                    else: st.warning("Could not fetch FASTA sequence.")
                        st.markdown("</div>", unsafe_allow_html=True)
                    
                    st.markdown("<h3 class='sub-header' style='font-size:1.4rem; margin-top: 20px; border-bottom: none;'>📥 Download Annotation Options</h3>", unsafe_allow_html=True) 
                    genes_to_consider_for_download = []
                    if len(gene_annotations) > 1:
                        st.markdown("**Select genes for annotation download:**")
                        for i, gene_item in enumerate(gene_annotations):
                            if st.checkbox(f"Include {gene_item['Gene Symbol']} ({gene_item['Gene ID']})", value=True, key=f"dl_sel_anno_home_{gene_item['Gene ID']}"):
                                genes_to_consider_for_download.append(gene_item)
                        if not genes_to_consider_for_download: st.info("No genes selected for annotation download. Please check at least one.")
                    else: genes_to_consider_for_download = gene_annotations

                    if genes_to_consider_for_download: 
                        dl_col1, dl_col2 = st.columns(2)
                        with dl_col1:
                            st.markdown("**Select fields for download:**")
                            download_all_fields = st.checkbox("All annotation fields", value=True, key="anno_dl_all_fields_home")
                            download_function_field, download_exons_field, download_organism_field = False, False, False
                            if not download_all_fields:
                                download_function_field = st.checkbox("Gene functions", value=True, key="anno_dl_func_home")
                                download_exons_field = st.checkbox("Exon information", value=False, key="anno_dl_exons_home")
                                download_organism_field = st.checkbox("Source organism", value=False, key="anno_dl_org_home")
                        with dl_col2:
                            st.markdown("**Download format:**")
                            download_format_anno = st.radio("Download Format", ["CSV", "JSON", "TXT"], key="anno_dl_format_home", label_visibility="collapsed", horizontal=True)
                        
                        if st.button("Download Selected Annotations", key="dl_anno_button_home", type="primary"):
                            download_data_anno = []
                            for gene_to_dl in genes_to_consider_for_download:
                                if download_all_fields: download_data_anno.append(gene_to_dl)
                                else:
                                    gene_export = {"Gene Symbol": gene_to_dl["Gene Symbol"], "Gene ID": gene_to_dl["Gene ID"]} 
                                    if download_function_field: gene_export.update(extract_feature_data(gene_to_dl, "function"))
                                    if download_exons_field: gene_export.update(extract_feature_data(gene_to_dl, "exons"))
                                    if download_organism_field: gene_export.update(extract_feature_data(gene_to_dl, "organism"))
                                    download_data_anno.append(gene_export)
                            if not download_data_anno: st.warning("No data prepared for download. Check selections.")
                            else:
                                df_anno = pd.DataFrame(download_data_anno)
                                dl_filename_base = gene_name_input if len(genes_to_consider_for_download) == 1 else "selected_genes"
                                file_suffix = f"{dl_filename_base}_annotations"
                                if download_format_anno == "CSV": st.download_button(label="⬇️ Download CSV", data=df_anno.to_csv(index=False).encode("utf-8"), file_name=f"{file_suffix}.csv", mime="text/csv", key="csv_dl_anno_home")
                                elif download_format_anno == "JSON": st.download_button(label="⬇️ Download JSON", data=df_anno.to_json(orient="records", indent=2), file_name=f"{file_suffix}.json", mime="application/json", key="json_dl_anno_home")
                                else: st.download_button(label="⬇️ Download TXT", data=df_anno.to_string(index=False), file_name=f"{file_suffix}.txt", mime="text/plain", key="txt_dl_anno_home")
            else: st.warning(f"⚠️ No gene data found for '{gene_name_input}'. Please check spelling or try another gene.")
    else: # Batch search mode
        st.markdown("### Batch Gene Search")
        st.markdown("Enter multiple gene names separated by commas, spaces, or new lines.")
        batch_input_genes = st.text_area("Enter gene names:", height=100, key="batch_gene_input_home", placeholder="e.g. TP53, BRCA1, EGFR\nINS\nMYC, PTEN")
        if batch_input_genes:
            gene_names_batch = list(set([item.strip() for item in batch_input_genes.replace(",", " ").replace("\n", " ").split() if item.strip()]))
            if gene_names_batch:
                st.markdown(f"Found {len(gene_names_batch)} unique gene names to search.")
                if st.button("Search Batch", key="search_batch_button_home", type="primary"):
                    all_results_batch = []
                    progress_bar = st.progress(0.0)
                    status_text = st.empty()
                    total_genes = len(gene_names_batch)
                    for i, gene_name_item in enumerate(gene_names_batch):
                        status_text.text(f"Searching for {gene_name_item}... ({i+1}/{total_genes})")
                        time.sleep(0.35) 
                        gene_annotations_item = fetch_gene_annotation(gene_name_item, selected_organism)
                        if gene_annotations_item:
                            filtered_annotations = gene_annotations_item
                            if database_option == "Fully Annotated Only": filtered_annotations = [g for g in gene_annotations_item if g["Annotation Status"] == "Complete"]
                            elif database_option == "Partially Annotated": filtered_annotations = [g for g in gene_annotations_item if g["Annotation Status"] == "Partial"]
                            all_results_batch.extend(filtered_annotations)
                        progress_bar.progress(float(i + 1) / total_genes)
                    
                    if all_results_batch:
                        status_text.success(f"Batch search complete! Found {len(all_results_batch)} gene records matching your criteria.")
                        df_display_batch = pd.DataFrame([{"Symbol": g.get("Gene Symbol", "N/A"), "Name": g.get("Official Name", "N/A"), "Organism": g.get("Source Organism", "N/A"), "Chr.": g.get("Chromosome", "N/A"), "Exons": g.get("Exons", "N/A"), "Status": g.get("Annotation Status", "N/A")} for g in all_results_batch])
                        st.dataframe(df_display_batch, height=300)
                        st.markdown("<h3 class='sub-header' style='font-size:1.4rem; margin-top: 20px; border-bottom: none;'>📥 Batch Annotation Download Options</h3>", unsafe_allow_html=True)
                        b_col1, b_col2 = st.columns(2)
                        with b_col1:
                            st.markdown("**Select fields for download:**")
                            batch_dl_all = st.checkbox("All fields", value=True, key="batch_dl_all_fields_home")
                            batch_dl_func, batch_dl_exons, batch_dl_org = False, False, False
                            if not batch_dl_all:
                                batch_dl_func = st.checkbox("Gene functions", value=True, key="batch_dl_func_home")
                                batch_dl_exons = st.checkbox("Exon information", value=False, key="batch_dl_exons_home")
                                batch_dl_org = st.checkbox("Source organism", value=False, key="batch_dl_org_home")
                        with b_col2:
                            st.markdown("**Download format:**")
                            batch_dl_format = st.radio("Download Format", ["CSV", "JSON", "TXT"], key="batch_dl_format_select_home", label_visibility="collapsed", horizontal=True)
                        
                        if st.button("Download Batch Annotations", key="dl_batch_anno_button_home", type="primary"):
                            batch_dl_data = []
                            for gene_res in all_results_batch:
                                if batch_dl_all: batch_dl_data.append(gene_res)
                                else:
                                    gene_export_batch = {"Gene Symbol": gene_res.get("Gene Symbol", "N/A"), "Gene ID": gene_res.get("Gene ID", "N/A")}
                                    if batch_dl_func: gene_export_batch.update(extract_feature_data(gene_res, "function"))
                                    if batch_dl_exons: gene_export_batch.update(extract_feature_data(gene_res, "exons"))
                                    if batch_dl_org: gene_export_batch.update(extract_feature_data(gene_res, "organism"))
                                    batch_dl_data.append(gene_export_batch)
                            if not batch_dl_data: st.warning("No data prepared for batch download.")
                            else:
                                batch_df_dl = pd.DataFrame(batch_dl_data)
                                file_suffix_batch = "batch_gene_annotations"
                                if batch_dl_format == "CSV": st.download_button(label="⬇️ Download CSV", data=batch_df_dl.to_csv(index=False).encode("utf-8"), file_name=f"{file_suffix_batch}.csv", mime="text/csv", key="csv_dl_batch_home")
                                elif batch_dl_format == "JSON": st.download_button(label="⬇️ Download JSON", data=batch_df_dl.to_json(orient="records", indent=2), file_name=f"{file_suffix_batch}.json", mime="application/json", key="json_dl_batch_home")
                                else: st.download_button(label="⬇️ Download TXT", data=batch_df_dl.to_string(index=False), file_name=f"{file_suffix_batch}.txt", mime="text/plain", key="txt_dl_batch_home")
                    else: status_text.warning("⚠️ No genes found matching your criteria in the batch search.")
            else: st.info("Enter gene names to begin batch search.")
    # --- END HOME TAB CONTENT ---

with tab_about:
    st.markdown("<h2 class='sub-header'>About Annotrax</h2>", unsafe_allow_html=True)
    st.markdown("<div class='about-content-wrapper'>", unsafe_allow_html=True)

    st.markdown("<h3 class='about-section-header'>What is Annotrax & Its Purpose?</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p>Annotrax is a user-friendly web application designed to streamline the process of exploring, retrieving, and analyzing gene annotations. 
    It serves as an efficient bridge to the vast genomic data available in the National Center for Biotechnology Information (NCBI) databases, 
    presenting it in an accessible and actionable format.</p>
    <p>The primary purpose of Annotrax is to empower students, researchers, and bioinformatics
    enthusiasts by saving them valuable time and effort in their quest for genomic insights. We aim to make complex gene data more manageable, 
    thereby accelerating research, facilitating learning, and contributing to the broader understanding of genomics.</p>
    """, unsafe_allow_html=True)

    st.markdown("<h3 class='about-section-header'>Key Features</h3>", unsafe_allow_html=True)
    st.markdown("""
    - **Comprehensive Gene Search:** Search for single genes or perform batch searches for multiple genes simultaneously using standard gene symbols.
    - **Detailed Gene Information:** View critical details such as official gene names, aliases, chromosomal location, exon counts, and functional summaries.
    - **Organism-Specific Filtering:** Easily filter results by common model organisms like *Homo sapiens*, *Mus musculus*, *Drosophila melanogaster*, and more.
    - **Annotation Status Filter:** Narrow down searches to "Fully Annotated" or "Partially Annotated" genes based on the completeness of their NCBI summaries.
    - **Sequence Data Access:** View GenBank sequence details (ID, length, features) and download gene sequences in FASTA format for further analysis.
    - **Flexible Data Export:** Download curated annotation data in user-friendly formats (CSV, JSON, TXT) with options to select specific fields for targeted analysis.
    - **User-Friendly Interface:** Intuitive navigation and clear presentation of complex data, designed for both novice and experienced users in the field of bioinformatics.
    - **Efficient & Responsive:** Leverages caching and NCBI's Entrez API for up-to-date and relatively fast data retrieval, while respecting API usage guidelines.
    """) 

    with st.expander("📖 User Manual"):
        st.markdown("""
        This guide will help you navigate and utilize the features of Annotrax effectively.

        #### Getting Started
        1.  **NCBI Email**: For NCBI Entrez API use, an email is set by the developer. This ensures compliance with NCBI policies.
        2.  **Sidebar Options (Left Panel)**:
            *   **Filter by Annotation Status**: Filter genes based on their annotation status: "All Genes", "Fully Annotated Only", or "Partially Annotated".
            *   **Select Organism**: Choose the organism of interest from the dropdown list (e.g., Homo sapiens, Mus musculus).
            *   **Select Search Mode**:
                *   **Single Gene**: For searching one gene at a time.
                *   **Multiple Genes (Batch)**: For searching a list of genes.

        #### Single Gene Search (in Home Tab)
        1.  Select "Single Gene" in the sidebar search mode.
        2.  Enter the gene name (e.g., TP53, BRCA1) in the search box on the Home tab.
        3.  Results will display gene details: official name, aliases, function, location, exon count, etc.
        4.  If multiple gene records match your query (e.g., due to different gene IDs for the same symbol or related genes), they will be displayed in separate sub-tabs for clarity.
        5.  **View Sequence Data**: For each gene with associated Nucleotide IDs, click the "View Sequence" button. This will fetch and display GenBank details including Accession ID, sequence length, a preview of the sequence (first 1000 bp), and key features like CDS, gene, exon, mRNA.
        6.  **Download FASTA**: Click the "Download FASTA" button to retrieve the gene sequence in FASTA format for the primary Nucleotide ID linked to the gene.
        7.  **Download Annotations**:
            *   If multiple gene results are displayed, you can select which specific genes to include in your annotation download.
            *   Choose to download "All annotation fields" or select specific fields like "Gene functions," "Exon information," or "Source organism."
            *   Select your preferred download format: CSV, JSON, or TXT.
            *   Click the "Download Selected Annotations" button.

        #### Batch Gene Search (in Home Tab)
        1.  Select "Multiple Genes (Batch)" in the sidebar search mode.
        2.  In the text area on the Home tab, enter multiple gene names. You can separate them by commas, spaces, or new lines.
        3.  Click "Search Batch." A progress bar will indicate the search status as Annotrax queries NCBI for each gene.
        4.  A table summarizing the results for all found genes will be displayed.
        5.  **Download Batch Annotations**:
            *   Similar to single gene downloads, choose to download "All fields" or select specific annotation fields.
            *   Select your desired download format (CSV, JSON, TXT).
            *   Click the "Download Batch Annotations" button to save the compiled data.

        #### Tips for Effective Use
        *   Use standard gene symbols (e.g., HUGO nomenclature for human genes) for the most accurate results.
        *   Be mindful of NCBI API usage limits. The app includes a small delay between requests in batch mode to comply with these limits.
        *   Data is cached for 1 hour. This means if you search for the same gene in the same organism within an hour, the results will load much faster and reduce calls to the NCBI API.
        """) 

    st.markdown("<h3 class='about-section-header'>Future Aspects</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p>Annotrax is an evolving project. Potential future enhancements include:</p>
    """, unsafe_allow_html=True)
    st.markdown("""
    - **Protein Information:** Integration of protein-specific data from UniProt or NCBI Protein.
    - **Expanded Organism List:** Option to search for genes in a wider range of organisms.
    - **Pathway & Interaction Data:** Links to relevant pathway databases (e.g., KEGG, Reactome) or protein-protein interaction networks.
    - **Advanced Visualization:** Graphical representation of gene structures or annotation distributions.
    - **User Accounts & Saved Searches:** Functionality for users to save their common searches or preferred settings.
    - **Direct Variant Information:** Links or integration with variant databases like dbSNP or ClinVar for human genes.
    """)

    st.markdown("<h3 class='about-section-header'>Acknowledgements</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p>This project, Annotrax, stands on the shoulders of giants and is fueled by the spirit of collaborative learning and mentorship. 
    I extend my deepest gratitude to all who have contributed, directly or indirectly, to its development.</p>
    <p>A very special and heartfelt acknowledgement goes to 
    <span class='acknowledgement-highlight'>Dr. Kushagra Kashyap</span>. 
    His invaluable guidance, profound expertise in bioinformatics, and unwavering support were instrumental 
    throughout the conception and development of this application. Dr. Kashyap's insightful feedback, 
    encouragement to explore innovative solutions, and dedication to fostering a rich learning environment 
    have not only shaped Annotrax but have also significantly contributed to my growth as a developer and researcher. 
    His vision for leveraging computational tools to advance biological research served as a constant source of inspiration. 
    Thank you, Dr. Kashyap, for your exceptional mentorship.</p>
    <p>I am also grateful to:</p>
    <ul>
        <li>The creators and maintainers of the <a href="https://www.ncbi.nlm.nih.gov/" target="_blank">National Center for Biotechnology Information (NCBI)</a> and its Entrez API, which provides the foundational data for Annotrax.</li>
        <li>The developers of <a href="https://biopython.org/" target="_blank">Biopython</a>, for their powerful and versatile bioinformatics library.</li>
        <li>The <a href="https://streamlit.io/" target="_blank">Streamlit</a> team, for creating an amazing framework that makes building data applications so intuitive.</li>
        <li>The open-source community, for the countless tools and resources that make projects like this possible.</li>
    </ul>
    """, unsafe_allow_html=True)


    st.markdown("<h3 class='about-section-header'>Data & Technology</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p>
        Annotrax retrieves all its genomic data in real-time from the National Center for Biotechnology Information (NCBI) 
        using their Entrez Programming Utilities (E-utilities) API. Key Python libraries enabling this application include 
        <a href="https://biopython.org/" target="_blank">Biopython</a> for interacting with NCBI and parsing biological data formats, 
        and <a href="https://pandas.pydata.org/" target="_blank">Pandas</a> for data manipulation and table generation. 
        The user interface is built with <a href="https://streamlit.io/" target="_blank">Streamlit</a>, 
        allowing for rapid development of interactive data applications.
    </p>
    """, unsafe_allow_html=True)

    st.markdown("<h3 class='about-section-header'>Contact & Support</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p>If you have any questions, encounter issues, or have suggestions for improving Annotrax, please feel free to reach out. 
    Your feedback is valuable!</p>
    <p><strong>Email:</strong> <a href="mailto:kaletejal05@mail.com">kaletejal05@mail.com</a></p>
    <p><strong>LinkedIn:</strong> <a href="https://www.linkedin.com/in/tejal-kale-0658bb209/" target="_blank">Connect with Tejal Kale on LinkedIn</a></p>
    <p>We will do our best to respond to your queries as soon as possible.</p>
    """, unsafe_allow_html=True)

    st.markdown("</div>", unsafe_allow_html=True) 

with tab_dev_desk:
    st.markdown("<h2 class='sub-header'>Developer's Desk</h2>", unsafe_allow_html=True)
    st.markdown("<div class='developer-desk-container'>", unsafe_allow_html=True)
    col_img, col_text = st.columns([1, 2.5]) 
    with col_img:
        image_url = "https://media.licdn.com/dms/image/v2/D5603AQFfI1KWVSWl0Q/profile-displayphoto-shrink_400_400/B56ZRRpRfTH0Ag-/0/1736536563970?e=1752105600&v=beta&t=6o_NtvfTuXoxDGHXmXaLYpbi3-blU9E8eHv-OO1gScU"
        
        st.markdown(
            f"""
            <div style="display: flex; justify-content: center; flex-direction: column; align-items: center;">
                <img src="{image_url}" alt="Tejal Kale" class="developer-profile-image" style="width: 200px;"> 
                <p class="developer-name-caption">Tejal Kale</p>
            </div>
            """,
            unsafe_allow_html=True
        )

    with col_text:
        st.markdown("<p class='developer-role'>M.Sc. Bioinformatics Student, DES Pune University</p>", unsafe_allow_html=True)
        st.markdown("""
        <div class='developer-intro-text'>
        As someone who's spent countless hours navigating the complexities of gene annotation, I know firsthand the frustration of juggling multiple tools and databases. 
        That's why I created Annotrax - a platform born out of my own struggles and passion for bioinformatics. 
        I wanted to build a tool that would save others the time and effort I've wasted, and instead, empower them to focus on what really matters: discovering new insights and advancing our understanding of the genetic code. 
        With Annotrax, I've aimed to create a seamless, intuitive experience that brings together the best of gene annotation and analysis in one place. 
        My hope is that it becomes an indispensable companion for researchers, students, and scientists, helping them unlock the secrets of the genome and drive innovation.
        </div>
        <p class="developer-signature">— Tejal Kale, Developer of Annotrax</p>
        """, unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div class="info-text" style="text-align: center; padding-top: 10px;">
<small>Annotrax retrieves data from NCBI's Entrez API. Please use responsibly and respect NCBI's usage guidelines. Cache active for 1 hour.</small><br>
<small>Developed by Tejal Kale.</small>
</div>
""", unsafe_allow_html=True)
