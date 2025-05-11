import streamlit as st
import pandas as pd
import numpy as np
import base64
from Bio import Entrez, SeqIO
import io
import json
import time
import os

# Set your email address for NCBI Entrez API
Entrez.email = "kaletejal05@mail.com"

# Page configuration
st.set_page_config(
    page_title="Annotrax",
    layout="wide"
)

# Custom CSS styles
st.markdown("""
<style>
    /* Main app background */
    .stApp {
        background-color: #1a1a1a;
    }

    /* Sidebar styling */
    [data-testid=stSidebar] {
        background-color: #262626 !important;
    }

    /* Text colors */
    * {
        color: #4CAF50 !important;
    }
    
    /* Headers */
    h1, h2, h3, h4, h5, h6 {
        color: #4CAF50 !important;
    }

    /* Containers */
    .main-container {
        background-color: #333333;
        padding: 2rem;
        border-radius: 10px;
        margin: 1rem 0;
    }

    /* Cards */
    .feature-card {
        background-color: #404040 !important;
        border-radius: 10px;
        padding: 1.5rem;
        margin: 1rem 0;
        border: 1px solid #4CAF50;
    }

    /* Result boxes */
    .result-box {
        background-color: #404040;
        border-radius: 5px;
        padding: 20px;
        margin-bottom: 20px;
        border: 1px solid #4CAF50;
    }

    /* Input fields */
    .stTextInput input {
        background-color: #ffffff !important;
        border: 1px solid #4CAF50 !important;
        color: #4CAF50 !important;
    }

    /* Buttons - MODIFIED */
    .stButton>button {
        background-color: #ffffff !important;  /* White background */
        color: #4CAF50 !important; /* Green text */
        border: 1px solid #4CAF50 !important; /* Green border */
    }
    
    .stButton>button:hover {
        background-color: #4CAF50 !important;  /* Green background on hover */
        color: white !important; /* White text on hover */
    }

    /* Dataframes */
    .stDataFrame {
        background-color: #333333 !important;
    }
</style>
""", unsafe_allow_html=True)

# Navigation tabs
TABS = ["🏠 Home", "🔍 Gene Search", "📚 Batch Search", "📖 About"]
st.sidebar.markdown("Annotrax")
selected_tab = st.sidebar.radio("Navigation", TABS)

# Cache functions
@st.cache_data(ttl=3600)
def fetch_gene_annotation(gene_name, organism="Homo sapiens"):
    try:
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            return None
        
        results = []
        for gene_id in record["IdList"]:
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(summary_handle)
            summary_handle.close()
            
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            genomic_info = docsum.get("GenomicInfo", [{}])
            
            # Get nucleotide IDs
            nuccore_ids = []
            try:
                link_handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id)
                link_results = Entrez.read(link_handle)
                link_handle.close()
                if link_results[0]["LinkSetDb"]:
                    for link_info in link_results[0]["LinkSetDb"]:
                        if link_info["DbTo"] == "nuccore":
                            nuccore_ids.extend([link["Id"] for link in link_info["Link"]])
                            break
            except Exception:
                pass

            gene_data = {
                "Gene ID": gene_id,
                "Gene Symbol": docsum.get("Name", gene_name),
                "Official Name": docsum.get("Description", "N/A"),
                "Aliases": ", ".join(docsum.get("OtherAliases", "").split(", ")[:5]),
                "Function": docsum.get("Summary", "No function description available"),
                "Chromosome": genomic_info[0].get("ChrLoc", "N/A") if genomic_info else "N/A",
                "Position": f"{genomic_info[0].get('ChrStart', 'N/A')}-{genomic_info[0].get('ChrStop', 'N/A')}" if genomic_info else "N/A",
                "Exons": genomic_info[0].get("ExonCount", "N/A") if genomic_info else "N/A",
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
                    "type": feature.type,
                    "location": str(feature.location),
                    "qualifiers": {k: str(v[0]) if isinstance(v, list) and v else str(v) 
                                  for k, v in feature.qualifiers.items()}
                }
                features.append(feature_info)
        
        return {
            "id": record.id,
            "name": record.name,
            "description": record.description,
            "length": len(record.seq),
            "features": features[:10], # Limiting features to avoid too much data
            "sequence": str(record.seq)[:1000] + ("..." if len(record.seq) > 1000 else "") # Sequence preview
        }
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

# Home Tab
if selected_tab == "🏠 Home":
    st.markdown("<div class='main-container'>", unsafe_allow_html=True)
    
    # Main header
    st.markdown("<h1 class='main-header'> Annotrax</h1>", unsafe_allow_html=True)
    
    # Title box
    st.markdown("""
    <div class='title-box' style='background-color: #2a2a2a; padding: 1.5rem; border-radius: 8px; margin-bottom: 1.5rem; border-left: 5px solid #4CAF50;'>
        <h3 style='color: #4CAF50!important; margin-bottom: 0.5rem;'>Your Comprehensive Gene Annotation Platform</h3>
        <p style='color: #e0e0e0;'>Annotrax provides researchers and students with fast, reliable access to gene annotations 
        and sequence data from trusted biological databases.</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    with col1:
        st.markdown("""
        <div class='feature-card'>
            <h4 style='color: #4CAF50!important;'>🌟 Key Features</h4>
            <ul>
                <li>Instant gene annotation retrieval</li>
                <li>Batch processing for multiple genes</li>
                <li>Multiple download formats (CSV, JSON, FASTA)</li>
                <li>Comprehensive sequence analysis tools</li>
                <li>Cross-species compatibility</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class='feature-card'>
            <h4 style='color: #4CAF50!important;'>📈 Quick Start</h4>
            <ol>
                <li>Select 'Gene Search' tab</li>
                <li>Enter gene name (e.g., TP53)</li>
                <li>Select organism</li>
                <li>Click 'Search' & Explore results!</li>
            </ol>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class='feature-card'>
            <h4 style='color: #4CAF50!important;'>📊 Database Power</h4>
            <ul>
                <li>Access to NCBI Entrez</li>
                <li>Supports major model organisms</li>
                <li>Data updated as per NCBI</li>
                <li>High query success rate</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("</div>", unsafe_allow_html=True)

# Gene Search Tab
elif selected_tab == "🔍 Gene Search":
    st.markdown("<h1 class='main-header'>Gene Annotation Search</h1>", unsafe_allow_html=True)
    
    with st.sidebar.expander("🔎 Search Filters", expanded=True):
        database_option = st.radio(
            "Annotation Status Filter",
            ["All", "Fully Annotated Only", "Partially Annotated"],
            key="db_filter_gene_search"
        )
        organism_options = ["Homo sapiens", "Mus musculus", "Drosophila melanogaster", 
                          "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Danio rerio", "Arabidopsis thaliana"]
        selected_organism = st.selectbox("Select Organism", organism_options, key="organism_gene_search")

    gene_name_input = st.text_input("🔍 Enter gene name (e.g., INS, BRCA1, TP53):", key="gene_input_main")
    
    # Initialize session state for annotations if not present
    if 'gene_annotations' not in st.session_state:
        st.session_state.gene_annotations = None
    if 'current_gene_search_term' not in st.session_state:
        st.session_state.current_gene_search_term = ""
    if 'current_organism_search' not in st.session_state:
        st.session_state.current_organism_search = ""

    if st.button("Search", key="search_button_gene_search", use_container_width=True, type="primary"):
        if gene_name_input:
            with st.spinner(f"Searching for '{gene_name_input}' in '{selected_organism}'..."):
                st.session_state.gene_annotations = fetch_gene_annotation(gene_name_input, selected_organism)
                st.session_state.current_gene_search_term = gene_name_input
                st.session_state.current_organism_search = selected_organism
        else:
            st.warning("Please enter a gene name to search.")
            st.session_state.gene_annotations = None # Clear previous results if search is empty
            
    if st.session_state.gene_annotations is not None:
        gene_annotations_to_display = st.session_state.gene_annotations
        
        if database_option == "Fully Annotated Only":
            gene_annotations_to_display = [g for g in gene_annotations_to_display if g["Annotation Status"] == "Complete"]
        elif database_option == "Partially Annotated":
            gene_annotations_to_display = [g for g in gene_annotations_to_display if g["Annotation Status"] == "Partial"]
            
        if not gene_annotations_to_display:
            if st.session_state.current_gene_search_term: # Only show warning if a search was made
                st.warning(f"⚠️ No genes matching '{database_option}' criteria were found for '{st.session_state.current_gene_search_term}' in '{st.session_state.current_organism_search}'.")
        else:
            st.success(f"✅ Found {len(gene_annotations_to_display)} result(s) for '{st.session_state.current_gene_search_term}' (Filter: {database_option})")
            
            # Determine if tabs are needed
            if len(gene_annotations_to_display) > 1:
                tab_titles = [f"{g['Gene Symbol']} ({g['Gene ID']})" for g in gene_annotations_to_display]
                tabs = st.tabs(tab_titles)
                
                for i, tab_content in enumerate(tabs):
                    with tab_content:
                        gene_data = gene_annotations_to_display[i]
                        st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                        
                        col1, col2 = st.columns([3, 1])
                        with col1:
                            st.markdown(f"### {gene_data['Gene Symbol']} - {gene_data['Official Name']}")
                            st.markdown(f"**Aliases:** {gene_data['Aliases']}")
                            st.markdown(f"**Location:** Chromosome {gene_data['Chromosome']}, Position {gene_data['Position']}")
                            st.markdown(f"**Exons:** {gene_data['Exons']}")
                            st.markdown(f"**Function:**")
                            st.markdown(f"{gene_data['Function']}")
                        
                        with col2:
                            st.markdown(f"**Organism:** {gene_data['Source Organism']}")
                            st.markdown(f"**Gene ID:** {gene_data['Gene ID']}")
                            st.markdown(f"**Annotation:** {gene_data['Annotation Status']}")
                            st.markdown(f"**Nucleotide IDs:** {gene_data['Nucleotide IDs']}")

                        st.markdown("---")
                        
                        seq_col1, seq_col2, seq_col3 = st.columns([1, 1, 2])
                        
                        nuccore_id_to_fetch = gene_data['Nucleotide IDs'].split(", ")[0] if gene_data['Nucleotide IDs'] != "N/A" else None

                        with seq_col1:
                            if nuccore_id_to_fetch:
                                if st.button("View Sequence", key=f"view_seq_{gene_data['Gene ID']}", help="Fetch and display GenBank sequence details"):
                                    with st.spinner("Fetching sequence details..."):
                                        st.session_state[f"seq_details_{gene_data['Gene ID']}"] = fetch_gene_sequence(nuccore_id_to_fetch)
                            else:
                                st.caption("No Nucleotide ID for sequence view.")
                        
                        with seq_col2:
                            if nuccore_id_to_fetch:
                                if st.button("Download FASTA", key=f"dl_fasta_btn_{gene_data['Gene ID']}", help="Download gene sequence in FASTA format"):
                                    with st.spinner("Generating FASTA..."):
                                        st.session_state[f"fasta_data_{gene_data['Gene ID']}"] = fetch_fasta_sequence_data(nuccore_id_to_fetch)
                            else:
                                st.caption("No Nucleotide ID for FASTA.")

                        if f"seq_details_{gene_data['Gene ID']}" in st.session_state and st.session_state[f"seq_details_{gene_data['Gene ID']}"]:
                            seq_data = st.session_state[f"seq_details_{gene_data['Gene ID']}"]
                            with st.expander("Sequence Details", expanded=False): # Default to collapsed
                                st.markdown(f"**NCBI ID:** {seq_data['id']}")
                                st.markdown(f"**Description:** {seq_data['description']}")
                                st.markdown(f"**Length:** {seq_data['length']} bp")
                                st.text_area("Sequence Preview (first 1000bp)", 
                                           value=seq_data['sequence'], 
                                           height=150, 
                                           disabled=True,
                                           key=f"text_seq_preview_{gene_data['Gene ID']}")
                                if seq_data['features']:
                                    st.markdown("**Features (first 10):**")
                                    for feat in seq_data['features']:
                                        st.caption(f"- **Type:** {feat['type']}, **Location:** {feat['location']}")


                        if f"fasta_data_{gene_data['Gene ID']}" in st.session_state and st.session_state[f"fasta_data_{gene_data['Gene ID']}"]:
                            st.download_button(
                                label="⬇️ FASTA Ready",
                                data=st.session_state[f"fasta_data_{gene_data['Gene ID']}"],
                                file_name=f"{gene_data['Gene Symbol']}_{nuccore_id_to_fetch}.fasta",
                                mime="text/plain",
                                key=f"dl_fasta_final_{gene_data['Gene ID']}",
                                help="Click to download the FASTA file"
                            )
                        
                        with seq_col3:
                            st.markdown("**Download Annotation Data**")
                            dl_fields = st.multiselect(
                                "Select fields for download:",
                                options=list(gene_data.keys()),
                                default=["Gene Symbol", "Official Name", "Function", "Chromosome", "Position", "Exons", "Source Organism", "Gene ID"],
                                key=f"fields_dl_{gene_data['Gene ID']}"
                            )
                            dl_format = st.radio("Select download format:", ["CSV", "JSON", "TXT"], 
                                               horizontal=True, 
                                               key=f"format_dl_{gene_data['Gene ID']}")
                            
                            if st.button("Prepare Download", key=f"prep_dl_{gene_data['Gene ID']}", type="secondary"):
                                filtered_data_dl = {k: gene_data[k] for k in dl_fields if k in gene_data}
                                if dl_format == "CSV":
                                    csv_data = pd.DataFrame([filtered_data_dl]).to_csv(index=False)
                                    st.session_state[f"dl_data_{gene_data['Gene ID']}"] = csv_data
                                    st.session_state[f"dl_fname_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.csv"
                                    st.session_state[f"dl_mime_{gene_data['Gene ID']}"] = "text/csv"
                                elif dl_format == "JSON":
                                    json_data_dl = json.dumps(filtered_data_dl, indent=2)
                                    st.session_state[f"dl_data_{gene_data['Gene ID']}"] = json_data_dl
                                    st.session_state[f"dl_fname_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.json"
                                    st.session_state[f"dl_mime_{gene_data['Gene ID']}"] = "application/json"
                                else: # TXT
                                    txt_data = "\n".join([f"{k}: {v}" for k, v in filtered_data_dl.items()])
                                    st.session_state[f"dl_data_{gene_data['Gene ID']}"] = txt_data
                                    st.session_state[f"dl_fname_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.txt"
                                    st.session_state[f"dl_mime_{gene_data['Gene ID']}"] = "text/plain"
                                
                                st.success(f"{dl_format} file prepared. Click below to download.")

                            if f"dl_data_{gene_data['Gene ID']}" in st.session_state:
                                st.download_button(
                                    label=f"⬇️ Download {st.session_state[f'dl_fname_{gene_data['Gene ID']}']}",
                                    data=st.session_state[f"dl_data_{gene_data['Gene ID']}"],
                                    file_name=st.session_state[f"dl_fname_{gene_data['Gene ID']}"],
                                    mime=st.session_state[f"dl_mime_{gene_data['Gene ID']}"],
                                    key=f"actual_dl_btn_{gene_data['Gene ID']}"
                                )
                        st.markdown("</div>", unsafe_allow_html=True)
            
            elif len(gene_annotations_to_display) == 1: # Single result, no tabs
                gene_data = gene_annotations_to_display[0]
                st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.markdown(f"### {gene_data['Gene Symbol']} - {gene_data['Official Name']}")
                    st.markdown(f"**Aliases:** {gene_data['Aliases']}")
                    st.markdown(f"**Location:** Chromosome {gene_data['Chromosome']}, Position {gene_data['Position']}")
                    st.markdown(f"**Exons:** {gene_data['Exons']}")
                    st.markdown(f"**Function:**")
                    st.markdown(f"{gene_data['Function']}")
                
                with col2:
                    st.markdown(f"**Organism:** {gene_data['Source Organism']}")
                    st.markdown(f"**Gene ID:** {gene_data['Gene ID']}")
                    st.markdown(f"**Annotation:** {gene_data['Annotation Status']}")
                    st.markdown(f"**Nucleotide IDs:** {gene_data['Nucleotide IDs']}")

                st.markdown("---")
                seq_col1, seq_col2, seq_col3 = st.columns([1, 1, 2])
                
                nuccore_id_to_fetch = gene_data['Nucleotide IDs'].split(", ")[0] if gene_data['Nucleotide IDs'] != "N/A" else None

                with seq_col1:
                    if nuccore_id_to_fetch:
                        if st.button("View Sequence", key=f"view_seq_single_{gene_data['Gene ID']}", help="Fetch and display GenBank sequence details"):
                            with st.spinner("Fetching sequence details..."):
                                st.session_state[f"seq_details_single_{gene_data['Gene ID']}"] = fetch_gene_sequence(nuccore_id_to_fetch)
                    else:
                        st.caption("No Nucleotide ID for sequence view.")
                
                with seq_col2:
                    if nuccore_id_to_fetch:
                        if st.button("Download FASTA", key=f"dl_fasta_btn_single_{gene_data['Gene ID']}", help="Download gene sequence in FASTA format"):
                            with st.spinner("Generating FASTA..."):
                                st.session_state[f"fasta_data_single_{gene_data['Gene ID']}"] = fetch_fasta_sequence_data(nuccore_id_to_fetch)
                    else:
                        st.caption("No Nucleotide ID for FASTA.")

                if f"seq_details_single_{gene_data['Gene ID']}" in st.session_state and st.session_state[f"seq_details_single_{gene_data['Gene ID']}"]:
                    seq_data = st.session_state[f"seq_details_single_{gene_data['Gene ID']}"]
                    with st.expander("Sequence Details", expanded=False):
                        st.markdown(f"**NCBI ID:** {seq_data['id']}")
                        st.markdown(f"**Description:** {seq_data['description']}")
                        st.markdown(f"**Length:** {seq_data['length']} bp")
                        st.text_area("Sequence Preview (first 1000bp)", 
                                   value=seq_data['sequence'], 
                                   height=150, 
                                   disabled=True,
                                   key=f"text_seq_preview_single_{gene_data['Gene ID']}")
                        if seq_data['features']:
                            st.markdown("**Features (first 10):**")
                            for feat in seq_data['features']:
                                 st.caption(f"- **Type:** {feat['type']}, **Location:** {feat['location']}")

                if f"fasta_data_single_{gene_data['Gene ID']}" in st.session_state and st.session_state[f"fasta_data_single_{gene_data['Gene ID']}"]:
                    st.download_button(
                        label="⬇️ FASTA Ready",
                        data=st.session_state[f"fasta_data_single_{gene_data['Gene ID']}"],
                        file_name=f"{gene_data['Gene Symbol']}_{nuccore_id_to_fetch}.fasta",
                        mime="text/plain",
                        key=f"dl_fasta_final_single_{gene_data['Gene ID']}",
                        help="Click to download the FASTA file"
                    )
                
                with seq_col3:
                    st.markdown("**Download Annotation Data**")
                    dl_fields_single = st.multiselect(
                        "Select fields for download:",
                        options=list(gene_data.keys()),
                        default=["Gene Symbol", "Official Name", "Function", "Chromosome", "Position", "Exons", "Source Organism", "Gene ID"],
                        key=f"fields_dl_single_{gene_data['Gene ID']}"
                    )
                    dl_format_single = st.radio("Select download format:", ["CSV", "JSON", "TXT"], 
                                       horizontal=True, 
                                       key=f"format_dl_single_{gene_data['Gene ID']}")
                    
                    if st.button("Prepare Download", key=f"prep_dl_single_{gene_data['Gene ID']}", type="secondary"):
                        filtered_data_dl_single = {k: gene_data[k] for k in dl_fields_single if k in gene_data}
                        if dl_format_single == "CSV":
                            csv_data_s = pd.DataFrame([filtered_data_dl_single]).to_csv(index=False)
                            st.session_state[f"dl_data_single_{gene_data['Gene ID']}"] = csv_data_s
                            st.session_state[f"dl_fname_single_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.csv"
                            st.session_state[f"dl_mime_single_{gene_data['Gene ID']}"] = "text/csv"
                        elif dl_format_single == "JSON":
                            json_data_dl_s = json.dumps(filtered_data_dl_single, indent=2)
                            st.session_state[f"dl_data_single_{gene_data['Gene ID']}"] = json_data_dl_s
                            st.session_state[f"dl_fname_single_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.json"
                            st.session_state[f"dl_mime_single_{gene_data['Gene ID']}"] = "application/json"
                        else: # TXT
                            txt_data_s = "\n".join([f"{k}: {v}" for k, v in filtered_data_dl_single.items()])
                            st.session_state[f"dl_data_single_{gene_data['Gene ID']}"] = txt_data_s
                            st.session_state[f"dl_fname_single_{gene_data['Gene ID']}"] = f"{gene_data['Gene Symbol']}_annotation.txt"
                            st.session_state[f"dl_mime_single_{gene_data['Gene ID']}"] = "text/plain"
                        
                        st.success(f"{dl_format_single} file prepared. Click below to download.")

                    if f"dl_data_single_{gene_data['Gene ID']}" in st.session_state:
                        st.download_button(
                            label=f"⬇️ Download {st.session_state[f'dl_fname_single_{gene_data['Gene ID']}']}",
                            data=st.session_state[f"dl_data_single_{gene_data['Gene ID']}"],
                            file_name=st.session_state[f"dl_fname_single_{gene_data['Gene ID']}"],
                            mime=st.session_state[f"dl_mime_single_{gene_data['Gene ID']}"],
                            key=f"actual_dl_btn_single_{gene_data['Gene ID']}"
                        )
                st.markdown("</div>", unsafe_allow_html=True)

    elif st.session_state.current_gene_search_term: # No results found and a search term exists
        st.info(f"No results found for '{st.session_state.current_gene_search_term}' in '{st.session_state.current_organism_search}'. Please check the gene name and organism or try a different query.")


# Batch Search Tab
elif selected_tab == "📚 Batch Search":
    st.markdown("<h1 class='main-header'>Batch Gene Analysis</h1>", unsafe_allow_html=True)
    
    st.markdown("### Upload gene list (TXT) or enter manually:")
    
    # File uploader
    uploaded_file = st.file_uploader("Upload a .txt file with one gene name per line", type="txt", key="batch_file_uploader")
    
    # Text area for manual input
    batch_input_area = st.text_area("Enter gene names (separated by commas or new lines):", height=150, key="batch_text_area")
    
    # Organism selection for batch search
    organism_options_batch = ["Homo sapiens", "Mus musculus", "Drosophila melanogaster", 
                              "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Danio rerio", "Arabidopsis thaliana"]
    selected_organism_batch = st.selectbox("Select Organism for Batch Search", organism_options_batch, key="organism_batch_search")

    if st.button("Process Batch", key="process_batch_button", use_container_width=True, type="primary"):
        gene_list = []
        if uploaded_file is not None:
            # To read file as string:
            stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            # Process lines from file
            gene_list.extend([line.strip() for line in stringio.readlines() if line.strip()])
        
        if batch_input_area:
            # Process lines from text area
            gene_list.extend([g.strip() for g in batch_input_area.replace("\n", ",").split(",") if g.strip()])
        
        if not gene_list:
            st.warning("Please provide gene names either by uploading a file or entering them manually.")
        else:
            # Remove duplicates and empty strings
            gene_list = sorted(list(set(filter(None, gene_list)))) 
            
            if not gene_list: # Check again after filtering
                st.warning("No valid gene names found after processing input.")
            else:
                st.info(f"Processing {len(gene_list)} unique gene(s) for organism: {selected_organism_batch}")
                results = []
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                for idx, gene in enumerate(gene_list):
                    status_text.text(f"Processing {gene} ({idx+1}/{len(gene_list)})...")
                    annotations = fetch_gene_annotation(gene, selected_organism_batch)
                    if annotations: # fetch_gene_annotation returns a list of results
                        results.extend(annotations) # Use extend for list of dicts
                    else:
                        # Add a placeholder if no annotation found for a specific gene, to keep track
                        results.append({
                            "Gene ID": "N/A", "Gene Symbol": gene, "Official Name": "Not Found", 
                            "Aliases": "N/A", "Function": "N/A", "Chromosome": "N/A", 
                            "Position": "N/A", "Exons": "N/A", "Source Organism": selected_organism_batch,
                            "Annotation Status": "Not Found", "Nucleotide IDs": "N/A"
                        })
                    progress_bar.progress((idx+1)/len(gene_list))
                    time.sleep(0.33) # Respect NCBI API rate limits (3 requests per second without API key)
                
                status_text.success("Batch processing complete!")
                
                if results:
                    df = pd.DataFrame(results)
                    st.success(f"Processed {len(gene_list)} gene(s). Found annotations for {len(df[df['Gene ID'] != 'N/A'])} of them.")
                    st.dataframe(df)
                    
                    csv = df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="Download Batch Results as CSV",
                        data=csv,
                        file_name=f"batch_gene_annotations_{selected_organism_batch.replace(' ', '_')}.csv",
                        mime="text/csv",
                        key="download_batch_csv"
                    )
                else: # Should not happen if placeholders are added
                    st.warning("No results found or an error occurred during batch search.")

# About Tab
elif selected_tab == "📖 About":
    st.markdown("<h1 class='main-header'>About Annotrax</h1>", unsafe_allow_html=True)
    
   # Developer's Desk Section
    with st.container():
        st.markdown("""
        <div class='main-container'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Developer's Desk</h2>
            <div style='display: flex; flex-direction: column; align-items: center; margin: 2rem 0;'>
                <img src='https://media.licdn.com/dms/image/v2/D5603AQFfI1KWVSWl0Q/profile-displayphoto-shrink_400_400/B56ZRRpRfTH0Ag-/0/1736536563970?e=1752710400&v=beta&t=4UtS-OUEjdd0dV7BMRQFA2CxjhOd-tgFPmYbeNAzqB4' 
                     style='width: 180px; height: 180px; object-fit: cover; border-radius: 50%; border: 3px solid #4CAF50; margin-bottom: 1rem;'>
                <h3>Tejal Kale</h3>
                <p style='color: #4CAF50!important; margin-bottom: 0.5rem;'>Developer of Annotrax</p>
                <div class='feature-card' style='max-width: 800px;'>
                    <p style='text-align: justify; line-height: 1.6;'>
                    Hello Bioinformaticians! I am Tejal Kale, a Master's student in Bioinformatics at Deccan Education Society's Pune University. 
                    As someone who's spent countless hours navigating the complexities of gene annotation, I know firsthand the frustration 
                    of juggling multiple tools and databases. That's why I created Annotrax - a platform born out of my own struggles 
                    and passion for bioinformatics.<br><br>
                    I wanted to build a tool that would save others the time and effort I've wasted, and instead empower them to focus 
                    on what really matters: discovering new insights and advancing our understanding of the genetic code. With Annotrax, 
                    I've aimed to create a seamless, intuitive experience that brings together the best of gene annotation and analysis 
                    in one place. My hope is that it becomes an indispensable companion for researchers, students, and scientists, 
                    helping them unlock the secrets of the genome and drive innovation.
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # Mentor Section
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Meet My Mentor</h2>
            <div style='display: flex; flex-direction: column; align-items: center; margin: 2rem 0;'>
                <img src='https://despu.edu.in/assets/admin/img/faculty_images/1738156003_375696d99642e3a89742.jpg' 
                     style='width: 180px; height: 180px; object-fit: cover; border-radius: 50%; border: 3px solid #4CAF50; margin-bottom: 1rem;'>
                <h3>Dr. Kushagra Kashyap</h3>
                <p style='color: #4CAF50!important; margin-bottom: 0.5rem;'>Assistant Professor - DES Pune University</p>
                <div class='feature-card' style='max-width: 800px;'>
                    <p style='text-align: justify; line-height: 1.6;'>
                    Meet Dr. Kushagra Kashyap, a passionate and enthusiastic educator with a PhD from CSIR-CDRI Lucknow. 
                    As an Assistant Professor at Deccan Education Society's Pune University, Dr. Kashyap brings his expertise 
                    in Structural Bioinformatics and Cheminformatics to the classroom, inspiring students to explore the 
                    fascinating world of computational biology and drug discovery. With a strong commitment to academic excellence, 
                    Dr. Kashyap's teaching style is characterized by clarity, enthusiasm, and a genuine passion for mentoring 
                    the next generation of scientists.
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)

  # Acknowledgements Section
    # Acknowledgements Section
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Acknowledgements</h2>
            <div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
                <p style='color: #e0e0e0; line-height: 1.6; font-size: 16px;'>
                I would like to extend my heartfelt gratitude to:<br><br>
                1. Dr. Kushagra Kashyap, for his invaluable guidance and support throughout this project. His expertise and enthusiasm were instrumental in shaping my work.<br><br>
                2. Dr. Poonam Deshpande, for her exceptional guidance, mentorship, and unwavering support. Her insights and feedback were crucial to the project's success.<br><br>
                3. DES Pune University, for providing me with the necessary resources, infrastructure, and opportunities to complete this project.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)

 # Purpose Section
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Purpose</h2>
            <div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
                <p style='color: #e0e0e0; line-height: 1.6; font-size: 16px;'>
                Annotrax was developed to streamline gene annotation workflows, providing researchers with:
                <ul style='color: #e0e0e0; padding-left: 20px;'>
                    <li>Quick access to comprehensive gene information</li>
                    <li>Batch processing capabilities</li>
                    <li>User-friendly interface for both experts and students</li>
                    <li>Centralized biological data from trusted sources</li>
                </ul>
                The tool aims to reduce time spent on data gathering and increase focus on analysis and discovery.
                </p>
            </div>
        </div>
        """, unsafe_allow_html=True)

##Contact Section
    with st.container():
        st.markdown("""
        <div class='main-container' style='margin-top: 2rem;'>
            <h2 style='color: #4CAF50!important; border-bottom: 2px solid #4CAF50; padding-bottom: 0.5rem;'>Contact & Feedback</h2>
<div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
    <h3 style='color: #4CAF50!important;'>Tejal Kale</h3>
    <p style='color: #e0e0e0 ; line-height: 1.6; font-size: 16px;'>
        Master's Student in Bioinformatics<br>
        DES Pune University<br>
        <span style='color: #4CAF50;'>📧 </span><a href="mailto:kaletejal05@mail.com" style='color: #4CAF50;'>kaletejal05@mail.com</a>
    </p>
    <div style='margin-top: 1rem;'>
        <a href='https://www.linkedin.com/in/tejal-kale-0b2195299?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app' target='_blank' style='text-decoration: none; margin-right: 1rem;'>
            <button style='background-color: #e0e0e0; color: white; border: none; padding: 8px 20px; border-radius: 5px; cursor: pointer;'>
                Connect on LinkedIn
            </button>
        </a>
        <a href='https://github.com/tejal-kale' target='_blank' style='text-decoration: none;'>
            <button style='background-color: #e0e0e0; color: white; border: none; padding: 8px 20px; border-radius: 5px; cursor: pointer;'>
                View GitHub Profile
            </button>
        </a>
    </div>
</div>

<!-- Dr. Kushagra Kashyap Contact -->
<div class='feature-card' style='margin: 2rem 0; padding: 1.5rem;'>
    <h3 style='color: #4CAF50!important;'>Dr. Kushagra Kashyap</h3>
    <p style='color: #e0e0e0; line-height: 1.6; font-size: 16px;'>
        Assistant Professor<br>
        DES Pune University
    </p>
    <div style='margin-top: 1rem;'>
        <a href='https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb' target='_blank' style='text-decoration: none;'>
            <button style='background-color: #e0e0e0; color: white; border: none; padding: 8px 20px; border-radius: 5px; cursor: pointer;'>
                Connect on LinkedIn
            </button>
        </a>
    </div>
</div>
""", unsafe_allow_html=True)
# Sidebar footer
st.sidebar.markdown("---")
with st.sidebar.expander("ℹ️ App Information"):
    st.markdown(""" 
    **Last Updated: May 2025** 

    **Developed with ❤ by Tejal Kale for all the Bioinformatics Enthusiasts** 
    
    This tool uses the NCBI Entrez API for fetching gene and sequence data. 
    Please be mindful of API usage limits.
    """)
