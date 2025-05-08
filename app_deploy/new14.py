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
    .main-header {
        font-size: 2.5rem;
        color: #1E6091;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #1E6091;
    }
    .result-box {
        background-color: #f0f2f6;
        border-radius: 5px;
        padding: 20px;
        margin-bottom: 20px;
    }
    .info-text {
        color: #31708f;
    }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 class='main-header'>🧬 Annotrax: Gene Annotation Explorer</h1>", unsafe_allow_html=True)
st.markdown("Search for genes and download annotations easily - saving time for students and researchers")

# Create a cache for gene data to avoid repeated API calls
@st.cache_data(ttl=3600)  # Cache for 1 hour
def fetch_gene_annotation(gene_name, organism="Homo sapiens"):
    """Fetch gene annotations from NCBI using Entrez."""
    try:
        # Search NCBI Gene
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_term, retmax=5) # Limit to 5 to avoid too many results
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            return None
        
        results = []
        for gene_id in record["IdList"]:
            # Get gene summary
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(summary_handle)
            summary_handle.close()
            
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
            
            # Get nucleotide sequence info if available
            nuccore_ids = []
            try:
                link_handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id)
                link_results = Entrez.read(link_handle)
                link_handle.close()
                
                if link_results[0]["LinkSetDb"]:
                    for link_info in link_results[0]["LinkSetDb"]: # Iterate through LinkSetDb entries
                        if link_info["DbTo"] == "nuccore":
                            for link in link_info["Link"]:
                                nuccore_ids.append(link["Id"])
                            break # Assuming we only care about the first set of nuccore links
            except Exception: # Broad except to catch any linking issues
                pass # nuccore_ids will remain empty
                
            # Gene location info
            genomic_info = docsum.get("GenomicInfo", [{}])
            chromosome = genomic_info[0].get("ChrLoc", "N/A") if genomic_info else "N/A"
            chr_start = genomic_info[0].get("ChrStart", "N/A") if genomic_info else "N/A"
            chr_stop = genomic_info[0].get("ChrStop", "N/A") if genomic_info else "N/A"
            
            # Gather annotation data
            gene_data = {
                "Gene ID": gene_id,
                "Gene Symbol": docsum.get("Name", gene_name),
                "Official Name": docsum.get("Description", "N/A"),
                "Aliases": ", ".join(docsum.get("OtherAliases", "").split(", ")[:5]), # Limit aliases
                "Function": docsum.get("Summary", "No function description available"),
                "Chromosome": chromosome,
                "Position": f"{chr_start}-{chr_stop}" if chr_start != "N/A" and chr_stop != "N/A" else "N/A",
                "Exons": docsum.get("GenomicInfo", [{}])[0].get("ExonCount", "N/A") if genomic_info and docsum.get("GenomicInfo")[0] else "N/A",
                "Source Organism": organism,
                "Annotation Status": "Complete" if len(docsum.get("Summary", "")) > 50 else "Partial",
                "Nucleotide IDs": ", ".join(nuccore_ids[:3]) if nuccore_ids else "N/A" # Limit displayed IDs
            }
            
            results.append(gene_data)
            
        return results
    
    except Exception as e:
        st.error(f"Error fetching from NCBI: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_gene_sequence(nuccore_id):
    """Fetch gene sequence (GenBank format) and features from NCBI Nucleotide database."""
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
        
        sequence_data = {
            "id": record.id,
            "name": record.name,
            "description": record.description,
            "length": len(record.seq),
            "features": features[:10],  # Limit to first 10 features for performance
            "sequence": str(record.seq)[:1000] + ("..." if len(record.seq) > 1000 else "")
        }
        
        return sequence_data
    
    except Exception as e:
        st.error(f"Error fetching GenBank sequence: {e}")
        return None

@st.cache_data(ttl=3600)
def fetch_fasta_sequence_data(nuccore_id):
    """Fetch gene sequence in FASTA format from NCBI Nucleotide database."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="fasta", retmode="text")
        fasta_sequence = handle.read()
        handle.close()
        return fasta_sequence
    except Exception as e:
        st.error(f"Error fetching FASTA sequence: {e}")
        return None

# Function to extract specific features from gene data
def extract_feature_data(gene_data, feature_type):
    result = {}
    if feature_type == "function":
        result["Function"] = gene_data.get("Function", "")
    elif feature_type == "exons":
        result["Exons"] = gene_data.get("Exons", "")
        result["Position"] = gene_data.get("Position", "")
    elif feature_type == "organism":
        result["Source Organism"] = gene_data.get("Source Organism", "")
    # "Gene Symbol" is added by the calling function if needed
    # If "all", the entire gene_data dict is used directly
    return result

# Sidebar for search options
st.sidebar.markdown("<h2 class='sub-header'>Search Options</h2>", unsafe_allow_html=True)

# Database selection
database_option = st.sidebar.radio(
    "Select Database",
    ["All Genes", "Fully Annotated Only", "Partially Annotated"]
)

# Organism selection
organism_options = ["Homo sapiens", "Mus musculus", "Drosophila melanogaster", 
                   "Caenorhabditis elegans", "Saccharomyces cerevisiae"]
selected_organism = st.sidebar.selectbox("Select Organism", organism_options)

# Search mode
search_mode = st.sidebar.radio(
    "Search Mode",
    ["Single Gene", "Multiple Genes (Batch)"]
)

# Main search interface
if search_mode == "Single Gene":
    gene_name_input = st.text_input("🔍 Enter gene name (e.g., INS, BRCA1, TP53):") # Renamed to avoid conflict
    
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
                    tabs = st.tabs([f"{g['Gene Symbol']} ({g['Gene ID']})" for g in gene_annotations])
                    
                    for i, tab in enumerate(tabs):
                        with tab:
                            gene_data_tab = gene_annotations[i]
                            st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                            
                            col1, col2 = st.columns([3, 1])
                            
                            with col1:
                                st.markdown(f"### {gene_data_tab['Gene Symbol']} - {gene_data_tab['Official Name']}")
                                st.markdown(f"**Aliases:** {gene_data_tab['Aliases']}")
                                st.markdown(f"**Location:** Chromosome {gene_data_tab['Chromosome']}, Position {gene_data_tab['Position']}")
                                st.markdown(f"**Exons:** {gene_data_tab['Exons']}")
                                st.markdown(f"**Function:**")
                                st.markdown(f"{gene_data_tab['Function']}")
                            
                            with col2:
                                st.markdown(f"**Organism:** {gene_data_tab['Source Organism']}")
                                st.markdown(f"**Gene ID:** {gene_data_tab['Gene ID']}")
                                st.markdown(f"**Annotation:** {gene_data_tab['Annotation Status']}")
                                
                                if gene_data_tab['Nucleotide IDs'] != "N/A":
                                    first_nuccore_id = gene_data_tab['Nucleotide IDs'].split(", ")[0]
                                    gene_id_key_part = gene_data_tab['Gene ID']

                                    if st.button(f"View Sequence Data", key=f"seq_view_tab_{gene_id_key_part}_{i}"):
                                        with st.spinner("Fetching GenBank details..."):
                                            seq_data = fetch_gene_sequence(first_nuccore_id)
                                        if seq_data:
                                            st.markdown("#### Sequence Information (GenBank)")
                                            st.markdown(f"**ID:** {seq_data['id']}")
                                            st.markdown(f"**Length:** {seq_data['length']} bp")
                                            with st.expander("View Sequence (First 1000 bp)"):
                                                st.text(seq_data['sequence'])
                                            with st.expander("View Features (first 10)"):
                                                for feat in seq_data['features']:
                                                    st.markdown(f"**{feat['type']}** at {feat['location']}")
                                                    for k, v_list in feat['qualifiers'].items():
                                                        if k in ['product', 'gene', 'note', 'protein_id']:
                                                            st.markdown(f"- {k}: {v_list}")
                                        else:
                                            st.warning("Could not fetch GenBank sequence data.")
                                    
                                    if st.button(f"Download FASTA", key=f"fasta_dl_tab_{gene_id_key_part}_{i}"):
                                        with st.spinner(f"Fetching FASTA for {first_nuccore_id}..."):
                                            fasta_content = fetch_fasta_sequence_data(first_nuccore_id)
                                        if fasta_content:
                                            st.download_button(
                                                label="⬇️ FASTA Ready",
                                                data=fasta_content,
                                                file_name=f"{gene_data_tab['Gene Symbol']}_{first_nuccore_id}.fasta",
                                                mime="text/fasta",
                                                key=f"actual_fasta_dl_tab_{gene_id_key_part}_{i}"
                                            )
                                        else:
                                            st.warning("Could not fetch FASTA sequence.")
                            st.markdown("</div>", unsafe_allow_html=True)
                else: # Single result display
                    gene_data_single = gene_annotations[0]
                    st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.markdown(f"### {gene_data_single['Gene Symbol']} - {gene_data_single['Official Name']}")
                        st.markdown(f"**Aliases:** {gene_data_single['Aliases']}")
                        st.markdown(f"**Location:** Chromosome {gene_data_single['Chromosome']}, Position {gene_data_single['Position']}")
                        st.markdown(f"**Exons:** {gene_data_single['Exons']}")
                        st.markdown(f"**Function:**")
                        st.markdown(f"{gene_data_single['Function']}")
                    with col2:
                        st.markdown(f"**Organism:** {gene_data_single['Source Organism']}")
                        st.markdown(f"**Gene ID:** {gene_data_single['Gene ID']}")
                        st.markdown(f"**Annotation:** {gene_data_single['Annotation Status']}")
                        
                        if gene_data_single['Nucleotide IDs'] != "N/A":
                            first_nuccore_id = gene_data_single['Nucleotide IDs'].split(", ")[0]
                            gene_id_key_part = gene_data_single['Gene ID']

                            if st.button(f"View Sequence Data", key=f"seq_view_single_{gene_id_key_part}"):
                                with st.spinner("Fetching GenBank details..."):
                                    seq_data = fetch_gene_sequence(first_nuccore_id)
                                if seq_data:
                                    st.markdown("#### Sequence Information (GenBank)")
                                    st.markdown(f"**ID:** {seq_data['id']}")
                                    st.markdown(f"**Length:** {seq_data['length']} bp")
                                    with st.expander("View Sequence (First 1000 bp)"):
                                        st.text(seq_data['sequence'])
                                    with st.expander("View Features (first 10)"):
                                        for feat in seq_data['features']:
                                            st.markdown(f"**{feat['type']}** at {feat['location']}")
                                            for k, v_list in feat['qualifiers'].items():
                                                if k in ['product', 'gene', 'note', 'protein_id']:
                                                    st.markdown(f"- {k}: {v_list}")
                                else:
                                    st.warning("Could not fetch GenBank sequence data.")

                            if st.button(f"Download FASTA", key=f"fasta_dl_single_{gene_id_key_part}"):
                                with st.spinner(f"Fetching FASTA for {first_nuccore_id}..."):
                                    fasta_content = fetch_fasta_sequence_data(first_nuccore_id)
                                if fasta_content:
                                    st.download_button(
                                        label="⬇️ FASTA Ready",
                                        data=fasta_content,
                                        file_name=f"{gene_data_single['Gene Symbol']}_{first_nuccore_id}.fasta",
                                        mime="text/fasta",
                                        key=f"actual_fasta_dl_single_{gene_id_key_part}"
                                    )
                                else:
                                    st.warning("Could not fetch FASTA sequence.")
                    st.markdown("</div>", unsafe_allow_html=True)
                
                # Download options for annotations
                st.markdown("<h3 class='sub-header'>📥 Download Annotation Options</h3>", unsafe_allow_html=True)
                
                download_genes_for_anno = gene_annotations # Default to all found & filtered genes
                if len(gene_annotations) > 1:
                    selected_indices_anno = []
                    st.markdown("Select genes for annotation download (defaults to all displayed):")
                    for i, gene in enumerate(gene_annotations):
                        # Use gene ID for key to ensure stability if list changes
                        if st.checkbox(f"Select {gene['Gene Symbol']} ({gene['Gene ID']}) for download", value=True, key=f"select_anno_{gene['Gene ID']}"):
                            selected_indices_anno.append(i)
                    
                    # If user unchecks all, this logic might need adjustment or provide a message.
                    # For now, if some are selected, use them. If selection leads to empty, could default to all or warn.
                    # This example assumes if checkboxes are shown, selection drives the download.
                    # A simpler approach might be to always download all displayed results, or selected if any are *unchecked*.
                    # Let's refine: download selected if any are explicitly checked, otherwise all.
                    # This requires a bit more logic to see if any checkbox was touched.
                    # For simplicity, let's make it: if checkboxes appear, user MUST select. Or, default to all.
                    # The current code: if selected_indices_anno is empty, download_genes uses all. If populated, it uses selection.
                    # This is slightly confusing. Let's make it explicit.
                    
                    # Let's make the selection more straightforward:
                    # Create a list of gene data to be downloaded based on checkboxes
                    temp_download_genes = []
                    any_selected = False
                    for i, gene_item in enumerate(gene_annotations):
                         # Checkboxes should default to True if user wants to download all initially.
                         # Let's assume for multiple results, we provide checkboxes and user picks.
                         # If no selection, maybe disable download button or warn.
                         # For this example, let's assume the original logic: if selected_indices exists, use it.
                         # The st.checkbox needs a default. If we want "select all by default", then value=True.
                         pass # Using the original logic with gene_annotations / selected_indices_anno

                    # Re-evaluating selection logic for annotation download:
                    # If multiple genes are shown, allow user to pick which ones to include in the annotation download file.
                    # Default to all selected.
                    
                    genes_to_consider_for_download = []
                    if len(gene_annotations) > 1:
                        st.markdown("Select genes to include in the annotation download file:")
                        for i, gene_item in enumerate(gene_annotations):
                            if st.checkbox(f"Include {gene_item['Gene Symbol']} ({gene_item['Gene ID']})", value=True, key=f"dl_sel_{gene_item['Gene ID']}"):
                                genes_to_consider_for_download.append(gene_item)
                        if not genes_to_consider_for_download:
                            st.warning("No genes selected for annotation download.")
                            # Potentially disable download button here or handle this case in download logic
                    else: # Single gene result
                        genes_to_consider_for_download = gene_annotations
                
                else: # Single gene already filtered
                    genes_to_consider_for_download = gene_annotations

                if genes_to_consider_for_download: # Proceed if there are genes to download
                    dl_col1, dl_col2 = st.columns(2)
                    with dl_col1:
                        st.markdown("### Select fields for annotation download:")
                        download_all_fields = st.checkbox("Download all annotation fields", value=True, key="anno_dl_all_fields")
                        
                        download_function_field = False
                        download_exons_field = False
                        download_organism_field = False
                        if not download_all_fields:
                            download_function_field = st.checkbox("Gene functions", value=True, key="anno_dl_func")
                            download_exons_field = st.checkbox("Exon information", value=False, key="anno_dl_exons")
                            download_organism_field = st.checkbox("Source organism", value=False, key="anno_dl_org")
                    
                    with dl_col2:
                        st.markdown("### Annotation download format:")
                        download_format_anno = st.radio("", ["CSV", "JSON", "TXT"], key="anno_dl_format")
                    
                    if st.button("Download Selected Annotations", key="dl_anno_button"):
                        download_data_anno = []
                        for gene_to_dl in genes_to_consider_for_download:
                            if download_all_fields:
                                download_data_anno.append(gene_to_dl)
                            else:
                                gene_export = {"Gene Symbol": gene_to_dl["Gene Symbol"], "Gene ID": gene_to_dl["Gene ID"]} # Always include symbol and ID
                                if download_function_field:
                                    gene_export.update(extract_feature_data(gene_to_dl, "function"))
                                if download_exons_field:
                                    gene_export.update(extract_feature_data(gene_to_dl, "exons"))
                                if download_organism_field:
                                    gene_export.update(extract_feature_data(gene_to_dl, "organism"))
                                # Ensure all keys from extract_feature_data are present if not already
                                for key in ["Function", "Exons", "Position", "Source Organism"]:
                                    if key not in gene_export: # if not selected, fill with N/A or empty
                                         if key in gene_to_dl and (download_function_field or download_exons_field or download_organism_field) : # only add if that category was chosen at all
                                            pass # it's already handled by extract_feature_data
                                         # else: gene_export[key] = "N/A" # or decide not to add the column
                                download_data_anno.append(gene_export)
                        
                        if not download_data_anno:
                             st.warning("No data prepared for download. Check selections.")
                        else:
                            df_anno = pd.DataFrame(download_data_anno)
                            dl_filename_base = gene_name_input if len(genes_to_consider_for_download) == 1 else "selected_genes"

                            if download_format_anno == "CSV":
                                csv_data = df_anno.to_csv(index=False).encode("utf-8")
                                st.download_button(label="⬇️ Download CSV", data=csv_data, file_name=f"{dl_filename_base}_annotations.csv", mime="text/csv")
                            elif download_format_anno == "JSON":
                                json_data = df_anno.to_json(orient="records", indent=2)
                                st.download_button(label="⬇️ Download JSON", data=json_data, file_name=f"{dl_filename_base}_annotations.json", mime="application/json")
                            else:  # TXT
                                txt_data = df_anno.to_string(index=False)
                                st.download_button(label="⬇️ Download TXT", data=txt_data, file_name=f"{dl_filename_base}_annotations.txt", mime="text/plain")
                elif len(gene_annotations) > 1 and not genes_to_consider_for_download : # Only show if multiple genes were possible but none selected
                     st.info("Select at least one gene to enable annotation download options.")

        else:
            st.warning(f"⚠️ No gene data found for '{gene_name_input}'. Please check spelling or try another gene.")

else: # Batch search mode
    st.markdown("### Batch Gene Search")
    st.markdown("Enter multiple gene names separated by commas, spaces, or new lines")
    
    batch_input_genes = st.text_area("Enter gene names:", height=100, key="batch_gene_input")
    
    if batch_input_genes:
        gene_names_batch = []
        for item in batch_input_genes.replace(",", " ").replace("\n", " ").split():
            if item.strip():
                gene_names_batch.append(item.strip())
        
        gene_names_batch = list(set(gene_names_batch)) # Remove duplicates
        
        if gene_names_batch:
            st.markdown(f"Found {len(gene_names_batch)} unique gene names to search.")
            
            if st.button("Search Batch", key="search_batch_button"):
                all_results_batch = []
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                for i, gene_name_item in enumerate(gene_names_batch):
                    status_text.text(f"Searching for {gene_name_item}... ({i+1}/{len(gene_names_batch)})")
                    time.sleep(0.34) # Adhere to NCBI rate limits (max 3/sec without key)
                    gene_annotations_item = fetch_gene_annotation(gene_name_item, selected_organism)
                    
                    if gene_annotations_item:
                        filtered_annotations = gene_annotations_item
                        if database_option == "Fully Annotated Only":
                            filtered_annotations = [g for g in gene_annotations_item if g["Annotation Status"] == "Complete"]
                        elif database_option == "Partially Annotated":
                            filtered_annotations = [g for g in gene_annotations_item if g["Annotation Status"] == "Partial"]
                        all_results_batch.extend(filtered_annotations)
                    
                    progress_bar.progress((i + 1) / len(gene_names_batch))
                
                if all_results_batch:
                    status_text.success(f"Batch search complete! Found {len(all_results_batch)} gene records matching your criteria.")
                    
                    df_display_batch = pd.DataFrame([{
                        "Gene Symbol": g.get("Gene Symbol", "N/A"),
                        "Name": g.get("Official Name", "N/A"),
                        "Organism": g.get("Source Organism", "N/A"),
                        "Chromosome": g.get("Chromosome", "N/A"),
                        "Exons": g.get("Exons", "N/A"),
                        "Status": g.get("Annotation Status", "N/A")
                    } for g in all_results_batch])
                    st.dataframe(df_display_batch)
                    
                    st.markdown("<h3 class='sub-header'>📥 Batch Annotation Download Options</h3>", unsafe_allow_html=True)
                    
                    b_col1, b_col2 = st.columns(2)
                    with b_col1:
                        st.markdown("### Select fields for download:")
                        batch_dl_all = st.checkbox("Download all fields", value=True, key="batch_dl_all_fields")
                        
                        batch_dl_func, batch_dl_exons, batch_dl_org = False, False, False
                        if not batch_dl_all:
                            batch_dl_func = st.checkbox("Gene functions", value=True, key="batch_dl_func")
                            batch_dl_exons = st.checkbox("Exon information", value=False, key="batch_dl_exons")
                            batch_dl_org = st.checkbox("Source organism", value=False, key="batch_dl_org")
                    
                    with b_col2:
                        st.markdown("### Download format:")
                        batch_dl_format = st.radio("", ["CSV", "JSON", "TXT"], key="batch_dl_format_select")
                    
                    if st.button("Download Batch Annotations", key="dl_batch_anno_button"):
                        batch_dl_data = []
                        for gene_res in all_results_batch:
                            if batch_dl_all:
                                batch_dl_data.append(gene_res)
                            else:
                                gene_export_batch = {"Gene Symbol": gene_res.get("Gene Symbol", "N/A"), "Gene ID": gene_res.get("Gene ID", "N/A")}
                                if batch_dl_func: gene_export_batch.update(extract_feature_data(gene_res, "function"))
                                if batch_dl_exons: gene_export_batch.update(extract_feature_data(gene_res, "exons"))
                                if batch_dl_org: gene_export_batch.update(extract_feature_data(gene_res, "organism"))
                                batch_dl_data.append(gene_export_batch)
                        
                        if not batch_dl_data:
                            st.warning("No data prepared for batch download.")
                        else:
                            batch_df_dl = pd.DataFrame(batch_dl_data)
                            if batch_dl_format == "CSV":
                                csv_batch = batch_df_dl.to_csv(index=False).encode("utf-8")
                                st.download_button(label="⬇️ Download CSV", data=csv_batch, file_name="batch_gene_annotations.csv", mime="text/csv")
                            elif batch_dl_format == "JSON":
                                json_batch = batch_df_dl.to_json(orient="records", indent=2)
                                st.download_button(label="⬇️ Download JSON", data=json_batch, file_name="batch_gene_annotations.json", mime="application/json")
                            else:  # TXT
                                txt_batch = batch_df_dl.to_string(index=False)
                                st.download_button(label="⬇️ Download TXT", data=txt_batch, file_name="batch_gene_annotations.txt", mime="text/plain")
                else:
                    status_text.warning("⚠️ No genes found matching your criteria in the batch search.")
        else:
            st.info("Enter gene names to begin batch search.")

# Footer
st.markdown("---")
st.markdown("""
<div class="info-text">
<small>Annotrax retrieves data from NCBI's Entrez API. Please use responsibly and respect NCBI's usage guidelines. Cache active for 1 hour.</small>
</div>
""", unsafe_allow_html=True)

# About section
with st.sidebar.expander("About Annotrax"):
    st.markdown("""
    **Annotrax** is a tool for exploring gene annotations.
    - Search single or multiple genes.
    - View details like function, location, exons.
    - Filter by organism and annotation status.
    - Download annotation data in CSV, JSON, or TXT.
    - View sequence details (GenBank) and download FASTA sequences.
    
    Powered by NCBI Entrez API.
    """)