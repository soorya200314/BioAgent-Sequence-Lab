"""Streamlit Web Application for DNA Sequence Analysis and Phylogeny Visualization."""
import io
import streamlit as st
import pandas as pd
import numpy as np
# Import our core bioinformatics engine
from bio_engine.alignment import needleman_wunsch, calculate_identity_and_distance
from visualization.tree_generator import generate_tree

# Set up page configurations
st.set_page_config(
    page_title="BioAgent Sequence Lab",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Mock dataset (HBB - Hemoglobin Beta subunit segment)
MOCK_FASTA = """>Human
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTGA
>Neanderthal
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGACGAAGTTGGTGGTGAGGCCCTGGGCAGGTGA
>Chimpanzee
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTATTGCCCTGTGGGGCAAGGTGAACGTGGAAGAAGTTGGTGGTGAGGCCCTGGGCAGGTGA
>Gorilla
ATGGTGCACCTGACTCCTGAGGAGAAGTCCGCCGTTATTGCCCTGTGGGGCAAGGTGAACGTGGAAGAAGTTGGTGGTCAGGCCCTGGGCAGGTGA
>Mouse
ATGGTCCTTACTGAGGATAAGTCTGCCGTTACTGCCGCCCTGTGGGGCATGGTGAACGTGGATGAAGTTGCTGGTGAGGCCCTGGGCAGGTGA
"""

def parse_fasta_string(fasta_str: str) -> dict[str, str]:
    """Parse FASTA contents from a string into a dictionary of headers to sequences."""
    sequences = {}
    current_header = None
    current_seq = []
    for line in fasta_str.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                sequences[current_header] = "".join(current_seq).upper()
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)
    if current_header:
        sequences[current_header] = "".join(current_seq).upper()
    return sequences

# App Custom Styled Banner
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap');
    
    /* Apply Outfit font to headers */
    .css-10trblm, .st-emotion-cache-10trblm, h1, h2, h3 {
        font-family: 'Outfit', sans-serif !important;
    }
    
    .banner {
        background: linear-gradient(135deg, #1A365D 0%, #2B6CB0 100%);
        padding: 1.8rem;
        border-radius: 12px;
        color: white;
        margin-bottom: 2rem;
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.15);
    }
    .banner h1 {
        margin: 0;
        font-size: 2.4rem;
        font-weight: 700;
        letter-spacing: -0.5px;
    }
    .banner p {
        margin: 0.4rem 0 0 0;
        font-size: 1.1rem;
        opacity: 0.9;
        font-weight: 300;
    }
    </style>
    <div class="banner">
        <h1>BioAgent Sequence Lab 🧬</h1>
        <p>Interactive Genomic Sequence Alignment, Mutation Detection, and Phylogenetic Inference</p>
    </div>
    """,
    unsafe_allow_html=True
)

# Sidebar UI
st.sidebar.image(
    "https://img.icons8.com/external-flatart-icons-flat-flatarticons/128/000000/external-dna-medical-health-flatart-icons-flat-flatarticons.png",
    width=80
)
st.sidebar.markdown("### Data Upload & Parameters")

# Sidebar File Uploader
uploaded_file = st.sidebar.file_uploader(
    "Upload FASTA File (.fasta, .fa)",
    type=["fasta", "fa"],
    help="Upload a multi-sequence DNA FASTA file."
)

# Load sequences from uploaded file or fall back to mock data
if uploaded_file is not None:
    try:
        fasta_content = uploaded_file.getvalue().decode("utf-8")
        sequences = parse_fasta_string(fasta_content)
        st.sidebar.success(f"Successfully loaded {len(sequences)} sequences from uploaded file.")
    except Exception as e:
        st.sidebar.error(f"Error parsing FASTA: {e}")
        st.sidebar.info("Falling back to default mock dataset.")
        sequences = parse_fasta_string(MOCK_FASTA)
else:
    st.sidebar.info("Using built-in mock hemoglobin (HBB) dataset.")
    sequences = parse_fasta_string(MOCK_FASTA)

# Error out if too few sequences loaded
if len(sequences) < 2:
    st.error("The loaded dataset must contain at least 2 sequences to perform analysis.")
    st.stop()

headers = list(sequences.keys())

# Select wild-type reference
ref_name = st.sidebar.selectbox(
    "Select Wild-Type (Baseline) Reference:",
    options=headers,
    index=0,
    help="The baseline sequence against which mutations will be compared."
)

st.sidebar.markdown("---")
st.sidebar.markdown(
    """
    **Alignment Parameters (Global)**
    - Match Score: `+2`
    - Mismatch Penalty: `-1`
    - Gap Penalty: `-2`
    """
)

# Main App Layout (Tabs)
tab1, tab2 = st.tabs([
    "🔍 Micro-Mutational Analysis Dashboard", 
    "🌳 Macro-Evolutionary Phylogenetic Tree"
])

# Tab 1: Micro-Mutational Analysis
with tab1:
    st.subheader("Point Mutation Scanning Report")
    st.markdown(
        f"Detecting nucleotide-level variations in all query sequences against **'{ref_name}'**."
    )
    ref_seq = sequences[ref_name]
    point_mutations = []

    # Run alignments and collect mutations
    with st.spinner("Analyzing sequence variations..."):
        for query_name in headers:
            if query_name == ref_name:
                continue
            query_seq = sequences[query_name]
            # Perform global alignment
            aligned_ref, aligned_query, _ = needleman_wunsch(ref_seq, query_seq)
            
            # Map columns to reference codon positions
            alignment_len = len(aligned_ref)
            codon_map = [0] * alignment_len
            ref_char_count = 0
            for col in range(alignment_len):
                if aligned_ref[col] != "-":
                    codon_map[col] = (ref_char_count // 3) + 1
                    ref_char_count += 1
            
            # Fill in gap columns in reference alignment
            last_codon = 1
            for col in range(alignment_len):
                if codon_map[col] > 0:
                    last_codon = codon_map[col]
                else:
                    codon_map[col] = last_codon
            
            # Detect point mutations column-by-column
            for col in range(alignment_len):
                ref_base = aligned_ref[col]
                mut_base = aligned_query[col]
                if ref_base != mut_base:
                    # Classify mutation type
                    if ref_base != "-" and mut_base != "-":
                        mut_type = "Substitution"
                    elif ref_base == "-":
                        mut_type = "Insertion"
                    else:
                        mut_type = "Deletion"
                    
                    # FIXED: Added the proper math definition before appending to the list
                    codon_pos_in_triplet = (col % 3) + 1
                    
                    point_mutations.append({
                        "Species Name": query_name,
                        "Sequence Alignment Position": col + 1,
                        "Codon Number": codon_map[col],
                        "Codon Position": codon_pos_in_triplet,
                        "Reference Base": ref_base,
                        "Mutated Base": mut_base,
                        "Mutation Type": mut_type
                    })

    if point_mutations:
        df_mut = pd.DataFrame(point_mutations)
        
        # Display summary metrics
        m_col1, m_col2, m_col3 = st.columns(3)
        m_col1.metric("Total Point Mutations", len(df_mut))
        m_col2.metric("Substitutions", len(df_mut[df_mut["Mutation Type"] == "Substitution"]))
        m_col3.metric("Indels (Ins/Del)", len(df_mut[df_mut["Mutation Type"] != "Substitution"]))
        
        # Sort options
        st.markdown("### Mutations Data Table")
        st.dataframe(
            df_mut.sort_values(by=["Species Name", "Sequence Alignment Position"]), # FIXED: Added missing comma
            use_container_width=True,
            hide_index=True
        )
        
        # Download CSV option
        csv_data = df_mut.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Download Mutation Report (CSV)",
            data=csv_data,
            file_name=f"mutation_report_ref_{ref_name}.csv",
            mime="text/csv",
            help="Download the mutations table as a CSV file."
        )
    else:
        st.success("All comparison sequences are 100% identical to the reference. No mutations detected.")

# Tab 2: Macro-Evolutionary Phylogenetic Tree
with tab2:
    st.subheader("Evolutionary Distance & Hierarchical Inference")
    st.markdown(
        "Pairwise evolutionary distances are calculated using $1.0 - \\text{Sequence Identity}$ "
        "obtained from optimal Needleman-Wunsch alignments. Tree construction uses the UPGMA method."
    )
    num_seqs = len(headers)
    dist_matrix = [[0.0] * num_seqs for _ in range(num_seqs)]
    
    # Compute distance matrix
    with st.spinner("Computing distance matrix..."):
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                seq_i = sequences[headers[i]]
                seq_j = sequences[headers[j]]
                aligned_i, aligned_j, _ = needleman_wunsch(seq_i, seq_j)
                _, distance = calculate_identity_and_distance(aligned_i, aligned_j)
                dist_matrix[i][j] = distance
                dist_matrix[j][i] = distance
                
    # Display distance matrix in UI
    st.markdown("### Pairwise Distance Matrix")
    df_dist = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    st.dataframe(df_dist.style.background_gradient(cmap="Blues"), use_container_width=True)
    
    # Plot and display the tree
    st.markdown("### Reconstructed UPGMA Phylogenetic Tree")
    with st.spinner("Generating tree plot..."):
        fig = generate_tree(dist_matrix, headers)
        st.pyplot(fig)
        
        # Save plot to in-memory bytes for downloading
        img_buf = io.BytesIO()
        fig.savefig(img_buf, format="png", dpi=300, facecolor="#FAFAFA")
        img_data = img_buf.getvalue()
        st.download_button(
            label="Download Phylogenetic Tree (PNG)",
            data=img_data,
            file_name="phylogeny_tree.png",
            mime="image/png",
            help="Download this high-resolution tree visualization."
        )