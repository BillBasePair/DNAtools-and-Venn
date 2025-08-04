import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os  # For dynamic script name
from matplotlib_venn import venn2, venn3
from venny4py.venny4py import venny4py
from Bio import SeqIO
from itertools import combinations
import tempfile  # For temporary files
import base64  # For embedding images in HTML
from Bio.SeqUtils import MeltingTemp as mt  # For melting temperature calculation

# Function to complement a sequence (5' to 3' direction)
def complement(seq):
    comp = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return seq.translate(comp)

# Function to reverse complement a sequence
def reverse_complement(seq):
    return complement(seq)[::-1]

# Function to check fuzzy match (simplified Levenshtein distance check)
def is_fuzzy_match(seq1, seq2, max_mismatches):
    if len(seq1) != len(seq2):
        return False
    mismatches = sum(a != b for a, b in zip(seq1, seq2))
    return mismatches <= max_mismatches

# Function to extract FASTQ sequences
def extract_fastq_sequences(lines):
    sequences = []
    for i in range(1, len(lines), 4):
        if i + 1 <= len(lines):
            seq = lines[i].strip().upper()
            if seq and not seq.islower():  # Assuming uppercase sequences are valid
                sequences.append(seq)
    return sequences

# Function to highlight all query substrings in full sequence with red text
def highlight_match(full_seq, query):
    result = full_seq
    start = 0
    while True:
        start = result.find(query, start)
        if start == -1:
            break
        end = start + len(query)
        result = result[:start] + '<span style="color: red">' + query + '</span>' + result[end:]
        start += len('<span style="color: red">' + query + '</span>')  # Move past the highlighted part
    return result

# Function to create and return HTML content with both Venn diagram and legend
def create_venn_html(sets, keys):
    # Create the data-dependent Venn diagram
    fig = plt.figure(figsize=(8, 8))
    plt.clf()
    venny4py(sets=sets)
    ax = plt.gca()
    ax.annotate('A', xy=((51 % 10) * 0.1, (51 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('B', xy=((63 % 10) * 0.1, (63 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('C', xy=((74 % 10) * 0.1, (74 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('D', xy=((76 % 10) * 0.1, (76 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('E', xy=((68 % 10) * 0.1, (68 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('F', xy=((47 % 10) * 0.1, (47 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('G', xy=((27 % 10) * 0.1, (27 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('H', xy=((14 % 10) * 0.1, (14 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('I', xy=((43 % 10) * 0.1, (43 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('J', xy=((23 % 10) * 0.1, (23 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('K', xy=((24 % 10) * 0.1, (24 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('L', xy=((26 % 10) * 0.1, (26 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('M', xy=((34 % 10) * 0.1, (30 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('N', xy=((65 % 10) * 0.1, 0.64), color='red', fontsize=10, xycoords='axes fraction')
    ax.annotate('O', xy=((59 % 10) * 0.1, (59 // 10) * 0.1), color='red', fontsize=10, xycoords='axes fraction')

    # Save the data-dependent Venn diagram as PNG
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
        plt.savefig(tmpfile.name, format='png', dpi=300, bbox_inches='tight')
        venn_path = tmpfile.name

    # Convert the data-dependent Venn diagram to base64
    with open(venn_path, "rb") as venn_file:
        venn_base64 = base64.b64encode(venn_file.read()).decode()

    # Convert the legend image to base64
    legend_path = os.path.join(os.path.dirname(__file__), '4way_venn_legend.png')
    with open(legend_path, "rb") as legend_file:
        legend_base64 = base64.b64encode(legend_file.read()).decode()

    # Create HTML content with both images
    html_content = f"""
    <html>
    <body>
        <h2>4-Set Venn Diagram</h2>
        <img src="data:image/png;base64,{venn_base64}" alt="Data-Dependent 4-Set Venn Diagram" style="width:70%;height:auto;">
        <h2>Legend</h2>
        <img src="data:image/png;base64,{legend_base64}" alt="4-Set Venn Legend" style="width:30%;height:auto;">
    </body>
    </html>
    """

    return html_content

# Main Streamlit app
st.set_page_config(page_title="DNA String Tools", layout="wide")
# ----- BEGIN BasePair Secure Password Gate -----
PASSWORD = "IceNine9&"
password = st.text_input("Enter password (same as BasePair wifi password):", type="password")
if password != PASSWORD:
    st.warning("Incorrect password. Please try again.")
    st.stop()
# ----- END Password Gate -----

st.title("🧬 DNA String Tools")
tab1, tab2, tab3, tab4 = st.tabs(["🔍 Sequence Finder", "🔁 Reverse Complement", "🌡️ Melting Temp Calculator", "🔗 Venn Diagrams"])

# Tab 1: Sequence Finder
with tab1:
    with st.expander("About This Tab"):
        st.write("Search for specific DNA/RNA subsequences (min 4 bases) within uploaded files, highlighting all occurrences in the full sequence for easy identification. Use 'Run Search' to compare exact, partial, or fuzzy matches, with results downloadable as CSV.")

    st.header("Sequence Finder")
    query_input_method = st.radio("Input method", ["Upload file", "Paste sequences"])
    query_seqs = []

    if query_input_method == "Upload file":
        file = st.file_uploader("Upload query file", type=["csv", "xlsx", "fastq", "fq", "txt", "fasta"])
        if file:
            name = file.name.lower()
            if name.endswith(("fastq", "fq")):
                lines = file.read().decode("utf-8").splitlines()
                n = st.number_input("Max reads", 1000, 1000000, 10000, 1000)
                query_seqs = extract_fastq_sequences(lines, n)
            elif name.endswith("csv"):
                df = pd.read_csv(file)
                col = st.selectbox("Sequence column", df.columns, index=df.columns.get_loc("Trimmed") if "Trimmed" in df.columns else 0)
                query_seqs = df[col].dropna().astype(str).str.upper().tolist()
            elif name.endswith("xlsx"):
                df = pd.read_excel(file)
                col = st.selectbox("Sequence column", df.columns, index=df.columns.get_loc("Trimmed") if "Trimmed" in df.columns else 0)
                query_seqs = df[col].dropna().astype(str).str.upper().tolist()
            else:
                text = file.read().decode("utf-8")
                query_seqs = [line.strip().upper() for line in text.splitlines() if line and not line.startswith(">")]
    else:
        pasted = st.text_area("Paste sequences")
        query_seqs = [line.strip().upper() for line in pasted.splitlines() if line.strip()]
        st.write(f"Debug: Query sequences loaded: {query_seqs}")  # Debug output

    search_files = st.file_uploader("Upload search files", type=["csv", "xlsx", "fastq", "fq"], accept_multiple_files=True)
    match_type = st.radio("Match type", ["Exact", "Partial", "Fuzzy"])
    allow_rc = st.checkbox("Match reverse complement", value=True)
    mismatches = st.slider("Mismatches (fuzzy only)", 1, 5, 1) if match_type == "Fuzzy" else 0

    if st.button("Run Search"):
        if query_seqs and search_files:
            results = []
            for f in search_files:
                fname = f.name
                ext = fname.lower().split(".")[-1]
                if ext in ["csv", "xlsx"]:
                    df = pd.read_csv(f) if ext == "csv" else pd.read_excel(f)
                    st.write(f"**{fname} preview:**")
                    st.dataframe(df.head())
                    col = st.selectbox(f"Column in {fname}", df.columns, key=f"col_{fname}", index=df.columns.get_loc("Trimmed") if "Trimmed" in df.columns else 0)
                    id_col = st.selectbox(f"Optional ID column in {fname}", ["(None)"] + list(df.columns), key=f"id_{fname}")
                    for idx, row in df.iterrows():
                        val = str(row[col]).strip().upper()
                        seqs = [val, reverse_complement(val)] if allow_rc else [val]
                        for q in query_seqs:
                            for s in seqs:
                                if match_type == "Exact" and q == s or \
                                   match_type == "Partial" and q in s or \
                                   match_type == "Fuzzy" and is_fuzzy_match(q, s, mismatches):
                                    highlighted = s  # Use full sequence
                                    start = 0
                                    while True:
                                        start = highlighted.find(q, start)
                                        if start == -1:
                                            break
                                        end = start + len(q)
                                        highlighted = highlighted[:start] + '<span style="color: red">' + q + '</span>' + highlighted[end:]
                                        start += len('<span style="color: red">' + q + '</span>')  # Move past the highlighted part
                                    results.append({"Query": q, "Match": s, "Highlighted Match": highlighted, "File": fname, "Row": idx,
                                                    "ID": row[id_col] if id_col != "(None)" else ""})
                                    break
                elif ext in ["fastq", "fq"]:
                    lines = f.read().decode("utf-8").splitlines()
                    for i in range(1, len(lines), 4):
                        val = lines[i].strip().upper()
                        seqs = [val, reverse_complement(val)] if allow_rc else [val]
                        for q in query_seqs:
                            for s in seqs:
                                if match_type == "Exact" and q == s or \
                                   match_type == "Partial" and q in s or \
                                   match_type == "Fuzzy" and is_fuzzy_match(q, s, mismatches):
                                    highlighted = s  # Use full sequence
                                    start = 0
                                    while True:
                                        start = highlighted.find(q, start)
                                        if start == -1:
                                            break
                                        end = start + len(q)
                                        highlighted = highlighted[:start] + '<span style="color: red">' + q + '</span>' + highlighted[end:]
                                        start += len('<span style="color: red">' + q + '</span>')  # Move past the highlighted part
                                    results.append({"Query": q, "Match": s, "Highlighted Match": highlighted, "File": fname, "Row": i // 4, "ID": ""})
                                    break

            df_out = pd.DataFrame(results)
            html = df_out.to_html(index=False, escape=False)  # Render HTML for red text
            st.markdown(html, unsafe_allow_html=True)
            if not df_out.empty:
                st.download_button("Download results", df_out.to_csv(index=False), "matches.csv", "text/csv")
        else:
            st.warning("Please upload query sequences and search files before running the search.")

# Tab 2: Reverse Complement
with tab2:
    st.header("Reverse Complement")
    seq_input = st.text_area("DNA or RNA Sequence, 5' to 3'")
    if seq_input:
        # Use st.table with HTML for bold red styling in column headers
        data = [{"Input (5' to 3')": s, "Reverse (5' to 3')": s[::-1], "Complement <b style='color: red'>(3' to 5')</b>": complement(s), "Reverse Complement (5' to 3')": reverse_complement(s)} for s in seq_input.splitlines() if s.strip()]
        df = pd.DataFrame(data)
        # Convert DataFrame to HTML with styled header
        html = df.to_html(index=False, escape=False)  # escape=False to render HTML
        st.markdown(html, unsafe_allow_html=True)

# Tab 3: Melting Temp Calculator
with tab3:
    with st.expander("About This Tab"):
        st.write("**Nearest Neighbor Tables:**")
        st.write("- DNA_NN3: SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proceedings of the National Academy of Sciences of the United States of America, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460")
        st.write("- DNA_NN4: SantaLucia, J., & Hicks, D. (2004). The thermodynamics of DNA structural motifs. Annual Review of Biophysics and Biomolecular Structure, 33, 415–440. https://doi.org/10.1146/annurev.biophys.32.110601.141800")
        st.write("- RNA_NN2: Freier, S. M., Kierzek, R., Jaeger, J. A., Sugimoto, N., Caruthers, M. H., Neilson, T., & Turner, D. H. (1986). Improved free-energy parameters for predictions of RNA duplex stability. Proceedings of the National Academy of Sciences of the United States of America, 83(24), 9373–9377. https://doi.org/10.1073/pnas.83.24.9373")
        st.write("- RNA_NN3: Xia, T., SantaLucia, J., Burkard, M. E., Kierzek, R., Schroeder, S. J., Jiao, X., Cox, C., & Turner, D. H. (1998). Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs. Biochemistry, 37(41), 14719–14735. https://doi.org/10.1021/bi9809425")
        st.write("**Salt Correction Methods:**")
        st.write("- Schildkraut (1965): Schildkraut, C. (1965). Dependence of the melting temperature of DNA on salt concentration. Biopolymers, 3(2), 195–208. https://doi.org/10.1002/bip.360030207")
        st.write("- SantaLucia (1998): SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proceedings of the National Academy of Sciences of the United States of America, 95(4), 1460–1465. https://doi.org/10.1073/pnas.95.4.1460")
        st.write("- Owczarzy (2004): Owczarzy, R., You, Y., Moreira, B. G., Owczarzy, J. A., Grollman, L. G., Behlke, M. A., & Walder, J. A. (2004). Effects of sodium ions on DNA duplex oligomers: improved predictions of melting temperatures. Biochemistry, 43(12), 3537–3554. https://doi.org/10.1021/bi034621r")
        st.write("- Owczarzy (2008): Owczarzy, R., Tataurov, A. V., Wu, Y., Manthey, J. A., McQuisten, K. A., Almabrazi, H. G., Pedersen, K. F., Lin, Y., Garretson, J., McEntaggart, N. G., Sailor, C. A., Dawson, R. B., & Peek, A. S. (2008). IDT SciTools: a suite for analysis and design of oligonucleotide based molecular diagnostics. Nucleic Acids Research, 36(Web Server issue), W163–W169. https://doi.org/10.1093/nar/gkn161")

    st.header("Melting Temperature Calculator")
    seq_input = st.text_area("Enter sequence (5' to 3')").upper()
    duplex_type = st.selectbox("Duplex type", ["DNA/DNA", "DNA/RNA", "RNA/RNA"])
    nn_table = st.selectbox("Nearest Neighbor Table", ["DNA_NN3 (SantaLucia 1998)", "DNA_NN4 (SantaLucia 2004)", "RNA_NN2 (Freier 1986)", "RNA_NN3 (Xia 1998)"])
    Na = st.number_input("Na+ concentration (mM)", value=137.0)  # Updated to 137 mM for PBS
    K = st.number_input("K+ concentration (mM)", value=0.0)
    Tris = st.number_input("Tris concentration (mM)", value=0.0)
    Mg = st.number_input("Mg2+ concentration (mM)", value=0.0)
    dNTPs = st.number_input("dNTPs concentration (mM)", value=0.0)
    saltcorr = st.selectbox("Salt correction method", ["Schildkraut (1965)", "SantaLucia (1998)", "Owczarzy (2004)", "Owczarzy (2008)"])

    if st.button("Calculate Tm"):
        if duplex_type == "DNA/DNA":
            check = 'dna'
        elif duplex_type == "DNA/RNA":
            check = 'dna'
        elif duplex_type == "RNA/RNA":
            check = 'rna'
        if nn_table == "DNA_NN3 (SantaLucia 1998)":
            nn_table = mt.DNA_NN3
        elif nn_table == "DNA_NN4 (SantaLucia 2004)":
            nn_table = mt.DNA_NN4
        elif nn_table == "RNA_NN2 (Freier 1986)":
            nn_table = mt.RNA_NN2
        elif nn_table == "RNA_NN3 (Xia 1998)":
            nn_table = mt.RNA_NN3
        if saltcorr == "Schildkraut (1965)":
            saltcorr = 1
        elif saltcorr == "SantaLucia (1998)":
            saltcorr = 4
        elif saltcorr == "Owczarzy (2004)":
            saltcorr = 5
        elif saltcorr == "Owczarzy (2008)":
            saltcorr = 7

        Tm = mt.Tm_NN(seq_input, check=check, nn_table=nn_table, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=saltcorr)
        st.write(f"Melting Temperature (Tm): {Tm:.2f} °C")
        st.write("Reference: SantaLucia J Jr (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci U S A 95(3):1460-1465.")

# Tab 4: Venn Diagrams with Sequence Extraction
with tab4:
    with st.expander("About This Tab"):
        st.write("Visualize overlaps between 2-4 sequence sets from uploaded files, with downloadable 2-, 3-, and 4-set diagrams and selectable intersection regions for export of those sequences.  Use 'Run Comparison' to analyze exact or fuzzy matches.  Sequences for any region in a 2-, 3-, or 4-way Venn Diagram (corresponding to two, three, or four files) can be extracted by selecting the appropriate region in the drop-down menu on the left side pane.  More than a single region can be selected. Please don't forget to scroll all the way down after hitting the button for \"Run X-way Venn Comparison\". Four-way Venn diagrams will not render in this same browser tab so you must open the separate HTML tab that gets generated for 4-way Venn diagrams.  Region selections don't show up the first time without making at least ONE column selection in the file previews, OR run it a 2nd time.  A selection will populate the dropdown with a red bar showing the intended selection then hit run again.  What can I say, I'm not a real programmer - Bill Jackson, July 2025.")
        st.write("[See detailed guide](https://github.com/your-repo/README.md#venn-diagram)")

    st.header("Overlap Viewer")
    # Sidebar for controls and intersection selection
    with st.sidebar:
        script_name = os.path.basename(__file__)  # Dynamically get the script filename
        st.write(f"Name of script: {script_name}", style={"font-size": "12px"})  # Updated to dynamic name
        st.header("Venn Diagram Controls")
        venn_files = st.file_uploader("Upload files for overlap", type=["csv", "xlsx", "fastq", "fq"], accept_multiple_files=True)
        if venn_files and 2 <= len(venn_files) <= 4:
            match_type = st.radio("Match type", ["Exact", "Fuzzy"], key="vennmatch")
            mismatches = st.slider("Allowed mismatches", 1, 5, 1, key="mismatch") if match_type == "Fuzzy" else 0
            allow_rc = st.checkbox("Match reverse complement", value=False)  # Disabled temporarily
        # Populate intersection selection dynamically
        intersection_labels = st.session_state.get("intersection_labels", [])
        selected_intersections = st.multiselect("Select intersections to extract sequences", intersection_labels, key="intersection_select")
# ----- BEGIN 4-WAY LEGEND IN SIDEBAR -----
        if len(venn_files) == 4:
            st.markdown("**Legend for selecting intersections - RED LETTERS - from 4-diagrams. Numbers are from an arbitrary example.**")
            legend_path = os.path.join(os.path.dirname(__file__), '4way_venn_legend.png')
            if os.path.exists(legend_path):
                st.image(legend_path, width=300)
        # ----- END 4-WAY LEGEND IN SIDEBAR -----

    if venn_files and 2 <= len(venn_files) <= 4:
        sets = {}
        for i, file in enumerate(venn_files):
            ext = file.name.split(".")[-1].lower()
            label = st.text_input(f"Label for file {i+1}", value=file.name, key=f"label_{i}")
            if ext in ["csv", "xlsx"]:
                df = pd.read_csv(file) if ext == "csv" else pd.read_excel(file)
                st.dataframe(df.head())
                default_col = "Trimmed" if "Trimmed" in df.columns else df.columns[0]
                col = st.selectbox(f"Sequence column in {label}", df.columns, index=df.columns.get_loc(default_col))
                seqs = df[col].dropna().astype(str).str.upper().tolist()
            elif ext in ["fastq", "fq"]:
                lines = file.read().decode("utf-8").splitlines()
                seqs = extract_fastq_sequences(lines)
            else:
                seqs = []

            if allow_rc:
                seqs += [reverse_complement(s) for s in seqs]
            sets[label] = set(seqs)

        # Generate intersections and labels as soon as files are uploaded
        num_sets = len(venn_files)
        keys = list(sets.keys())
        if num_sets == 2:
            intersection_labels = ['1 (only ' + keys[0] + ')', '2 (only ' + keys[1] + ')', '3 (' + keys[0] + ' & ' + keys[1] + ')']
            intersections = {}
            # 1 = S1 only
            inter = sets[keys[0]].copy()
            others = set(keys) - set([keys[0]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[0]] = inter
            # 2 = S2 only
            inter = sets[keys[1]].copy()
            others = set(keys) - set([keys[1]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[1]] = inter
            # 3 = S1 & S2
            inter = sets[keys[0]] & sets[keys[1]]
            intersections[intersection_labels[2]] = inter
            fig = plt.figure()
            venn = venn = venn2([sets[keys[0]], sets[keys[1]]], set_labels=keys)
            from matplotlib.patches import Patch
            colors = ['tab:blue', 'tab:orange']
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            st.pyplot(fig)
        elif num_sets == 3:
            intersection_labels = ['1 (only ' + keys[0] + ')', '2 (only ' + keys[1] + ')', '3 (only ' + keys[2] + ')',
                                  '4 (' + keys[0] + ' & ' + keys[1] + ')', '5 (' + keys[0] + ' & ' + keys[2] + ')',
                                  '6 (' + keys[1] + ' & ' + keys[2] + ')', '7 (' + keys[0] + ' & ' + keys[1] + ' & ' + keys[2] + ')']
            intersections = {}
            # 1 = S1 only
            inter = sets[keys[0]].copy()
            others = set(keys) - set([keys[0]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[0]] = inter
            # 2 = S2 only
            inter = sets[keys[1]].copy()
            others = set(keys) - set([keys[1]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[1]] = inter
            # 3 = S3 only
            inter = sets[keys[2]].copy()
            others = set(keys) - set([keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[2]] = inter
            # 4 = S1 & S2
            inter = sets[keys[0]] & sets[keys[1]]
            others = set(keys) - set([keys[0], keys[1]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[3]] = inter
            # 5 = S1 & S3
            inter = sets[keys[0]] & sets[keys[2]]
            others = set(keys) - set([keys[0], keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[4]] = inter
            # 6 = S2 & S3
            inter = sets[keys[1]] & sets[keys[2]]
            others = set(keys) - set([keys[1], keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[5]] = inter
            # 7 = S1, 2, & 3
            inter = sets[keys[0]] & sets[keys[1]] & sets[keys[2]]
            intersections[intersection_labels[6]] = inter
            fig = plt.figure()
            venn = venn = venn3([sets[keys[0]], sets[keys[1]], sets[keys[2]]], set_labels=keys)
            from matplotlib.patches import Patch
            colors = ['tab:blue', 'tab:orange', 'tab:green']
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99', '#9999ff']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99', '#9999ff']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            st.pyplot(fig)
        elif num_sets == 4:
            intersection_labels = [f"A (only {keys[0]})", f"B ({keys[0]} & {keys[1]})", f"C (only {keys[1]})", f"D (only {keys[2]})", f"E ({keys[2]} & {keys[3]})", f"F ({keys[1]} & {keys[2]} & {keys[3]})", f"G ({keys[1]} & {keys[3]})", f"H ({keys[0]} & {keys[3]})", f"I ({keys[0]} & {keys[1]} & {keys[2]})", f"J ({keys[0]} & {keys[2]})", f"K ({keys[0]} & {keys[2]} & {keys[3]})", f"L ({keys[0]} & {keys[1]} & {keys[3]})", f"M ({keys[0]} & {keys[1]} & {keys[2]} & {keys[3]})", f"N ({keys[1]} & {keys[2]})", f"O (only {keys[3]})"]
            intersections = {}
            # A = S1 only
            inter = sets[keys[0]].copy()
            others = set(keys) - set([keys[0]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[0]] = inter
            # B = S1 & S2
            inter = sets[keys[0]] & sets[keys[1]]
            others = set(keys) - set([keys[0], keys[1]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[1]] = inter
            # C = S2 only
            inter = sets[keys[1]].copy()
            others = set(keys) - set([keys[1]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[2]] = inter
            # D = S3 only
            inter = sets[keys[2]].copy()
            others = set(keys) - set([keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[3]] = inter
            # E = S3 & S4
            inter = sets[keys[2]] & sets[keys[3]]
            others = set(keys) - set([keys[2], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[4]] = inter
            # F = S2, 3, & 4
            inter = sets[keys[1]] & sets[keys[2]] & sets[keys[3]]
            others = set(keys) - set([keys[1], keys[2], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[5]] = inter
            # G = S2 & S4
            inter = sets[keys[1]] & sets[keys[3]]
            others = set(keys) - set([keys[1], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[6]] = inter
            # H = S1 & S4
            inter = sets[keys[0]] & sets[keys[3]]
            others = set(keys) - set([keys[0], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[7]] = inter
            # I = S1, 2, & 3
            inter = sets[keys[0]] & sets[keys[1]] & sets[keys[2]]
            others = set(keys) - set([keys[0], keys[1], keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[8]] = inter
            # J = S1 & 3
            inter = sets[keys[0]] & sets[keys[2]]
            others = set(keys) - set([keys[0], keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[9]] = inter
            # K = S1, 3 & 4
            inter = sets[keys[0]] & sets[keys[2]] & sets[keys[3]]
            others = set(keys) - set([keys[0], keys[2], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[10]] = inter
            # L = S1, 2 & 4
            inter = sets[keys[0]] & sets[keys[1]] & sets[keys[3]]
            others = set(keys) - set([keys[0], keys[1], keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[11]] = inter
            # M = S1, 2, 3, & 4 (ALL)
            inter = set.intersection(*sets.values())
            intersections[intersection_labels[12]] = inter
            # N = S2 & 3
            inter = sets[keys[1]] & sets[keys[2]]
            others = set(keys) - set([keys[1], keys[2]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[13]] = inter
            # O = S4 only
            inter = sets[keys[3]].copy()
            others = set(keys) - set([keys[3]])
            inter = inter - set.union(*(sets[o] for o in others))
            intersections[intersection_labels[14]] = inter

        # Update session state with intersection labels
        st.session_state["intersection_labels"] = intersection_labels

        if st.button(f"Run {num_sets}-Set Venn Comparison"):
            all_seqs = set().union(*sets.values())
            membership_data = {s: {name: any(is_fuzzy_match(s, t, mismatches) if match_type == "Fuzzy" else s == t for t in sets[name]) for name in sets} for s in all_seqs}
            df_memberships = pd.DataFrame(membership_data).T
            present_combinations = df_memberships[df_memberships.any(axis=1)]

            st.subheader(f"{num_sets}-Set Venn Diagram")
            fig = plt.figure()
            if num_sets == 2:
                venn = venn = venn2([sets[keys[0]], sets[keys[1]]], set_labels=keys)
            elif num_sets == 3:
                venn = venn = venn3([sets[keys[0]], sets[keys[1]], sets[keys[2]]], set_labels=keys)
            elif num_sets == 4:
                html_content = create_venn_html(sets, keys)
                st.download_button("Download 4-Set Venn with Legend", data=html_content, file_name="4_set_venn_with_legend.html", mime="text/html")
            from matplotlib.patches import Patch
            colors = ['tab:blue', 'tab:orange']
            from matplotlib.patches import Patch
            colors = ['tab:blue', 'tab:orange', 'tab:green']
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99', '#9999ff']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            from matplotlib.patches import Patch
            venn_colors = ['#ff9999', '#99ff99', '#9999ff']
            if len(keys) >= 2:
                patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                plt.legend(handles=patches, loc='upper right')
            st.pyplot(fig)

            # Handle selection
            selected_intersections = st.session_state.get("intersection_select", [])
            if selected_intersections:
                fig = plt.figure()
                if num_sets == 2:
                    venn = venn = venn2([sets[keys[0]], sets[keys[1]]], set_labels=keys)
                elif num_sets == 3:
                    venn = venn = venn3([sets[keys[0]], sets[keys[1]], sets[keys[2]]], set_labels=keys)
                elif num_sets == 4:
                    html_content = create_venn_html(sets, keys)
                    st.download_button("Download Selected 4-Set Venn with Legend", data=html_content, file_name="selected_4_set_venn_with_legend.html", mime="text/html")
                from matplotlib.patches import Patch
                colors = ['tab:blue', 'tab:orange']
                from matplotlib.patches import Patch
                colors = ['tab:blue', 'tab:orange', 'tab:green']
                from matplotlib.patches import Patch
                venn_colors = ['#ff9999', '#99ff99']
                if len(keys) >= 2:
                    patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                    plt.legend(handles=patches, loc='upper right')
                from matplotlib.patches import Patch
                venn_colors = ['#ff9999', '#99ff99', '#9999ff']
                if len(keys) >= 2:
                    patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                    plt.legend(handles=patches, loc='upper right')
                from matplotlib.patches import Patch
                venn_colors = ['#ff9999', '#99ff99']
                if len(keys) >= 2:
                    patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                    plt.legend(handles=patches, loc='upper right')
                from matplotlib.patches import Patch
                venn_colors = ['#ff9999', '#99ff99', '#9999ff']
                if len(keys) >= 2:
                    patches = [Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:len(keys)])]
                    plt.legend(handles=patches, loc='upper right')
                st.pyplot(fig)

                # Extract and display sequences based on selected intersection
                selected_seqs = set()
                for label in selected_intersections:
                    selected_seqs.update(intersections.get(label, set()))
                if selected_seqs:
                    st.subheader("Sequences in Selected Intersections")
                    col1, col2 = st.columns([1, 2])
                    with col2:
                        df_seqs = pd.DataFrame({'Sequence': list(selected_seqs)})
                        st.dataframe(df_seqs)
                        st.download_button("Download sequences as CSV", df_seqs.to_csv(index=False), file_name="selected_sequences.csv")
                else:
                    st.write("No sequences in the selected intersections.")
            else:
                st.warning("This app supports 2 to 4 sets. Please upload the correct number of files.")
    else:
        if venn_files and (len(venn_files) < 2 or len(venn_files) > 4):
            st.warning("Please upload 2 to 4 files for Venn diagram analysis.")