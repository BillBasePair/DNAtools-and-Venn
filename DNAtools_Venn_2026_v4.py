import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os  # For dynamic script name
from matplotlib_venn import venn2, venn3
from Bio import SeqIO
from itertools import combinations
import tempfile  # For temporary files
import base64  # For embedding images in HTML
from Bio.SeqUtils import MeltingTemp as mt  # For melting temperature calculation

# Try to import upsetplot for 4/5-way overlap diagrams
try:
    from upsetplot import UpSet, from_memberships
    UPSETPLOT_AVAILABLE = True
except ImportError:
    UPSETPLOT_AVAILABLE = False

# Try to import pyvenn for true geometric 4/5/6-way Venn diagrams
try:
    from venn import venn
    PYVENN_AVAILABLE = True
except ImportError:
    PYVENN_AVAILABLE = False

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

# ── Intersection label definitions ────────────────────────────────────────────
# For N sets, there are 2^N - 1 possible non-empty intersections.
# These dicts map a frozenset of set-indices → letter label.
# 4-way: A–O (15 regions), matching the original legend image.
LABELS_4WAY = {
    frozenset([0]):          "A",
    frozenset([0,1]):        "B",
    frozenset([1]):          "C",
    frozenset([2]):          "D",
    frozenset([2,3]):        "E",
    frozenset([1,2,3]):      "F",
    frozenset([1,3]):        "G",
    frozenset([0,3]):        "H",
    frozenset([0,1,2]):      "I",
    frozenset([0,2]):        "J",
    frozenset([0,2,3]):      "K",
    frozenset([0,1,3]):      "L",
    frozenset([0,1,2,3]):    "M",
    frozenset([1,2]):        "N",
    frozenset([3]):          "O",
}

# 5-way: A–AE (31 regions) — letters A–Z then AA–AE
def _gen_5way_labels():
    from itertools import combinations
    letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ") + ["AA","AB","AC","AD","AE"]
    idx = 0
    d = {}
    for r in range(1, 6):
        for combo in combinations(range(5), r):
            d[frozenset(combo)] = letters[idx]
            idx += 1
    return d
LABELS_5WAY = _gen_5way_labels()


def _membership_to_label(member_indices, num_sets):
    """Return the letter label for a given frozenset of set indices."""
    key = frozenset(member_indices)
    if num_sets == 4:
        return LABELS_4WAY.get(key, "?")
    elif num_sets == 5:
        return LABELS_5WAY.get(key, "?")
    return "?"


def build_intersections(sets, keys):
    """Compute all 2^N-1 exact intersections and return as ordered dict
    keyed by the dropdown label string e.g. 'A (only file1.csv)'."""
    from itertools import combinations as _comb
    n = len(keys)
    label_map = LABELS_4WAY if n == 4 else (LABELS_5WAY if n == 5 else {})
    intersections = {}
    intersection_labels = []
    for r in range(1, n + 1):
        for combo in _comb(range(n), r):
            combo_set = frozenset(combo)
            letter = label_map.get(combo_set, str(combo))
            # Human-readable description
            names = [keys[i] for i in combo]
            if r == 1:
                desc = f"{letter} (only {names[0]})"
            else:
                desc = f"{letter} (" + " & ".join(names) + ")"
            # Exact exclusive intersection
            inter = set.intersection(*(sets[keys[i]] for i in combo))
            outside = [keys[i] for i in range(n) if i not in combo]
            if outside:
                inter = inter - set.union(*(sets[k] for k in outside))
            intersections[desc] = inter
            intersection_labels.append(desc)
    return intersections, intersection_labels


# ── UpSet plot with A–O / A–AE letter labels on bars ─────────────────────────
def create_upset_figure(sets, keys, intersections):
    """Build an UpSet plot and annotate each bar with its region letter.
    Letters are placed BELOW the dot matrix so they never collide with counts.
    Works for 4 or 5 sets. Returns matplotlib Figure or None."""
    if not UPSETPLOT_AVAILABLE:
        return None
    n = len(keys)
    label_map = LABELS_4WAY if n == 4 else (LABELS_5WAY if n == 5 else {})

    all_seqs = set().union(*sets.values())
    memberships = []
    for seq in all_seqs:
        member_of = tuple(k for k in keys if seq in sets[k])
        if member_of:
            memberships.append(member_of)

    data = from_memberships(memberships)
    # Extra bottom margin so the letter row fits under the dot matrix
    fig = plt.figure(figsize=(14, 7))
    fig.subplots_adjust(bottom=0.18)
    upset = UpSet(data, subset_size="count", show_counts=True, sort_by="cardinality")
    axes_dict = upset.plot(fig)
    plt.suptitle(f"{n}-Set Overlap (UpSet Plot) — red letters match dropdown",
                 fontsize=12, y=1.02)

    # ── Place letter labels in the dot-matrix axis, one row below the lowest dots ──
    # "matrix" axis holds the dot grid; we annotate along its x positions.
    matrix_ax = axes_dict.get("matrix")
    bar_ax    = axes_dict.get("intersections")

    if bar_ax is not None and matrix_ax is not None:
        # Derive the sorted column order from the data index (same sort as UpSet)
        counts = data.groupby(level=list(range(n))).size().sort_values(ascending=False)

        # Get x-centre of each bar from the bar axis
        bars = bar_ax.patches
        for bar, (idx_tuple, _) in zip(bars, counts.items()):
            member_indices = frozenset(i for i, v in enumerate(idx_tuple) if v)
            letter = label_map.get(member_indices, "?")
            # x in bar_ax coordinates
            x_bar = bar.get_x() + bar.get_width() / 2
            # Convert to matrix_ax coordinates (same figure, x-axes are aligned)
            # Place letter just below the bottom of the matrix axis (y < 0 in axes coords)
            matrix_ax.text(
                x_bar, -0.7, letter,
                ha="center", va="top",
                fontsize=10, fontweight="bold", color="red",
                transform=matrix_ax.get_xaxis_transform(),
                clip_on=False
            )

    return fig


# Fallback bar chart when upsetplot not installed
def create_summary_figure(intersections, keys):
    labels = list(intersections.keys())
    counts = [len(v) for v in intersections.values()]
    fig, ax = plt.subplots(figsize=(10, max(6, len(labels) * 0.4)))
    bars = ax.barh(labels, counts, color="steelblue")
    ax.bar_label(bars, padding=3)
    ax.set_xlabel("Sequence count")
    ax.set_title(f"{len(keys)}-Set Overlap Region Sizes")
    ax.invert_yaxis()
    plt.tight_layout()
    return fig


# pyvenn geometric diagram (4 or 5-way true Venn ellipses)
def create_pyvenn_figure(sets, keys):
    """Use pyvenn to draw a geometric 4- or 5-way Venn. Returns Figure or None."""
    if not PYVENN_AVAILABLE:
        return None
    fig, ax = plt.subplots(figsize=(10, 8))
    venn(sets, ax=ax)
    ax.set_title(f"{len(keys)}-Set Venn Diagram", fontsize=13)
    plt.tight_layout()
    return fig

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
        st.caption(f"Name of script: {script_name}")  # Updated to dynamic name
        st.header("Venn Diagram Controls")
        venn_files = st.file_uploader("Upload files for overlap", type=["csv", "xlsx", "fastq", "fq"], accept_multiple_files=True)
        if venn_files and 2 <= len(venn_files) <= 5:
            match_type = st.radio("Match type", ["Exact", "Fuzzy"], key="vennmatch")
            mismatches = st.slider("Allowed mismatches", 1, 5, 1, key="mismatch") if match_type == "Fuzzy" else 0
            allow_rc = st.checkbox("Match reverse complement", value=False)  # Disabled temporarily
        # Populate intersection selection dynamically
        intersection_labels = st.session_state.get("intersection_labels", [])
        selected_intersections = st.multiselect("Select intersections to extract sequences", intersection_labels, key="intersection_select")
# ----- BEGIN 4/5-WAY NOTE IN SIDEBAR -----
        if venn_files and len(venn_files) in (4, 5):
            n = len(venn_files)
            if PYVENN_AVAILABLE:
                st.info(f"{n}-way: geometric Venn diagram (pyvenn). Letters A–{'O' if n==4 else 'AE'} match dropdown.")
            elif UPSETPLOT_AVAILABLE:
                st.info(f"{n}-way: UpSet plot with letter labels. Letters A–{'O' if n==4 else 'AE'} match dropdown.")
            else:
                st.warning("Install **upsetplot** or **venn** (pyvenn) for richer diagrams. Bar-chart fallback active.")
        # ----- END 4/5-WAY NOTE IN SIDEBAR -----

    if venn_files and 2 <= len(venn_files) <= 5:
        sets = {}
        for i, file in enumerate(venn_files):
            ext = file.name.split(".")[-1].lower()
            label = st.text_input(f"Label for file {i+1}", value=file.name, key=f"label_{i}")
            if ext in ["csv", "xlsx"]:
                df = pd.read_csv(file) if ext == "csv" else pd.read_excel(file)
                st.dataframe(df.head())
                default_col = "Trimmed" if "Trimmed" in df.columns else df.columns[0]
                col = st.selectbox(f"Sequence column in {label}", df.columns,
                                   index=df.columns.get_loc(default_col))
                seqs = df[col].dropna().astype(str).str.upper().tolist()
            elif ext in ["fastq", "fq"]:
                lines = file.read().decode("utf-8").splitlines()
                seqs = extract_fastq_sequences(lines)
            else:
                seqs = []
            if allow_rc:
                seqs += [reverse_complement(s) for s in seqs]
            sets[label] = set(seqs)

        num_sets = len(venn_files)
        keys = list(sets.keys())

        # ── Compute all intersections via shared helper ─────────────────────
        intersections, intersection_labels = build_intersections(sets, keys)

        # Preview diagram before Run button (2- and 3-way only)
        venn_colors = ['#ff9999', '#99ff99', '#9999ff', '#ffcc99', '#cc99ff']
        if num_sets == 2:
            fig, ax = plt.subplots()
            venn2([sets[keys[0]], sets[keys[1]]], set_labels=keys, ax=ax)
            patches = [mpatches.Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:2])]
            ax.legend(handles=patches, loc='upper right')
            st.pyplot(fig);  plt.close(fig)
        elif num_sets == 3:
            fig, ax = plt.subplots()
            venn3([sets[keys[0]], sets[keys[1]], sets[keys[2]]], set_labels=keys, ax=ax)
            patches = [mpatches.Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:3])]
            ax.legend(handles=patches, loc='upper right')
            st.pyplot(fig);  plt.close(fig)

        # Update sidebar dropdown
        st.session_state["intersection_labels"] = intersection_labels

        if st.button(f"Run {num_sets}-Set Venn Comparison"):

            st.subheader(f"{num_sets}-Set Overlap Diagram")

            if num_sets == 2:
                fig, ax = plt.subplots()
                venn2([sets[keys[0]], sets[keys[1]]], set_labels=keys, ax=ax)
                patches = [mpatches.Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:2])]
                ax.legend(handles=patches, loc='upper right')
                st.pyplot(fig);  plt.close(fig)

            elif num_sets == 3:
                fig, ax = plt.subplots()
                venn3([sets[keys[0]], sets[keys[1]], sets[keys[2]]], set_labels=keys, ax=ax)
                patches = [mpatches.Patch(color=c, label=k) for k, c in zip(keys, venn_colors[:3])]
                ax.legend(handles=patches, loc='upper right')
                st.pyplot(fig);  plt.close(fig)

            elif num_sets in (4, 5):
                # Priority: pyvenn geometric > upsetplot > bar-chart fallback
                if PYVENN_AVAILABLE:
                    st.info("Showing geometric Venn (pyvenn). Letters on bars match the dropdown.")
                    fig_pv = create_pyvenn_figure(sets, keys)
                    if fig_pv is not None:
                        st.pyplot(fig_pv);  plt.close(fig_pv)
                    # Also show UpSet underneath if available for letter reference
                    if UPSETPLOT_AVAILABLE:
                        st.markdown("**UpSet plot with region letters (use these to pick from dropdown):**")
                        fig_up = create_upset_figure(sets, keys, intersections)
                        if fig_up is not None:
                            st.pyplot(fig_up);  plt.close(fig_up)
                elif UPSETPLOT_AVAILABLE:
                    fig_up = create_upset_figure(sets, keys, intersections)
                    if fig_up is not None:
                        st.pyplot(fig_up);  plt.close(fig_up)
                    else:
                        st.warning("UpSet plot failed to render.")
                else:
                    st.info("Install **upsetplot** or **venn** (pyvenn) for richer diagrams. "
                            "Showing bar-chart fallback.")
                    fig_fb = create_summary_figure(intersections, keys)
                    st.pyplot(fig_fb);  plt.close(fig_fb)

            # ── Sequence extraction from selected intersections ──────────────
            selected_intersections = st.session_state.get("intersection_select", [])
            if selected_intersections:
                selected_seqs = set()
                for lbl in selected_intersections:
                    selected_seqs.update(intersections.get(lbl, set()))
                if selected_seqs:
                    st.subheader("Sequences in Selected Intersections")
                    df_seqs = pd.DataFrame({"Sequence": sorted(selected_seqs)})
                    st.dataframe(df_seqs)
                    st.download_button("Download sequences as CSV",
                                       df_seqs.to_csv(index=False),
                                       file_name="selected_sequences.csv")
                else:
                    st.write("No sequences found in the selected intersections.")
    else:
        if venn_files and (len(venn_files) < 2 or len(venn_files) > 5):
            st.warning("Please upload 2 to 5 files for Venn diagram analysis.")