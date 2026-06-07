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

# v6: Try to import ViennaRNA for the new DNA & RNA Folding tab.
# Guarded so the rest of the app still works if ViennaRNA isn't installed.
try:
    import RNA as _RNA
    from math import sqrt as _sqrt
    VIENNARNA_AVAILABLE = True
except ImportError:
    VIENNARNA_AVAILABLE = False

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

    # ── Patch ax.scatter to sanitize edgecolors before upsetplot passes them ──
    # upsetplot 0.9.0 passes pandas Series as edgecolors; matplotlib 3.9+ rejects this.
    import matplotlib.axes._axes as _mpl_axes
    import numpy as _np
    import pandas as _pd
    _orig_scatter = _mpl_axes.Axes.scatter
    def _patched_scatter(self, *args, **kwargs):
        if "edgecolors" in kwargs:
            ec = kwargs["edgecolors"]
            if isinstance(ec, (_pd.Series, _pd.Index)):
                kwargs["edgecolors"] = ec.tolist()
            elif hasattr(ec, "values"):
                kwargs["edgecolors"] = list(ec)
        return _orig_scatter(self, *args, **kwargs)
    _mpl_axes.Axes.scatter = _patched_scatter

    axes_dict = upset.plot(fig)

    # Restore original scatter after plot
    _mpl_axes.Axes.scatter = _orig_scatter
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
st.set_page_config(page_title="DNA String Tools", page_icon="🧬", layout="wide")
# ----- BEGIN BasePair Secure Password Gate -----
PASSWORD = "IceNine9&"
password = st.text_input("Enter password (same as BasePair wifi password):", type="password")
if password != PASSWORD:
    st.warning("Incorrect password. Please try again.")
    st.stop()
# ----- END Password Gate -----

st.title("🧬 DNA String Tools")
tab1, tab2, tab3, tab4, tab5 = st.tabs(["🔍 Seq Search and/or ST Mapping", "🔁 Reverse Complement", "🧬 DNA & RNA Folding", "🌡️ Melting Temp Calculator", "🔗 Venn Diagrams"])

# Tab 1: Seq Search and/or ST Mapping
with tab1:
    with st.expander("About This Tab"):
        st.write(
            "Search for DNA/RNA sequences within one or more reference files, or map between "
            "paired values such as Stormtrooper tag codes and their corresponding DNA sequences. "
            "Upload a query file (or paste sequences directly) and one or more reference files, "
            "then choose your mode: **Find** searches for your query sequences within a reference "
            "column using exact, partial, or fuzzy matching; **Map** performs a key→value lookup "
            "returning the paired column entry for each query. Results are displayed inline and "
            "downloadable as Excel or CSV. Multiple reference files can be searched simultaneously. "
            "Tip: file previews appear immediately after upload to help you identify the right "
            "columns before running."
        )

    st.header("Seq Search and/or ST Mapping")

    # ── Step 1: Query ──────────────────────────────────────────────────────────
    st.markdown("### Step 1 — Your Query")
    t1_input_method = st.radio(
        "Input method",
        ["Upload a file", "Paste sequences"],
        key="t1_input_method",
        horizontal=True,
    )

    t1_query_seqs = []   # list of plain strings

    if t1_input_method == "Upload a file":
        t1_query_file = st.file_uploader(
            "Upload query file",
            type=["xlsx", "csv", "fastq", "fq", "fasta", "txt"],
            key="t1_query_file",
        )
        if t1_query_file:
            t1_qname = t1_query_file.name.lower()
            if t1_qname.endswith(("fastq", "fq")):
                t1_lines = t1_query_file.read().decode("utf-8").splitlines()
                t1_query_seqs = extract_fastq_sequences(t1_lines)
                st.caption(f"{len(t1_query_seqs):,} sequences loaded from FASTQ")
            elif t1_qname.endswith(("fasta", "fa")):
                t1_text = t1_query_file.read().decode("utf-8")
                t1_query_seqs = [
                    line.strip().upper()
                    for line in t1_text.splitlines()
                    if line.strip() and not line.startswith(">")
                ]
                st.caption(f"{len(t1_query_seqs):,} sequences loaded from FASTA")
            elif t1_qname.endswith("txt"):
                t1_text = t1_query_file.read().decode("utf-8")
                t1_query_seqs = [l.strip().upper() for l in t1_text.splitlines() if l.strip()]
                st.caption(f"{len(t1_query_seqs):,} sequences loaded from TXT")
            else:
                # xlsx or csv — sheet selector + column selector + preview
                if t1_qname.endswith("xlsx"):
                    t1_q_sheets = pd.ExcelFile(t1_query_file).sheet_names
                    t1_q_sheet = st.selectbox("Query sheet", t1_q_sheets, key="t1_q_sheet")
                    t1_query_file.seek(0)
                    t1_qdf = pd.read_excel(t1_query_file, sheet_name=t1_q_sheet)
                else:
                    t1_qdf = pd.read_csv(t1_query_file)
                st.dataframe(t1_qdf.head(5), use_container_width=True)
                t1_default_qcols = [
                    c for c in t1_qdf.columns
                    if any(k in str(c).lower() for k in ["seq", "5'", "tag", "reporter", "id"])
                ]
                t1_q_cols = st.multiselect(
                    "Column(s) to use as query",
                    t1_qdf.columns,
                    default=t1_default_qcols,
                    key="t1_q_cols",
                )
                if t1_q_cols:
                    for c in t1_q_cols:
                        t1_query_seqs += t1_qdf[c].dropna().astype(str).str.strip().tolist()
                    st.caption(f"{len(t1_query_seqs):,} query values loaded")
    else:
        t1_pasted = st.text_area(
            "Paste sequences (one per line)",
            height=150,
            key="t1_paste",
        )
        t1_query_seqs = [l.strip() for l in t1_pasted.splitlines() if l.strip()]
        if t1_query_seqs:
            st.caption(f"{len(t1_query_seqs):,} sequences loaded")

    # ── Step 2: Reference files ────────────────────────────────────────────────
    st.markdown("### Step 2 — Reference File(s)")
    t1_ref_files = st.file_uploader(
        "Upload one or more reference files",
        type=["xlsx", "csv"],
        accept_multiple_files=True,
        key="t1_ref_files",
    )

    # Per-file config: sheet, preview, column selectors
    t1_ref_configs = []   # list of dicts: {df, label, lookup_col, extra_cols}
    for i, rf in enumerate(t1_ref_files):
        st.markdown(f"**Reference file {i+1}: {rf.name}**")
        rf_name = rf.name.lower()
        if rf_name.endswith("xlsx"):
            rf_sheets = pd.ExcelFile(rf).sheet_names
            default_sheet_idx = (
                rf_sheets.index("Sequence List") if "Sequence List" in rf_sheets else 0
            )
            rf_sheet = st.selectbox(
                f"Sheet — {rf.name}",
                rf_sheets,
                index=default_sheet_idx,
                key=f"t1_rf_sheet_{i}",
            )
            rf.seek(0)
            rf_df = pd.read_excel(rf, sheet_name=rf_sheet)
        else:
            rf_df = pd.read_csv(rf)

        # Immediate preview
        st.dataframe(rf_df.head(5), use_container_width=True)

        rf_col1, rf_col2 = st.columns(2)
        with rf_col1:
            # Smart default: prefer sequence/tag columns
            def _best_col(df, hints):
                for h in hints:
                    for c in df.columns:
                        if h in str(c).lower():
                            return list(df.columns).index(c)
                return 0
            lookup_idx = _best_col(rf_df, ["seq", "5'", "reporter", "tag", "target"])
            t1_lookup_col = st.selectbox(
                f"Search / lookup in this column",
                rf_df.columns,
                index=lookup_idx,
                key=f"t1_lookup_col_{i}",
            )
        with rf_col2:
            t1_extra_cols = st.multiselect(
                f"Also return these columns (optional)",
                [c for c in rf_df.columns if c != t1_lookup_col],
                default=[],
                key=f"t1_extra_cols_{i}",
            )
        t1_ref_configs.append({
            "df": rf_df,
            "name": rf.name,
            "lookup_col": t1_lookup_col,
            "extra_cols": t1_extra_cols,
        })

    # ── Step 3: Operation ──────────────────────────────────────────────────────
    st.markdown("### Step 3 — Operation")
    t1_mode = st.radio(
        "What would you like to do?",
        [
            "🔍 Find — search for my query within the reference column",
            "🏷️ Map — look up a paired value for each query (e.g. tag ↔ sequence)",
        ],
        key="t1_mode",
    )

    if t1_mode.startswith("🔍"):
        # Find mode options
        t1_col_a, t1_col_b = st.columns(2)
        with t1_col_a:
            t1_match_type = st.radio(
                "Match type",
                ["Exact", "Partial", "Fuzzy"],
                key="t1_match_type",
                horizontal=True,
            )
        with t1_col_b:
            t1_allow_rc = st.checkbox("Match reverse complement", value=True, key="t1_rc")
        if t1_match_type == "Fuzzy":
            t1_mismatches = st.slider("Max mismatches (fuzzy)", 1, 5, 1, key="t1_mismatches")
        else:
            t1_mismatches = 0
        t1_map_return_col = None
        t1_only_matches = False
    else:
        # Map mode options
        t1_match_type = "Exact"
        t1_allow_rc = False
        t1_mismatches = 0
        # Return column selector — shown per reference file
        st.markdown("**Map settings**")
        # We collect these per ref file below at run time; show one global note
        st.caption(
            "For each reference file, select the column whose value you want returned for each match. "
            "This is set per file in Step 2 via 'Also return these columns', or use the selector below."
        )
        t1_only_matches = st.checkbox(
            "Show only successful matches (hide NOT FOUND)",
            value=False,
            key="t1_only_matches",
        )
        # Per-ref map target column
        for i, cfg in enumerate(t1_ref_configs):
            cfg["map_col"] = st.selectbox(
                f"Return value from this column — {cfg['name']}",
                cfg["df"].columns,
                index=min(1, len(cfg["df"].columns) - 1),
                key=f"t1_map_col_{i}",
            )

    t1_simple = st.checkbox(
        "Simple output: Input + Result only (no extra columns)",
        value=True,
        key="t1_simple",
    )

    # ── Step 4: Run ───────────────────────────────────────────────────────────
    st.markdown("### Step 4 — Run")
    if st.button("▶ Run", key="t1_run", type="primary"):
        if not t1_query_seqs:
            st.error("No query sequences loaded. Check Step 1.")
            st.stop()
        if not t1_ref_configs:
            st.error("No reference files uploaded. Check Step 2.")
            st.stop()

        all_results = []

        if t1_mode.startswith("🔍"):
            # ── Find mode ──
            for cfg in t1_ref_configs:
                rf_df      = cfg["df"]
                rf_name    = cfg["name"]
                lookup_col = cfg["lookup_col"]
                extra_cols = cfg["extra_cols"]

                for idx, row in rf_df.iterrows():
                    ref_val = str(row[lookup_col]).strip().upper()
                    ref_vals = [ref_val]
                    if t1_allow_rc:
                        ref_vals.append(reverse_complement(ref_val))

                    for q_raw in t1_query_seqs:
                        q = q_raw.strip().upper()
                        if not q:
                            continue
                        for rv in ref_vals:
                            matched = False
                            if t1_match_type == "Exact":
                                matched = q == rv
                            elif t1_match_type == "Partial":
                                matched = q in rv
                            elif t1_match_type == "Fuzzy":
                                matched = is_fuzzy_match(q, rv, t1_mismatches)
                            if matched:
                                highlighted = highlight_match(rv, q)
                                result_row = {
                                    "Query":      q_raw,
                                    "Match":      rv,
                                    "Highlighted": highlighted,
                                    "File":       rf_name,
                                    "Row":        idx,
                                }
                                if not t1_simple:
                                    for c in extra_cols:
                                        result_row[c] = row[c]
                                all_results.append(result_row)
                                break  # one match per ref row per query

        else:
            # ── Map mode ──
            for cfg in t1_ref_configs:
                rf_df      = cfg["df"]
                rf_name    = cfg["name"]
                lookup_col = cfg["lookup_col"]
                map_col    = cfg.get("map_col", lookup_col)
                extra_cols = cfg["extra_cols"]

                # Build lookup dict (case-insensitive, strip asterisks)
                t1_lookup = {}
                for _, row in rf_df.iterrows():
                    k = str(row[lookup_col]).strip().upper().replace("*", "")
                    if k and k != "NAN" and k not in t1_lookup:
                        entry = {"_result": str(row[map_col]).strip()}
                        for c in extra_cols:
                            entry[c] = row[c]
                        t1_lookup[k] = entry

                for q_raw in t1_query_seqs:
                    q_clean = q_raw.strip().upper().replace("*", "")
                    if not q_clean or q_clean == "NAN":
                        continue
                    match = t1_lookup.get(q_clean)
                    if match:
                        result_row = {
                            "Input":  q_raw,
                            "Result": match["_result"],
                            "File":   rf_name,
                        }
                        if not t1_simple:
                            for c in extra_cols:
                                result_row[c] = match.get(c, "")
                        all_results.append(result_row)
                    else:
                        all_results.append({
                            "Input":  q_raw,
                            "Result": "NOT FOUND",
                            "File":   rf_name,
                        })

            if t1_only_matches:
                all_results = [r for r in all_results if r.get("Result") != "NOT FOUND"]

        if not all_results:
            st.warning("No results found. Check your column selections and query values.")
            st.stop()

        t1_result_df = pd.DataFrame(all_results)

        # Summary
        if t1_mode.startswith("🔍"):
            st.success(f"**{len(t1_result_df):,}** matches found across {len(t1_ref_configs)} file(s).")
        else:
            total = len(t1_query_seqs) * len(t1_ref_configs)
            found = int((t1_result_df["Result"] != "NOT FOUND").sum())
            st.success(f"**{found:,} / {total:,}** items matched.")

        # Results table
        st.dataframe(t1_result_df, use_container_width=True)

        # Courier New rendering if sequences present
        seq_cols = [c for c in t1_result_df.columns if any(
            k in str(c).lower() for k in ["seq", "match", "result", "5'"]
        )]
        if seq_cols:
            with st.expander("View sequences in Courier New", expanded=False):
                display_cols = ["Query" if "Query" in t1_result_df.columns else "Input"] + seq_cols
                display_cols = [c for c in display_cols if c in t1_result_df.columns]
                seq_html = t1_result_df[display_cols].to_html(index=False, escape=False)
                seq_html = seq_html.replace(
                    "<td>",
                    "<td style='font-family: Courier New; font-size: 11pt;'>"
                )
                st.markdown(seq_html, unsafe_allow_html=True)

        # Downloads
        import io as _io
        dl_col1, dl_col2 = st.columns(2)
        with dl_col1:
            t1_xlsx_buf = _io.BytesIO()
            with pd.ExcelWriter(t1_xlsx_buf, engine="openpyxl") as writer:
                t1_result_df.to_excel(writer, index=False)
            t1_xlsx_buf.seek(0)
            st.download_button(
                "⬇️ Download as .xlsx",
                data=t1_xlsx_buf,
                file_name="search_map_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="t1_dl_xlsx",
            )
        with dl_col2:
            st.download_button(
                "⬇️ Download as .csv",
                data=t1_result_df.to_csv(index=False).encode("utf-8"),
                file_name="search_map_results.csv",
                mime="text/csv",
                key="t1_dl_csv",
            )


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

# ============================================================================
# Tab 3: DNA & RNA Folding (added in v6)
# ----------------------------------------------------------------------------
# Standalone secondary-structure folder using ViennaRNA, with correct DNA
# energy parameters (Mathews 2004/1999) and Owczarzy 2008 monovalent-salt
# correction for Mg2+: [Na+]_eff = [Na+] + 120 * sqrt([Mg2+]).
#
# This tab is fully self-contained and adds NO dependencies for any other
# tab. If ViennaRNA isn't installed, the tab degrades gracefully with a
# message instead of crashing the whole app.
# ============================================================================

_FOLDTYPE_LOADERS_V6 = {
    "DNA4.0 (Mathews 2004)": "params_load_DNA_Mathews2004",
    "DNA3.0 (Mathews 1999)": "params_load_DNA_Mathews1999",
    "RNA4.0 (Turner 2004)":  "params_load_RNA_Turner2004",
    "RNA3.0 (Turner 1999)":  "params_load_RNA_Turner1999",
}


def _effective_na_mM_v6(na_mM, mg_mM):
    """Owczarzy 2008: [Na+]_eff = [Na+] + 120 * sqrt([Mg2+]). All in mM."""
    try:
        return float(na_mM) + 120.0 * _sqrt(max(0.0, float(mg_mM)))
    except (TypeError, ValueError):
        return float(na_mM) if na_mM is not None else 137.0


def _fold_v6(seq, foldtype_label, temperature, dangles, no_gu, na_mM, mg_mM):
    """Fold a single sequence. Returns (processed_seq, dot_bracket, mfe, eff_na_mM).

    DNA params are loaded for DNA* foldtypes; T is kept as T.
    RNA params load for RNA* foldtypes; T is converted to U before folding.
    """
    loader_name = _FOLDTYPE_LOADERS_V6.get(foldtype_label, "params_load_RNA_Turner2004")
    loader = getattr(_RNA, loader_name, None)
    if loader is not None:
        loader()
    is_dna_params = foldtype_label.startswith("DNA")
    # Normalize input case + whitespace
    raw = "".join(seq.split()).upper()
    if is_dna_params:
        processed = raw  # keep T as-is; DNA params accept T natively
    else:
        processed = raw.replace("T", "U")
    md = _RNA.md()
    md.temperature = float(temperature)
    md.dangles = int(dangles)
    md.noGU = 1 if no_gu else 0
    eff_na_mM = _effective_na_mM_v6(na_mM, mg_mM)
    md.salt = eff_na_mM / 1000.0  # ViennaRNA expects molar
    fc = _RNA.fold_compound(processed, md)
    structure, mfe = fc.mfe()
    return processed, structure, mfe, eff_na_mM


def _arc_diagram_v6(seq, structure, title=""):
    """Render an arc diagram: sequence along x-axis, semicircle arcs for base pairs."""
    fig, ax = plt.subplots(figsize=(max(6, len(seq) * 0.25), 3.2))
    # Build pair list from dot-bracket
    stack = []
    pairs = []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    # Draw arcs
    max_h = 0
    for (i, j) in pairs:
        cx = (i + j) / 2.0
        r = (j - i) / 2.0
        theta = np.linspace(0, np.pi, 60)
        x = cx + r * np.cos(theta)
        y = r * np.sin(theta)
        max_h = max(max_h, r)
        ax.plot(x, y, color="#3a6ea5", linewidth=1.0, alpha=0.85)
    # Draw sequence
    for i, nt in enumerate(seq):
        color = {"A": "#1f9d55", "T": "#e3342f", "U": "#e3342f",
                 "G": "#f6993f", "C": "#3490dc"}.get(nt, "black")
        ax.text(i, -0.6, nt, ha="center", va="top",
                fontsize=9, family="monospace", color=color, weight="bold")
    ax.set_xlim(-1, len(seq))
    ax.set_ylim(-2, max(2, max_h + 1))
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=10, family="monospace")
    plt.tight_layout()
    return fig


with tab3:
    with st.expander("About This Tab"):
        st.write(
            "Fold one or more DNA or RNA sequences using ViennaRNA. "
            "DNA folding uses the Mathews 2004/1999 energy parameters "
            "(NOT the default RNA Turner parameters), and Mg²⁺ is folded "
            "into an effective monovalent [Na⁺] via the Owczarzy 2008 "
            "correction: [Na⁺]_eff = [Na⁺] + 120·√[Mg²⁺]. "
            "Useful when the DNA mfold server is down or for quick "
            "interactive folds outside the SELEX pipeline."
        )

    st.header("DNA & RNA Folding")

    if not VIENNARNA_AVAILABLE:
        st.error(
            "ViennaRNA Python bindings are not installed in this environment. "
            "Install with `pip install ViennaRNA` (or the appropriate conda "
            "package) and restart the app to use this tab."
        )
    else:
        # ---- Folding parameter controls ----
        st.subheader("Folding Settings")
        col_a, col_b, col_c = st.columns(3)
        with col_a:
            foldtype_v6 = st.selectbox(
                "Energy Rules",
                list(_FOLDTYPE_LOADERS_V6.keys()),
                index=0,
                key="v6_foldtype",
                help="DNA energy rules for DNA sequences; RNA energy rules for RNA sequences."
            )
            temp_v6 = st.number_input(
                "Folding Temperature (°C)",
                min_value=0.0, max_value=100.0, value=25.0, step=0.5,
                key="v6_temp"
            )
        with col_b:
            dangles_label_v6 = st.selectbox(
                "Dangling Ends",
                ["None (0)", "Some (1)", "All (2)"],
                index=2,
                key="v6_dangles"
            )
            dangles_v6 = int(dangles_label_v6.split("(")[1].split(")")[0])
            _wobble_default_v6 = not foldtype_v6.startswith("DNA")
            allow_wobbles_v6 = st.checkbox(
                "Allow G-U / G-T wobble pairs",
                value=_wobble_default_v6,
                key="v6_wobbles",
                help="Off by default for DNA (G-T wobbles are weak), on for RNA."
            )
            no_gu_v6 = not allow_wobbles_v6
        with col_c:
            na_v6 = st.number_input(
                "[Na⁺] Concentration (mM)",
                min_value=0.0, max_value=1000.0, value=137.0, step=10.0,
                key="v6_na"
            )
            mg_v6 = st.number_input(
                "[Mg²⁺] Concentration (mM)",
                min_value=0.0, max_value=100.0, value=1.0, step=0.5,
                key="v6_mg"
            )
            st.caption(
                f"Effective monovalent [Na⁺] = "
                f"{_effective_na_mM_v6(na_v6, mg_v6):.1f} mM"
            )

        st.divider()

        # ---- Input ----
        st.subheader("Sequence Input")
        input_mode_v6 = st.radio(
            "Input mode",
            ["Single sequence", "Multiple sequences (one per line, optional FASTA)"],
            key="v6_input_mode"
        )

        sequences_to_fold = []  # list of (name, seq)
        if input_mode_v6 == "Single sequence":
            single_name = st.text_input("Optional name/ID", value="seq1", key="v6_single_name")
            single_seq = st.text_area(
                "Sequence (DNA or RNA; A/C/G/T/U)",
                value="",
                height=100,
                key="v6_single_seq",
                help="Whitespace and case are ignored."
            )
            if single_seq.strip():
                sequences_to_fold = [(single_name.strip() or "seq1", single_seq)]
        else:
            multi_text = st.text_area(
                "Paste sequences (one per line, or FASTA format)",
                value="",
                height=180,
                key="v6_multi_seq",
                help="Plain: one sequence per line. FASTA: lines starting with '>' are names."
            )
            if multi_text.strip():
                # Parse plain or FASTA
                lines = [ln.rstrip() for ln in multi_text.splitlines() if ln.strip()]
                if any(ln.startswith(">") for ln in lines):
                    cur_name, cur_seq_parts = None, []
                    for ln in lines:
                        if ln.startswith(">"):
                            if cur_name is not None:
                                sequences_to_fold.append((cur_name, "".join(cur_seq_parts)))
                            cur_name = ln[1:].strip() or f"seq{len(sequences_to_fold)+1}"
                            cur_seq_parts = []
                        else:
                            cur_seq_parts.append(ln)
                    if cur_name is not None:
                        sequences_to_fold.append((cur_name, "".join(cur_seq_parts)))
                else:
                    for i, ln in enumerate(lines, 1):
                        sequences_to_fold.append((f"seq{i}", ln))

        show_arc_v6 = st.checkbox("Show arc diagram for each sequence", value=True, key="v6_show_arc")
        run_fold_v6 = st.button("🧬 Fold", key="v6_run_fold")

        # ---- Run ----
        if run_fold_v6:
            if not sequences_to_fold:
                st.warning("Please enter at least one sequence.")
            else:
                results = []
                valid_alphabet = set("ACGTU")
                with st.spinner(f"Folding {len(sequences_to_fold)} sequence(s)..."):
                    for name, raw_seq in sequences_to_fold:
                        cleaned = "".join(raw_seq.split()).upper()
                        if not cleaned:
                            results.append({
                                "Name": name, "Sequence": "", "Length": 0,
                                "Dot-Bracket": "", "MFE (kcal/mol)": None,
                                "Error": "empty sequence"
                            })
                            continue
                        bad = set(cleaned) - valid_alphabet
                        if bad:
                            results.append({
                                "Name": name, "Sequence": cleaned,
                                "Length": len(cleaned),
                                "Dot-Bracket": "", "MFE (kcal/mol)": None,
                                "Error": f"invalid characters: {','.join(sorted(bad))}"
                            })
                            continue
                        try:
                            proc, db, mfe, eff_na = _fold_v6(
                                cleaned, foldtype_v6, temp_v6, dangles_v6,
                                no_gu_v6, na_v6, mg_v6
                            )
                            results.append({
                                "Name": name,
                                "Sequence": proc,
                                "Length": len(proc),
                                "Dot-Bracket": db,
                                "MFE (kcal/mol)": round(float(mfe), 3),
                                "Error": ""
                            })
                        except Exception as ex:
                            results.append({
                                "Name": name, "Sequence": cleaned,
                                "Length": len(cleaned),
                                "Dot-Bracket": "", "MFE (kcal/mol)": None,
                                "Error": f"fold failed: {ex}"
                            })

                # Methods summary
                eff_na_display = _effective_na_mM_v6(na_v6, mg_v6)
                _wobble_text = "disabled" if no_gu_v6 else "enabled"
                st.info(
                    f"Folded with ViennaRNA — Energy rules: **{foldtype_v6}**, "
                    f"T = {temp_v6}°C, [Na⁺] = {na_v6} mM, [Mg²⁺] = {mg_v6} mM "
                    f"(effective [Na⁺] = {eff_na_display:.1f} mM via Owczarzy 2008), "
                    f"Dangling ends = {dangles_v6}, "
                    f"G-U/G-T wobbles {_wobble_text}."
                )

                # Per-sequence display
                for r in results:
                    st.markdown(f"### `{r['Name']}`")
                    if r["Error"]:
                        st.error(r["Error"])
                        continue
                    mfe_str = (
                        f"{r['MFE (kcal/mol)']:.2f} kcal/mol"
                        if r["MFE (kcal/mol)"] is not None else "n/a"
                    )
                    st.markdown(
                        f"**Length:** {r['Length']} nt  **MFE:** {mfe_str}"
                    )
                    st.code(
                        f"Seq: {r['Sequence']}\nDB:  {r['Dot-Bracket']}",
                        language="text"
                    )
                    if show_arc_v6 and r["Dot-Bracket"]:
                        try:
                            fig = _arc_diagram_v6(
                                r["Sequence"], r["Dot-Bracket"],
                                title=f"{r['Name']}  (MFE = {mfe_str})"
                            )
                            st.pyplot(fig)
                            plt.close(fig)
                        except Exception as ex:
                            st.warning(f"Could not render arc diagram: {ex}")

                # Combined results table + download
                df_results = pd.DataFrame(results)
                st.markdown("### Results Table")
                st.dataframe(
                    df_results.style.set_properties(**{"font-family": "Courier New"}),
                    use_container_width=True
                )
                st.download_button(
                    "Download results as CSV",
                    df_results.to_csv(index=False).encode("utf-8"),
                    file_name="folding_results.csv",
                    mime="text/csv",
                    key="v6_download"
                )
# Tab 4: Melting Temp Calculator
with tab4:
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

# Tab 5: Venn Diagrams with Sequence Extraction
with tab5:
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
