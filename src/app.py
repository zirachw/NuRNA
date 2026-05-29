import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st

from NuRNA import (
    ANIMATE_MAX_LEN,
    build_dp_matrix, build_dp_matrix_steps,
    traceback as nussinov_traceback,
    validate_sequence, parse_dot_bracket, gc_content,
    plot_arc_diagram, plot_dp_matrix, plot_rna_2d, plot_comparison,
    compute_metrics,
)

st.set_page_config(page_title="NuRNA", layout="wide", initial_sidebar_state="expanded")

_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Instrument+Serif&family=Instrument+Sans:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500&display=swap');

/* ── App shell ── */
[data-testid="stAppViewContainer"] { background-color: #F7F3EE !important; }
[data-testid="stMainBlockContainer"] {
    padding: 1rem 4rem 0 4rem !important;
}
[data-testid="stBottom"],
[data-testid="stToolbar"]  { display: none !important; }
[data-testid="stHeader"]   {
    background: transparent !important; border: none !important;
    position: fixed !important; top: 0 !important;
    z-index: 9999 !important; pointer-events: none !important;
}
[data-testid="stHeader"] button { pointer-events: auto !important; }

/* ── Global typography ── */
h1, h2, h3 { font-family: 'Instrument Serif', serif !important; }
p, label, div[data-testid="stMarkdown"] p { font-family: 'Instrument Sans', sans-serif !important; }
[data-testid="stTextArea"] textarea {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 12px !important;
    background-color: #F0EBE3 !important;
}

/* ── Global buttons ── */
[data-testid="stButton"] button[kind="primary"] {
    background-color: #2D5A4E !important;
    color: white !important; border: none !important;
    font-family: 'Instrument Sans', sans-serif !important;
}
[data-testid="stButton"] button[kind="primary"]:hover { background-color: #1F3D35 !important; }
[data-testid="stButton"] button[kind="secondary"] {
    background-color: #EDE8E0 !important; color: #1C1C1C !important;
    border: 1px solid #C8C0B4 !important;
    font-family: 'Instrument Sans', sans-serif !important;
}
[data-testid="stButton"] button[kind="secondary"]:hover { background-color: #DDD6CC !important; }

/* ── Tabs (pill style) ── */
[data-baseweb="tab-list"]  { gap: 6px !important; background: transparent !important; border-bottom: none !important; }
[data-baseweb="tab"]       {
    font-family: 'Instrument Sans', sans-serif !important; font-size: 13px !important;
    background: #EDE8E0 !important; border: 1px solid #C8C0B4 !important;
    border-radius: 999px !important; padding: 4px 14px !important; color: #444 !important;
}
[aria-selected="true"][data-baseweb="tab"] { background: #2D5A4E !important; border-color: #2D5A4E !important; color: white !important; }
[data-baseweb="tab-highlight"],
[data-baseweb="tab-border"] { display: none !important; }

/* ── Metric cards ── */
.metric-card  { background: #EDE8E0; border-radius: 8px; padding: clamp(8px,1vw,13px) clamp(10px,1.2vw,15px); margin-bottom: 8px; }
.metric-label { font-family: 'Instrument Sans', sans-serif; font-size: clamp(10px,1.1vw,13px); font-weight: 500; color: #1C1C1C; margin: 0 0 3px 0; }
.metric-value { font-family: 'Instrument Sans', sans-serif; font-size: clamp(9px,0.9vw,11px); color: #666; margin: 0; line-height: 1.3; }
.col-heading     { font-family: 'Instrument Sans', sans-serif; font-weight: 500; font-size: 14px; color: #1C1C1C; margin: 0 0 4px 0; padding-top: 44px; }
.col-heading-gap { padding-top: 16px; margin: 0 0 4px 0; }

/* ── Welcome screen ── */
.welcome-center { min-height: calc(100vh - 4rem); display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; }

/* ════════ SIDEBAR ════════ */
[data-testid="stSidebarHeader"],
[data-testid="stSidebarCollapseButton"] { display: none !important; }
section[data-testid="stSidebar"] [style*="col-resize"] { display: none !important; pointer-events: none !important; }
section[data-testid="stSidebar"] { height: 100vh !important; }
[data-testid="stSidebarContent"] {
    height: 100% !important;
    display: flex !important;
    flex-direction: column !important;
}
[data-testid="stSidebarUserContent"] {
    flex: 1 !important;
    display: flex !important;
    flex-direction: column !important;
    justify-content: center !important;
    padding-bottom: 0 !important;
}

[data-testid="stSidebar"] [data-testid="stCaptionContainer"] p { font-size: 10.5px !important; line-height: 1.3 !important; }

/* Tooltip */
[data-testid="stTooltipContent"],
[data-testid="stTooltipContent"] p { font-family: 'Instrument Sans', sans-serif !important; font-size: 11px !important; }

/* Toggle label */
[data-testid="stSidebar"] [data-testid="stCheckbox"] p {
    font-family: 'Instrument Sans', sans-serif !important;
    font-size: 12px !important;
}

/* Predict button */
[data-testid="stSidebar"] [data-testid="stBaseButton-primary"] {
    font-family: 'Instrument Sans', sans-serif !important;
    font-size: 13px !important;
    font-weight: 500 !important;
    border-radius: 8px !important;
    padding: 8px 12px !important;
}
[data-testid="stSidebar"] [data-testid="stBaseButton-primary"] p {
    font-family: 'Instrument Sans', sans-serif !important;
    font-size: 13px !important;
    font-weight: 500 !important;
}
[data-testid="stSidebar"] button[data-testid="stBaseButton-primary"][disabled] {
    background-color: #7A9E95 !important;
    cursor: not-allowed !important;
}

[data-testid="stSidebar"] .sb-heading { font-family: 'Instrument Sans', sans-serif !important; font-size: 11px !important; font-weight: 600 !important; color: #1C1C1C !important; margin: 0 !important; padding: 0 0 5px 0 !important; letter-spacing: 0.04em; text-transform: uppercase; }

/* Examples container gap */
[data-testid="stSidebar"] [data-testid="stLayoutWrapper"] [data-testid="stVerticalBlock"] {
    gap: 0.5rem !important;
}

/* Example buttons */
[data-testid="stSidebar"] [data-testid="stBaseButton-secondary"] {
    font-size: 11px !important;
    padding: 3px 8px !important;
    border: none !important;
    border-radius: 6px !important;
    text-align: left !important;
    background-color: #E4DDD3 !important;
}
[data-testid="stSidebar"] [data-testid="stBaseButton-secondary"] p {
    font-size: 11px !important;
    line-height: 1.4 !important;
}
[data-testid="stSidebar"] [data-testid="stBaseButton-secondary"]:hover {
    background-color: #D0C9BF !important;
}

/* Active example button */
[data-testid="stSidebar"] [data-testid="stLayoutWrapper"] [data-testid="stBaseButton-primary"] {
    font-size: 11px !important;
    padding: 3px 8px !important;
    border: none !important;
    border-radius: 6px !important;
    text-align: left !important;
    background-color: #2D5A4E !important;
    color: white !important;
}
[data-testid="stSidebar"] [data-testid="stLayoutWrapper"] [data-testid="stBaseButton-primary"] p {
    font-size: 11px !important;
    line-height: 1.4 !important;
    color: white !important;
}
[data-testid="stSidebar"] [data-testid="stLayoutWrapper"] [data-testid="stBaseButton-primary"]:hover {
    background-color: #1F3D35 !important;
}

/* Search bar */
[data-testid="stSidebar"] [data-testid="stTextInputRootElement"] {
    background-color: #F7F3EE !important;
    border: 1px solid #C8C0B4 !important;
    border-radius: 6px !important;
}
[data-testid="stSidebar"] [data-testid="stTextInputRootElement"]:focus-within {
    border: 1px solid #2D5A4E !important;
}
[data-testid="stSidebar"] [data-testid="stTextInputRootElement"] input {
    font-family: 'Instrument Sans', sans-serif !important;
    font-size: 12px !important;
    padding: 5px 8px 5px 10px !important;
    color: #1C1C1C !important;
}
</style>
"""
st.markdown(_CSS, unsafe_allow_html=True)


@st.cache_data
def _load_examples() -> list[dict]:
    path = Path(__file__).parent / "examples.json"
    with open(path) as f:
        return json.load(f)


EXAMPLES = _load_examples()

for _k, _v in {
    "seq_box": "", "ref_box": "", "result": None, "active_example": None,
    "is_animating": False, "anim_seq": "", "anim_speed": 1.0,
}.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("# NuRNA")
    st.caption("Nussinov RNA Secondary Structure Prediction")

    st.markdown("<p class='sb-heading'>Examples</p>", unsafe_allow_html=True)
    _search = st.text_input(
        "Search",
        placeholder="Search by name...",
        label_visibility="collapsed",
    )
    _filtered = [
        ex for ex in EXAMPLES
        if not _search or _search.lower() in f"{ex['family']}: {ex['label']}".lower()
    ]
    with st.container(height=90, border=False):
        for ex in _filtered:
            _is_active = st.session_state.active_example == ex["label"]
            if st.button(
                f"{ex['family']}: {ex['label']}",
                key=f"ex__{ex['label']}",
                use_container_width=True,
                type="primary" if _is_active else "secondary",
            ):
                st.session_state.seq_box = ex["sequence"]
                st.session_state.ref_box = ex.get("reference", "")
                st.session_state.result = None
                st.session_state.active_example = ex["label"]
                st.rerun()

    st.markdown("<p class='sb-heading'>Sequence</p>", unsafe_allow_html=True)
    seq: str = st.text_area(
        "Sequence",
        key="seq_box",
        height=80,
        placeholder="AUGCAUGCAUGC...",
        help="A, U, G, C only. Case-insensitive.",
        label_visibility="collapsed",
    ).upper().strip()
    st.markdown("<p class='sb-heading'>Reference Dot-Bracket</p>", unsafe_allow_html=True)
    ref_notation: str = st.text_area(
        "Reference",
        key="ref_box",
        height=80,
        placeholder="(((...))) – optional",
        label_visibility="collapsed",
    ).strip()

    _locked    = st.session_state.is_animating
    animatable = bool(seq) and len(seq) <= ANIMATE_MAX_LEN and not _locked
    speed = 1.0

    st.markdown("<p class='sb-heading'>Settings</p>", unsafe_allow_html=True)
    animate = st.toggle(
        "Animate DP", value=False, disabled=not animatable,
        help=f"Enter a sequence of at most {ANIMATE_MAX_LEN} characters to enable this.",
    ) and animatable
    if animate:
        speed = st.slider("Speed", 0.5, 3.0, 1.0, 0.5,
                          help="Animation playback speed multiplier.")

    predict_clicked = st.button(
        "Predict", type="primary", use_container_width=True,
        disabled=_locked or not seq,
    )

# ---------------------------------------------------------------------------
# Welcome screen
# ---------------------------------------------------------------------------
if not predict_clicked and st.session_state.result is None and not st.session_state.is_animating:
    st.markdown(
        "<div class='welcome-center'>"
        "<h1 style='font-size:56px;'>NuRNA</h1>"
        "<p style='color:#666;font-size:17px;'>"
        "Nussinov dynamic programming &middot; RNA secondary structure prediction"
        "</p>"
        "<p style='color:#999;font-size:13px;margin-top:28px;'>"
        "Choose an example or enter a sequence in the sidebar, then press Predict."
        "</p>"
        "</div>",
        unsafe_allow_html=True,
    )
    st.stop()

# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------
if predict_clicked:
    if not validate_sequence(seq):
        st.error("Sequence contains invalid characters. Only A, U, G, C are accepted.")
        st.stop()

    n = len(seq)
    do_animate = animate and n <= ANIMATE_MAX_LEN

    t0 = time.perf_counter()
    dp = build_dp_matrix(seq)
    pairs: set = set()
    nussinov_traceback(dp, seq, 0, n - 1, pairs)
    elapsed = time.perf_counter() - t0

    dot_bracket = ["."] * n
    for i, j in pairs:
        dot_bracket[i] = "("
        dot_bracket[j] = ")"

    _fam = next((ex["family"] for ex in EXAMPLES if ex["label"] == st.session_state.active_example), None)
    st.session_state.result = {
        "seq": seq, "db": "".join(dot_bracket), "pairs": pairs,
        "matrix": dp if do_animate else None,
        "animated": do_animate, "elapsed": elapsed,
        "family": _fam,
    }

    if do_animate:
        st.session_state.is_animating = True
        st.session_state.anim_seq     = seq
        st.session_state.anim_speed   = speed

    st.rerun()

# ---------------------------------------------------------------------------
# Animation — blocking loop inside tab layout; result already computed
# ---------------------------------------------------------------------------
if st.session_state.is_animating:
    _seq   = st.session_state.anim_seq
    _speed = st.session_state.anim_speed
    _n     = len(_seq)
    _gc    = gc_content(_seq)

    _total_steps  = _n * (_n - 1) // 2
    _render_every = max(1, _total_steps // 60)
    _delay        = 0.04 / _speed

    _col_l, _col_r = st.columns([3, 2], gap="large")
    with _col_l:
        _tab_dp, _tab_2d, _tab_arc = st.tabs(["DP Matrix", "2D Structure", "Arc Diagram"])
        with _tab_dp:
            _placeholder = st.empty()
        with _tab_2d:
            st.info("Available after animation completes.")
        with _tab_arc:
            st.info("Available after animation completes.")

    with _col_r:
        _has_ref_anim = bool(ref_notation) and len(ref_notation) == _n
        _fam_anim = st.session_state.result.get("family")
        st.markdown("<p class='col-heading'>Stats</p>", unsafe_allow_html=True)
        st.markdown(
            f"<div class='metric-card'>"
            f"<p class='metric-label'>Family</p>"
            f"<p class='metric-value'>{_fam_anim or '—'}</p>"
            f"</div>",
            unsafe_allow_html=True,
        )
        _sa1, _sa2 = st.columns(2)
        _sb1, _sb2 = st.columns(2)
        for _col, _lbl, _val in [
            (_sa1, "Length",         f"{_n} nt"),
            (_sa2, "GC Content",     f"{_gc:.1%}"),
            (_sb1, "Base Pairs",     "..."),
            (_sb2, "Execution Time", "..."),
        ]:
            _col.markdown(
                f"<div class='metric-card'>"
                f"<p class='metric-label'>{_lbl}</p>"
                f"<p class='metric-value'>{_val}</p>"
                f"</div>",
                unsafe_allow_html=True,
            )
        if _has_ref_anim:
            st.markdown("<p class='col-heading col-heading-gap'>Metrics</p>", unsafe_allow_html=True)
            _m1, _m2 = st.columns(2)
            _m3, _m4 = st.columns(2)
            for _col, _lbl in [(_m1, "Sensitivity"), (_m2, "PPV"), (_m3, "F1"), (_m4, "MCC")]:
                _col.markdown(
                    f"<div class='metric-card'>"
                    f"<p class='metric-label'>{_lbl}</p>"
                    f"<p class='metric-value'>...</p>"
                    f"</div>",
                    unsafe_allow_html=True,
                )

    _full_matrix = st.session_state.result["matrix"]
    _vmax = max(1, max(max(row) for row in _full_matrix))

    _dp = None
    for _idx, _snap in enumerate(build_dp_matrix_steps(_seq), start=1):
        _dp = _snap
        if _idx % _render_every == 0:
            _fig = plot_dp_matrix(_snap, _seq, vmax=_vmax)
            _placeholder.pyplot(_fig, use_container_width=True)
            plt.close(_fig)
            time.sleep(_delay)

    if _dp is not None:
        _fig = plot_dp_matrix(_dp, _seq, vmax=_vmax)
        _placeholder.pyplot(_fig, use_container_width=True)
        plt.close(_fig)

    st.session_state.is_animating = False
    st.rerun()

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
result = st.session_state.result
if result is None:
    st.stop()

rseq  = result["seq"]
pairs = result["pairs"]
db    = result["db"]
gc    = gc_content(rseq)

has_ref   = bool(ref_notation) and len(ref_notation) == len(rseq)
ref_pairs = parse_dot_bracket(ref_notation) if has_ref else None
metrics   = compute_metrics(pairs, ref_pairs) if has_ref else None

show_matrix  = result["animated"] and result["matrix"] is not None
chip_options = ["DP Matrix", "2D Structure", "Arc Diagram"] if show_matrix else ["2D Structure", "Arc Diagram"]

col_left, col_right = st.columns([3, 2], gap="large")

with col_left:
    tabs    = st.tabs(chip_options)
    tab_mat = tabs[0] if show_matrix else None
    tab_2d  = tabs[1] if show_matrix else tabs[0]
    tab_arc = tabs[2] if show_matrix else tabs[1]

    with tab_arc:
        if pairs:
            fig = plot_arc_diagram(rseq, pairs)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)
        else:
            st.info("No base pairs predicted for this sequence.")

    if tab_mat is not None:
        with tab_mat:
            fig = plot_dp_matrix(result["matrix"], rseq)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    with tab_2d:
        if has_ref and ref_pairs:
            fig = plot_comparison(rseq, pairs, ref_pairs)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)
        elif pairs:
            fig = plot_rna_2d(rseq, pairs)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)
        else:
            st.info("No base pairs predicted for this sequence.")

with col_right:
    elapsed = result["elapsed"]
    elapsed_str = f"{elapsed * 1000:.1f} ms" if elapsed < 1 else f"{elapsed:.2f} s"

    fam = result.get("family")
    st.markdown("<p class='col-heading'>Stats</p>", unsafe_allow_html=True)
    st.markdown(
        f"<div class='metric-card'>"
        f"<p class='metric-label'>Family</p>"
        f"<p class='metric-value'>{fam or '—'}</p>"
        f"</div>",
        unsafe_allow_html=True,
    )
    sa1, sa2 = st.columns(2)
    sb1, sb2 = st.columns(2)
    for col, lbl, val in [
        (sa1, "Length",         f"{len(rseq)} nt"),
        (sa2, "GC Content",     f"{gc:.1%}"),
        (sb1, "Base Pairs",     str(len(pairs))),
        (sb2, "Execution Time", elapsed_str),
    ]:
        col.markdown(
            f"<div class='metric-card'>"
            f"<p class='metric-label'>{lbl}</p>"
            f"<p class='metric-value'>{val}</p>"
            f"</div>",
            unsafe_allow_html=True,
        )

    if has_ref and metrics:
        st.markdown("<p class='col-heading col-heading-gap'>Metrics</p>", unsafe_allow_html=True)
        m1, m2 = st.columns(2)
        m3, m4 = st.columns(2)
        for col, lbl, val in [
            (m1, "Sensitivity", f"{metrics['sensitivity']:.3f}"),
            (m2, "PPV",         f"{metrics['ppv']:.3f}"),
            (m3, "F1",          f"{metrics['f1']:.3f}"),
            (m4, "MCC",         f"{metrics['mcc']:.3f}"),
        ]:
            col.markdown(
                f"<div class='metric-card'>"
                f"<p class='metric-label'>{lbl}</p>"
                f"<p class='metric-value'>{val}</p>"
                f"</div>",
                unsafe_allow_html=True,
            )
