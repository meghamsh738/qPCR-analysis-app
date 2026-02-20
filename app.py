# qPCR Analysis App (Streamlit) — Paste support + auto Cq detection + replicate averages
# Run: streamlit run app.py

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import re
import json
import os
from pathlib import Path
from datetime import datetime

from qpcr_core import (
    VALID_COLS,
    DEFAULT_TOP_CONC,
    DEFAULT_DILUTION_FACTOR,
    _guess_label_col,
    _guess_rep_col,
    _best_numeric_column_as_cq,
    coerce_columns,
    fit_standard_curve,
    suggest_standard_curve_exclusions,
    quantify_samples,
    mark_outliers,
    replicate_stats,
    normalize_to_ref,
    download_excel,
    read_from_paste,
    read_from_file,
    _serial_dilution,
    make_std_curve_figure,
)

SETTINGS_PATH = Path.home() / ".easylab" / "qpcr-analysis" / "settings.json"
PATH_FIELDS = [
    ("data_path", "Data folder", "Saved calculations, cached state, and metadata."),
    ("attachments_path", "Attachments folder", "Files generated or stored with this workspace."),
    ("exports_path", "Exports folder", "CSV / Excel export destination."),
    ("sync_path", "Sync folder", "Optional sync target for backups."),
]

def default_paths():
    base = Path.home() / "Documents" / "Easylab" / "qPCR Analysis"
    return {
        "data_path": str(base / "data"),
        "attachments_path": str(base / "attachments"),
        "exports_path": str(base / "exports"),
        "sync_path": str(base / "sync"),
    }

def load_settings():
    if SETTINGS_PATH.exists():
        try:
            return json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}

def save_settings(paths):
    SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    SETTINGS_PATH.write_text(json.dumps(paths, indent=2), encoding="utf-8")

def ensure_paths(paths):
    for key, value in paths.items():
        if value:
            Path(value).mkdir(parents=True, exist_ok=True)

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_WELLS_PATH = APP_DIR / "sample-data" / "qpcr_example.csv"

st.set_page_config(page_title="qPCR Analysis", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
        @import url("https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;600;700&family=Space+Mono:wght@400;700&display=swap");

        :root{
          --font-display: "Space Grotesk", "Segoe UI", system-ui, -apple-system, Arial, sans-serif;
          --font-body: "Space Grotesk", "Segoe UI", system-ui, -apple-system, Arial, sans-serif;
          --font-mono: "Space Mono", ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
          --tracking-display: 0.08em;
          --tracking-label: 0.12em;
          --space-1: 4px;
          --space-2: 8px;
          --space-3: 12px;
          --space-4: 16px;
          --space-5: 24px;
          --space-6: 32px;
          --radius-sm: 8px;
          --radius-md: 12px;
          --radius-lg: 12px;
          --bg: #F6F2EA;
          --panel: #FFFFFF;
          --surface-2: #FFF7EC;
          --surface-3: #F2EADF;
          --border: #15151A;
          --border-2: #0F0F14;
          --muted: #575762;
          --text: #111113;
          --text-subtle: #2F2F36;
          --accent: #1F5BFF;
          --accent-strong: #1F5BFF;
          --accent-soft: rgba(31, 91, 255, 0.18);
          --focus: rgba(31, 91, 255, 0.35);
          --success: #12b76a;
          --success-soft: rgba(18, 183, 106, 0.16);
          --warning: #f59e0b;
          --warning-soft: rgba(245, 158, 11, 0.18);
          --shadow: 10px 10px 0 rgba(17, 17, 20, 0.9);
          --shadow-soft: 3px 3px 0 rgba(17, 17, 20, 0.9);
          --app-gradient: radial-gradient(1200px 600px at 8% -10%, rgba(31, 91, 255, 0.12), transparent 60%),
            radial-gradient(900px 500px at 95% 10%, rgba(0, 0, 0, 0.08), transparent 70%),
            var(--bg);
        }

        *{
          box-sizing: border-box;
        }

        html, body, [data-testid="stAppViewContainer"]{
          background: var(--app-gradient);
          color: var(--text);
          font-family: var(--font-body);
          font-weight: 400;
          letter-spacing: 0;
          color-scheme: light;
        }

        body{
          margin: 0;
          min-height: 100vh;
        }

        h1, h2, h3, h4, h5{
          font-family: var(--font-display);
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: var(--tracking-display);
        }

        h1{
          font-size: 28px;
          margin-bottom: var(--space-2);
        }

        h2{
          font-size: 22px;
        }

        h3{
          font-size: 18px;
        }

        .block-container{
          max-width: min(1460px, 96vw);
          padding-top: 28px;
          padding-bottom: 64px;
        }

        [data-testid="stSidebar"]{
          background: var(--panel);
          border-right: 2px solid var(--border);
          color: var(--text);
        }

        [data-testid="stSidebar"] h1,
        [data-testid="stSidebar"] h2,
        [data-testid="stSidebar"] h3,
        [data-testid="stSidebar"] h4{
          color: var(--text);
        }

        [data-testid="stSidebar"] [data-testid="stWidgetLabel"] > label{
          color: var(--text);
          font-family: var(--font-mono);
          font-size: 12px;
          font-weight: 600;
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
        }

        [data-testid="stSidebar"] .stCaption,
        [data-testid="stSidebar"] .stMarkdown,
        [data-testid="stSidebar"] .stMarkdown p,
        [data-testid="stSidebar"] .stMarkdown span,
        [data-testid="stSidebar"] .stMarkdown li{
          color: var(--text-subtle);
        }

        [data-testid="stSidebar"] .stSlider > div > div > span{
          color: var(--text);
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"]{
          display: grid;
          gap: 6px;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label{
          border: 2px solid var(--border);
          background: var(--surface-2);
          border-radius: var(--radius-md);
          padding: 6px 10px;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label span{
          color: var(--text);
          font-weight: 600;
          opacity: 1;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label:has(input:checked){
          border-color: rgba(31, 91, 255, 0.45);
          background: var(--accent-soft);
          box-shadow: 0 0 0 2px var(--focus);
        }

        [data-testid="stSidebar"] input::placeholder,
        [data-testid="stSidebar"] textarea::placeholder{
          color: var(--muted);
        }

        .streamlit-expanderHeader{
          color: var(--text);
          font-weight: 600;
        }

        [data-testid="stMetric"]{
          background: var(--panel);
          border: 2px solid var(--border);
          border-radius: var(--radius-md);
          padding: 12px 14px;
          box-shadow: var(--shadow-soft);
        }

        [data-testid="stMetricLabel"]{
          color: var(--muted);
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          font-weight: 600;
        }

        [data-testid="stMetricValue"]{
          color: var(--text);
        }

        .stDataFrame,
        .stDataEditor{
          border-radius: var(--radius-md);
          overflow: hidden;
          border: 2px solid var(--border);
          box-shadow: var(--shadow-soft);
        }

        a{
          color: inherit;
          text-decoration: underline;
          text-underline-offset: 0.18em;
          text-decoration-color: var(--border);
        }

        a:hover{
          text-decoration-color: var(--accent);
        }

        code, pre{
          font-family: var(--font-mono);
        }

        .stButton > button,
        .stDownloadButton > button{
          height: 40px;
          border-radius: var(--radius-sm);
          border: 2px solid var(--border);
          background: var(--surface);
          color: var(--text);
          font-weight: 600;
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          transition: background-color 120ms ease-out, border-color 120ms ease-out, box-shadow 120ms ease-out, transform 120ms ease-out, filter 120ms ease-out;
          max-width: 100%;
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }

        .stButton > button:hover,
        .stDownloadButton > button:hover{
          background: var(--surface-2);
        }

        .stButton > button:active,
        .stDownloadButton > button:active{
          transform: translateY(1px);
        }

        .stButton > button:focus-visible,
        .stDownloadButton > button:focus-visible{
          outline: none;
          border-color: rgba(31, 91, 255, 0.6);
          box-shadow: 0 0 0 2px var(--focus);
        }

        button[data-testid="baseButton-primary"],
        button[kind="primary"]{
          background: var(--accent) !important;
          color: #fff !important;
          border-color: var(--border) !important;
        }

        button[data-testid="baseButton-primary"]:hover,
        button[kind="primary"]:hover{
          filter: brightness(0.96);
        }

        input, textarea, select{
          border-radius: var(--radius-sm) !important;
          border: 2px solid var(--border) !important;
          background: var(--surface-2) !important;
          color: var(--text) !important;
        }

        input:focus, textarea:focus, select:focus{
          border-color: rgba(31, 91, 255, 0.6) !important;
          box-shadow: 0 0 0 2px var(--focus) !important;
        }

        .hero{
          padding: 16px 20px;
          border-radius: var(--radius-md);
          background: var(--panel);
          color: var(--text);
          border: 2px solid var(--border);
          margin-bottom: 12px;
          box-shadow: var(--shadow-soft);
        }

        .hero h4{
          margin: 6px 0 6px 0;
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: var(--tracking-display);
        }

        .hero p{
          margin: 0;
          color: var(--text-subtle);
        }

        .pill{
          display: inline-flex;
          align-items: center;
          gap: 6px;
          padding: 4px 10px;
          border-radius: 999px;
          background: var(--accent-soft);
          border: 2px solid rgba(31, 91, 255, 0.45);
          color: var(--accent);
          font-weight: 600;
          font-size: 12px;
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          min-height: 32px;
        }

        .stMarkdown,
        .stCaption,
        [data-testid="stMetricValue"],
        [data-testid="stMetricLabel"]{
          overflow-wrap: anywhere;
        }

        ::selection{
          background: var(--accent-soft);
        }

        @keyframes fadeUp {
          from { opacity: 0; transform: translateY(6px); }
          to { opacity: 1; transform: translateY(0); }
        }

        @keyframes slideInLeft {
          from { opacity: 0; transform: translateX(-10px); }
          to { opacity: 1; transform: translateX(0); }
        }

        [data-testid="stSidebar"]{
          animation: slideInLeft 240ms ease-out;
        }

        div[data-testid="stAlert"]{
          animation: fadeUp 240ms ease-out;
        }

        .setup-card{
          background: var(--panel);
          border: 2px solid var(--border);
          border-radius: var(--radius-lg);
          padding: 24px;
          box-shadow: var(--shadow);
          animation: fadeUp 240ms ease-out;
          margin-bottom: 18px;
        }

        .setup-grid{
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 12px;
        }

        .setup-message{
          border: 2px dashed var(--border);
          border-radius: var(--radius-md);
          padding: 12px 14px;
          background: var(--surface-2);
          color: var(--text-subtle);
          margin-top: 12px;
        }

        .setup-message.error{
          border-style: solid;
          border-color: rgba(245, 158, 11, 0.45);
          background: var(--warning-soft);
          color: var(--text);
        }

        .setup-actions{
          display: flex;
          gap: 10px;
          justify-content: flex-end;
          flex-wrap: wrap;
          margin-top: 16px;
        }

        .tutorial-panel{
          border: 2px solid var(--border);
          border-radius: var(--radius-md);
          background: var(--surface-2);
          padding: 10px 12px;
          margin-top: 6px;
        }

        .tutorial-top-card{
          border: 2px solid var(--border);
          border-radius: var(--radius-md);
          background: var(--surface-2);
          padding: 10px 12px;
          margin-top: 0;
          box-shadow: var(--shadow-soft);
          position: fixed;
          right: 14px;
          bottom: 88px;
          width: min(420px, calc(100vw - 28px));
          max-width: 420px;
          max-height: 42vh;
          overflow: auto;
          z-index: 2147483650;
        }

        /* Keep Back/Next/Skip controls floating and always reachable. */
        .tutorial-top-card + div[data-testid="stHorizontalBlock"]{
          position: fixed;
          right: 14px;
          bottom: 14px;
          width: min(420px, calc(100vw - 28px));
          max-width: 420px;
          z-index: 2147483651;
          background: color-mix(in srgb, var(--surface) 88%, transparent);
          border: 1px solid var(--border);
          border-radius: var(--radius-md);
          padding: 8px;
          margin-top: 0;
          box-shadow: var(--shadow-soft);
        }

        .tutorial-top-card .kicker{
          font-family: var(--font-mono);
          text-transform: uppercase;
          letter-spacing: var(--tracking-label);
          font-size: 12px;
          color: var(--muted);
          margin-bottom: 6px;
        }

        .tutorial-panel strong{
          display: block;
          margin-bottom: 6px;
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          font-family: var(--font-mono);
          font-size: 12px;
        }

        .tutorial-callout{
          border: 3px solid rgba(31, 91, 255, 0.8);
          border-radius: var(--radius-md);
          background: var(--accent-soft);
          color: var(--text);
          padding: 10px 12px;
          margin: 6px 0 10px 0;
          box-shadow: 0 0 0 4px rgba(31, 91, 255, 0.22), var(--shadow-soft);
        }

        .tutorial-callout p{
          margin: 0;
        }

        .tutorial-anchor{
          width: 0;
          height: 0;
          opacity: 0;
          pointer-events: none;
        }

        .easylab-tutorial-focus{
          outline: 3px solid rgba(31, 91, 255, 0.82);
          outline-offset: 4px;
          border-radius: 12px;
          box-shadow: 0 0 0 8px rgba(31, 91, 255, 0.14);
          scroll-margin-top: 96px;
        }

        @media (max-width: 960px){
          .setup-grid{
            grid-template-columns: 1fr;
          }

          .tutorial-top-card{
            right: 10px;
            bottom: 82px;
            width: calc(100vw - 20px);
            max-width: none;
            max-height: 38vh;
          }

          .tutorial-top-card + div[data-testid="stHorizontalBlock"]{
            right: 10px;
            bottom: 10px;
            width: calc(100vw - 20px);
            max-width: none;
          }
        }

        @media (prefers-reduced-motion: reduce){
          *{
            animation-duration: 0s !important;
            animation-iteration-count: 1 !important;
            transition-duration: 0s !important;
            scroll-behavior: auto !important;
          }
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# ------------- Setup state -------------
stored_settings = load_settings()
defaults = default_paths()

if "paths_initialized" not in st.session_state:
    for key, _, _ in PATH_FIELDS:
        st.session_state[key] = stored_settings.get(key, defaults.get(key, ""))
    st.session_state.paths_initialized = True
    st.session_state.setup_done = bool(stored_settings)

def current_paths():
    return {key: st.session_state.get(key, "").strip() for key, _, _ in PATH_FIELDS}

def apply_paths(paths):
    ensure_paths(paths)
    save_settings(paths)

TUTORIAL_STEPS = [
    {
        "id": "input",
        "title": "Input source",
        "description": "Choose Example, Upload file, or Paste table in the sidebar.",
    },
    {
        "id": "input",
        "title": "Outlier ΔCq threshold",
        "description": "This sidebar threshold flags replicate disagreement in the clean-wells table.",
    },
    {
        "id": "input",
        "title": "Fit scope",
        "description": "Pick Gene (all plates) or Gene × Plate depending on whether plate-specific curves are needed.",
    },
    {
        "id": "input",
        "title": "Reference gene and quant mode",
        "description": "Reference gene is used for normalization. Quant mode switches between Absolute and ΔΔCt workflows.",
    },
    {
        "id": "mapping",
        "title": "Sample label column",
        "description": "Map the column that identifies biological sample labels.",
    },
    {
        "id": "mapping",
        "title": "Cq column",
        "description": "Map the measurement column containing Cq/Ct values for fitting and quantification.",
    },
    {
        "id": "mapping",
        "title": "Replicate column",
        "description": "Optional replicate-column mapping; leave none if replicates are already collapsed.",
    },
    {
        "id": "clean",
        "title": "Review & clean wells",
        "description": "Use keep checkboxes to exclude bad wells before averaging and fitting.",
    },
    {
        "id": "clean",
        "title": "Replicate averages",
        "description": "Validate replicate means/SD before proceeding to standards fitting.",
    },
    {
        "id": "standards",
        "title": "Label → concentration map",
        "description": "Fill all standard concentrations. Missing concentrations block proper fitting.",
    },
    {
        "id": "standards",
        "title": "Serial dilution auto-fill",
        "description": "Use top concentration + dilution factor + order to quickly prefill standard concentrations.",
    },
    {
        "id": "autoqc",
        "title": "Auto-QC thresholds",
        "description": "Tune R² target, min levels kept, allowed level drops, and efficiency bounds.",
    },
    {
        "id": "autoqc",
        "title": "Apply suggestions toggle",
        "description": "Apply Auto-QC exclusions to curve fit and downstream quantities when suggestions are acceptable.",
    },
    {
        "id": "autoqc",
        "title": "Curve review",
        "description": "Inspect slope, intercept, R², and efficiency; then view optional curve plots.",
    },
    {
        "id": "quantify",
        "title": "Quantify samples",
        "description": "Check computed quantities for sample wells using the selected fitted curve.",
    },
    {
        "id": "normalize",
        "title": "Normalize to reference gene",
        "description": "Review normalized quantities and confirm the selected reference gene is present.",
    },
    {
        "id": "export",
        "title": "Export final report",
        "description": "Download Excel once clean-wells, fit, quantification, and normalization outputs are validated.",
    },
]

if "tutorial_active" not in st.session_state:
    st.session_state["tutorial_active"] = False
if "tutorial_step_idx" not in st.session_state:
    st.session_state["tutorial_step_idx"] = 0

def _current_tutorial_step():
    idx = int(st.session_state.get("tutorial_step_idx", 0))
    idx = max(0, min(len(TUTORIAL_STEPS) - 1, idx))
    return TUTORIAL_STEPS[idx]

def _start_tutorial():
    st.session_state["tutorial_active"] = True
    st.session_state["tutorial_step_idx"] = 0

def _stop_tutorial():
    st.session_state["tutorial_active"] = False
    st.session_state["tutorial_step_idx"] = 0

def _next_tutorial_step():
    idx = int(st.session_state.get("tutorial_step_idx", 0))
    if idx >= len(TUTORIAL_STEPS) - 1:
        _stop_tutorial()
    else:
        st.session_state["tutorial_step_idx"] = idx + 1

def _prev_tutorial_step():
    idx = int(st.session_state.get("tutorial_step_idx", 0))
    st.session_state["tutorial_step_idx"] = max(0, idx - 1)

def _safe_rerun():
    rerun_fn = getattr(st, "rerun", None) or getattr(st, "experimental_rerun", None)
    if callable(rerun_fn):
        rerun_fn()

def render_tutorial_anchor(step_id, *, sidebar=False):
    safe_step = re.sub(r"[^a-z0-9_-]", "", str(step_id).lower())
    html = f"<div class='tutorial-anchor' data-step='{safe_step}' aria-hidden='true'></div>"
    if sidebar:
        st.sidebar.markdown(html, unsafe_allow_html=True)
    else:
        st.markdown(html, unsafe_allow_html=True)

def render_tutorial_focus():
    if not st.session_state.get("tutorial_active"):
        active_step = ""
    else:
        active_step = re.sub(r"[^a-z0-9_-]", "", str(_current_tutorial_step().get("id", "")).lower())

    components.html(
        f"""
<script>
(() => {{
  const activeStep = {json.dumps(active_step)};
  const host = window.parent;
  const doc = host?.document;
  if (!doc) return;

  const focusClass = 'easylab-tutorial-focus';
  const overlayId = 'easylab-tutorial-overlay';
  const minWidth = 160;
  const minHeight = 44;
  const minArea = 9000;
  let overlay = doc.getElementById(overlayId);
  if (!overlay) {{
    overlay = doc.createElement('div');
    overlay.id = overlayId;
    overlay.style.position = 'fixed';
    overlay.style.border = '3px solid rgba(31, 91, 255, 0.88)';
    overlay.style.borderRadius = '12px';
    overlay.style.boxShadow = '0 0 0 9999px rgba(10, 12, 22, 0.26), 0 0 0 6px rgba(31, 91, 255, 0.2)';
    overlay.style.pointerEvents = 'none';
    overlay.style.zIndex = '2147483600';
    overlay.style.opacity = '0';
    overlay.style.transition = 'opacity 120ms ease-out';
    doc.body.appendChild(overlay);
  }}

  const isVisible = (node) => {{
    if (!node || !(node instanceof host.HTMLElement)) return false;
    const style = host.getComputedStyle(node);
    if (style.display === 'none' || style.visibility === 'hidden' || style.opacity === '0') return false;
    const rect = node.getBoundingClientRect();
    return rect.width > 1 && rect.height > 1;
  }};

  const isMeaningfulTarget = (node) => {{
    if (!isVisible(node)) return false;
    const rect = node.getBoundingClientRect();
    if (rect.width < minWidth || rect.height < minHeight || rect.width * rect.height < minArea) return false;
    const hasAnchor = Boolean(node.querySelector('.tutorial-anchor'));
    const text = (node.textContent || '').trim();
    if (hasAnchor && text.length === 0 && !node.querySelector('input, select, textarea, button, [role="button"]')) {{
      return false;
    }}
    return true;
  }};

  const collectNeighbors = (node, maxCount = 12) => {{
    const out = [];
    if (!node || !node.parentElement) return out;
    let cur = node.nextElementSibling;
    while (cur && out.length < maxCount) {{
      out.push(cur);
      cur = cur.nextElementSibling;
    }}
    cur = node.previousElementSibling;
    while (cur && out.length < maxCount * 2) {{
      out.push(cur);
      cur = cur.previousElementSibling;
    }}
    return out;
  }};

  const dedupe = (nodes) => {{
    const seen = new Set();
    return nodes.filter((node) => {{
      if (!node || seen.has(node)) return false;
      seen.add(node);
      return true;
    }});
  }};

  const scoreTarget = (candidate, anchorRect) => {{
    const rect = candidate.getBoundingClientRect();
    const anchorCenter = anchorRect.top + anchorRect.height / 2;
    const candidateCenter = rect.top + rect.height / 2;
    const distance = Math.abs(candidateCenter - anchorCenter);
    const area = rect.width * rect.height;
    const largePenalty = area > 240000 ? 5000 : 0;
    return distance + largePenalty;
  }};

  const pickFocusTarget = (anchor) => {{
    const candidates = [];
    const anchorRect = anchor.getBoundingClientRect();
    const anchorContainer =
      anchor.closest('div[data-testid="stElementContainer"]') ||
      anchor.closest('div[data-testid="stVerticalBlock"]') ||
      anchor.parentElement;

    const callout = doc.querySelector(`.tutorial-callout[data-step="${{activeStep}}"]`);
    const calloutContainer =
      callout?.closest('div[data-testid="stElementContainer"]') ||
      callout?.closest('div[data-testid="stVerticalBlock"]') ||
      null;

    if (calloutContainer && isVisible(calloutContainer)) return calloutContainer;
    if (anchorContainer) {{
      candidates.push(...collectNeighbors(anchorContainer));
      candidates.push(anchorContainer);
      const parentBlock = anchorContainer.closest('div[data-testid="stVerticalBlock"]');
      if (parentBlock) candidates.push(parentBlock);
    }}

    const meaningful = dedupe(candidates).filter((candidate) => isMeaningfulTarget(candidate));
    meaningful.sort((a, b) => scoreTarget(a, anchorRect) - scoreTarget(b, anchorRect));
    if (meaningful.length > 0) return meaningful[0];

    if (calloutContainer && isVisible(calloutContainer)) return calloutContainer;
    return anchorContainer || callout || null;
  }};

  const updateOverlayFromFocus = () => {{
    const focused = doc.querySelector(`.${{focusClass}}`);
    if (!focused || !isVisible(focused)) {{
      overlay.style.opacity = '0';
      return;
    }}
    const rect = focused.getBoundingClientRect();
    if (rect.width < 8 || rect.height < 8) {{
      overlay.style.opacity = '0';
      return;
    }}
    const pad = 8;
    const top = Math.max(6, rect.top - pad);
    const left = Math.max(6, rect.left - pad);
    const width = Math.max(0, rect.width + pad * 2);
    const height = Math.max(0, rect.height + pad * 2);
    overlay.style.top = `${{top}}px`;
    overlay.style.left = `${{left}}px`;
    overlay.style.width = `${{width}}px`;
    overlay.style.height = `${{height}}px`;
    overlay.style.opacity = '1';
  }};

  if (!host.__easylabTutorialOverlayBound) {{
    host.addEventListener(
      'scroll',
      () => {{
        if (typeof host.__easylabTutorialOverlayUpdate === 'function') host.__easylabTutorialOverlayUpdate();
      }},
      true
    );
    host.addEventListener('resize', () => {{
      if (typeof host.__easylabTutorialOverlayUpdate === 'function') host.__easylabTutorialOverlayUpdate();
    }});
    host.__easylabTutorialOverlayBound = true;
  }}
  host.__easylabTutorialOverlayUpdate = updateOverlayFromFocus;

  doc.querySelectorAll(`.${{focusClass}}`).forEach((node) => node.classList.remove(focusClass));
  if (!activeStep) {{
    host.__easylabTutorialActiveStep = '';
    updateOverlayFromFocus();
    return;
  }}

  const anchor = doc.querySelector(`.tutorial-anchor[data-step="${{activeStep}}"]`);
  if (!anchor) {{
    updateOverlayFromFocus();
    return;
  }}
  const target = pickFocusTarget(anchor);
  if (!target) {{
    updateOverlayFromFocus();
    return;
  }}

  target.classList.add(focusClass);
  const targetRect = target.getBoundingClientRect();
  const targetOffscreen = targetRect.top < 96 || targetRect.bottom > host.innerHeight - 72;
  const stepChanged = host.__easylabTutorialActiveStep !== activeStep;
  if (stepChanged || targetOffscreen) {{
    target.scrollIntoView({{ behavior: 'smooth', block: 'center', inline: 'nearest' }});
    host.setTimeout(() => {{
      if (typeof host.__easylabTutorialOverlayUpdate === 'function') host.__easylabTutorialOverlayUpdate();
    }}, 260);
  }}
  host.__easylabTutorialActiveStep = activeStep;
  updateOverlayFromFocus();
}})();
</script>
        """,
        height=0,
        width=0,
    )

def render_tutorial_controls(*, sidebar=False):
    panel = st.sidebar if sidebar else st
    key_suffix = "sidebar" if sidebar else "top"
    card_class = "tutorial-panel" if sidebar else "tutorial-top-card"

    if not st.session_state.get("tutorial_active"):
        panel.caption("Guided run with focused highlights and Next/Back/Skip controls.")
        if panel.button("Start tutorial", key=f"tutorial_start_{key_suffix}", type="primary", use_container_width=True):
            _start_tutorial()
            _safe_rerun()
        return

    step = _current_tutorial_step()
    panel.markdown(
        (
            f"<div class='{card_class}'>"
            "<div class='kicker'>Guided tutorial</div>"
            f"<strong>Step {st.session_state.get('tutorial_step_idx', 0) + 1}/{len(TUTORIAL_STEPS)} · {step['title']}</strong>"
            f"<div>{step['description']}</div>"
            "</div>"
        ),
        unsafe_allow_html=True,
    )
    c_back, c_next, c_skip = panel.columns(3)
    back = c_back.button(
        "Back",
        key=f"tutorial_back_{key_suffix}",
        use_container_width=True,
        disabled=st.session_state["tutorial_step_idx"] == 0,
    )
    next_label = "Finish" if st.session_state["tutorial_step_idx"] >= len(TUTORIAL_STEPS) - 1 else "Next"
    nxt = c_next.button(next_label, key=f"tutorial_next_{key_suffix}", use_container_width=True, type="primary")
    skip = c_skip.button("Skip", key=f"tutorial_skip_{key_suffix}", use_container_width=True)
    if back:
        _prev_tutorial_step()
        _safe_rerun()
    if nxt:
        _next_tutorial_step()
        _safe_rerun()
    if skip:
        _stop_tutorial()
        _safe_rerun()

def render_tutorial_callout(step_id, *, sidebar=False):
    if not st.session_state.get("tutorial_active"):
        return
    step = _current_tutorial_step()
    if step["id"] != step_id:
        return
    safe_step = re.sub(r"[^a-z0-9_-]", "", str(step_id).lower())
    html = (
        f"<div class='tutorial-callout' data-step='{safe_step}'>"
        f"<p><strong>Step {st.session_state.get('tutorial_step_idx', 0) + 1}/{len(TUTORIAL_STEPS)} · {step['title']}</strong></p>"
        f"<p>{step['description']}</p>"
        "</div>"
    )
    if sidebar:
        st.sidebar.markdown(html, unsafe_allow_html=True)
    else:
        st.markdown(html, unsafe_allow_html=True)

if os.environ.get("E2E") == "1" and not st.session_state.get("setup_done"):
    try:
        apply_paths(defaults)
        st.session_state.setup_done = True
    except Exception:
        pass

# ------------- UI -------------
if not st.session_state.get("setup_done"):
    st.markdown("## First-run setup")
    st.markdown("Choose where qPCR Analysis stores data, attachments, exports, and sync content. These can be changed later in Settings.")

    with st.form("setup_form"):
        cols = st.columns(2)
        for idx, (key, label, helper) in enumerate(PATH_FIELDS):
            with cols[idx % 2]:
                st.text_input(label, key=key)
                st.caption(helper)

        error_message = st.session_state.get("setup_error", "")
        if error_message:
            st.markdown(f"<div class='setup-message error'>{error_message}</div>", unsafe_allow_html=True)

        use_defaults = st.form_submit_button("Use defaults")
        finish = st.form_submit_button("Finish setup")

    if use_defaults:
        for key, _, _ in PATH_FIELDS:
            st.session_state[key] = defaults.get(key, "")
        _safe_rerun()

    if finish:
        paths = current_paths()
        missing = [label for key, label, _ in PATH_FIELDS if not paths.get(key)]
        if missing:
            st.session_state.setup_error = "Please fill all paths before finishing setup."
            _safe_rerun()
        try:
            apply_paths(paths)
            st.session_state.setup_done = True
            st.session_state.setup_error = ""
            _safe_rerun()
        except Exception as exc:
            st.session_state.setup_error = f"Unable to create folders: {exc}"
            _safe_rerun()

    st.stop()

# ------------- UI -------------
st.sidebar.title("qPCR Analysis")
st.sidebar.caption("Load example wells, upload a file, or paste a table.")
render_tutorial_controls(sidebar=True)

with st.sidebar.expander("Settings", expanded=False):
    st.markdown("**Storage paths**")
    for key, label, helper in PATH_FIELDS:
        st.text_input(label, key=key)
        st.caption(helper)
    if st.button("Save settings"):
        try:
            apply_paths(current_paths())
            st.success("Settings saved and folders created.")
        except Exception as exc:
            st.error(f"Unable to save settings: {exc}")
    st.caption("License: All Rights Reserved.")

render_tutorial_anchor("input", sidebar=True)
render_tutorial_callout("input", sidebar=True)

input_mode = st.sidebar.radio("Input source", ["Example (sample-data/qpcr_example.csv)", "Upload file", "Paste table"], index=0)

uploaded_file = None
pasted_text = ""
if input_mode == "Upload file":
    uploaded_file = st.sidebar.file_uploader("CSV / TSV / Excel", type=["csv", "tsv", "txt", "xlsx"])
elif input_mode == "Paste table":
    st.sidebar.caption("Paste Excel/Sheets selection or CSV/TSV (first row must be headers).")
    pasted_text = st.sidebar.text_area(
        "Paste table",
        height=220,
        placeholder="Plate\tWell\tGene\tType\tLabel\tReplicate\tGroup\tCq\nPlate 1\tA1\thcar2\tSample\tA3\t1\tUNDETERMINED\t30.23937102\n..."
    )

st.sidebar.markdown("---")
outlier_thresh = st.sidebar.slider("Outlier ΔCq threshold", 0.1, 3.0, 0.75, 0.05)
plate_scope = st.sidebar.selectbox("Fit standard curve by", ["Gene (all plates)", "Gene × Plate"])
st.sidebar.markdown("---")
ref_gene = st.sidebar.text_input("Reference gene", value="gapdh")
quant_mode = st.sidebar.radio(
    "Quantification mode",
    ["Absolute (std curve)", "ΔΔCt (relative)"],
    horizontal=True,
    help="Absolute uses standard curves; ΔΔCt uses 2^-ΔΔCt relative expression."
)
if quant_mode != "Absolute (std curve)":
    st.info("Standard-curve fitting and Auto-QC are shown only in `Absolute (std curve)` mode.")
header_left, header_right = st.columns([0.72, 0.28], gap="medium")
with header_left:
    st.title("qPCR Analysis")
with header_right:
    render_tutorial_controls()

st.markdown(
    """
    <div class="hero">
        <div class="pill">Standard-curve workflow</div>
        <h4>Clean → replicate QC → fit → quantify → export</h4>
        <p>Load example wells, paste instrument output, flag outliers, prefill serial dilutions, and export a tidy Excel report.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.expander("Need to reformat your wells? AI helper prompt"):
    st.markdown(
        """
[ChatGPT](https://chat.openai.com/) · [Gemini](https://gemini.google.com/app) · [Grok](https://grok.com/)

**Prompt**
```text
Convert my table to CSV with headers: Plate, Well, Gene, Type, Label, Replicate, Cq. Normalize Well to A1-style, Replicate numeric, Cq numeric or NA, keep all rows, no invented data. Output CSV only.
```
Paste the returned CSV into the uploader or sidebar paste box.
        """
    )

# Load
render_tutorial_focus()
try:
    if input_mode.startswith("Example"):
        if not EXAMPLE_WELLS_PATH.exists():
            st.error(f"Missing example file: {EXAMPLE_WELLS_PATH}")
            st.stop()
        df_raw = pd.read_csv(EXAMPLE_WELLS_PATH)
    elif input_mode == "Upload file":
        if uploaded_file is None:
            st.info("Upload a CSV/Excel with columns: Plate, Well, Gene, Type, Label, Replicate, (Group optional), Cq.")
            st.stop()
        df_raw = read_from_file(uploaded_file)
    else:
        if not pasted_text.strip():
            st.info("Paste your table in the sidebar. Your last column (even if unnamed) will be auto-detected as Cq.")
            st.stop()
        df_raw = read_from_paste(pasted_text)
except Exception as e:
    st.error(f"Could not read the data: {e}")
    st.stop()

# Column mapping UI (allows manual selection of label/Cq/replicate headers)
st.sidebar.markdown("---")
render_tutorial_callout("mapping", sidebar=True)
render_tutorial_anchor("mapping", sidebar=True)
st.sidebar.subheader("Column mapping")

cols_available = list(df_raw.columns)
label_default = _guess_label_col(df_raw)
cq_default = "Cq" if "Cq" in df_raw.columns else _best_numeric_column_as_cq(df_raw) or cols_available[-1]
rep_guess = _guess_rep_col(df_raw)

label_col = st.sidebar.selectbox("Sample label column", options=cols_available, index=cols_available.index(label_default) if label_default in cols_available else 0)
cq_col = st.sidebar.selectbox("Cq column", options=cols_available, index=cols_available.index(cq_default) if cq_default in cols_available else len(cols_available)-1)
rep_options = ["<none>"] + cols_available
rep_col_choice = st.sidebar.selectbox("Replicate column (optional)", options=rep_options, index=rep_options.index(rep_guess) if rep_guess in cols_available else 0)

df_mapped = df_raw.copy()
df_mapped["Label"] = df_raw[label_col]
df_mapped["Cq"] = df_raw[cq_col]
if rep_col_choice != "<none>":
    df_mapped["Replicate"] = df_raw[rep_col_choice]
else:
    df_mapped["Replicate"] = 1

raw_df = coerce_columns(df_mapped)
render_tutorial_focus()

# Keep column + outliers
raw_df["keep"] = True
raw_df = mark_outliers(raw_df, threshold=outlier_thresh)

k1, k2, k3, k4 = st.columns(4)
k1.metric("Wells loaded", raw_df.shape[0])
k2.metric("Genes", raw_df["Gene"].nunique())
k3.metric("Standards", int((raw_df["Type"].str.lower()=="standard").sum()))
k4.metric("Outliers flagged", int(raw_df["Outlier"].sum()))

st.subheader("1) Review & clean wells")
render_tutorial_anchor("clean")
render_tutorial_callout("clean")
st.caption("Toggle **keep** to drop bad replicates. Outliers by ΔCq are flagged.")
edited = st.data_editor(
    raw_df,
    key="editor",
    num_rows="dynamic",
    column_config={
        "keep": st.column_config.CheckboxColumn("keep", default=True),
        "Outlier": st.column_config.CheckboxColumn("Outlier", disabled=True),
        "DeltaCq": st.column_config.NumberColumn("ΔCq", format="%.3f", disabled=True),
        "Cq": st.column_config.NumberColumn("Cq", format="%.3f")
    },
    hide_index=True,
    width="stretch",
)
clean_df = edited.copy()
clean_df.loc[~clean_df["keep"], "Cq"] = np.nan

# Collapse exact duplicate wells (same Plate/Well/Gene/Type/Label); average Cq across kept replicates.
extra_cols = [c for c in clean_df.columns if c not in VALID_COLS + ["keep","Outlier","DeltaCq"]]
well_keys = ["Plate","Well","Gene","Type","Label"]
agg_map = {"Cq": "mean", "keep": "any", "DeltaCq": "mean", "Outlier": "any"}
agg_map.update({c: "first" for c in extra_cols})
collapsed_df = clean_df.groupby(well_keys, as_index=False, dropna=False).agg(agg_map)

# ------------- replicate averages -------------
st.subheader("2) Replicate averages (auto-calculated)")
rep_stats = replicate_stats(collapsed_df)
st.dataframe(rep_stats, width="stretch", height=440)

# ------------- standards mapping -------------
if quant_mode == "Absolute (std curve)":
    st.subheader("3) Standards map (Label → Concentration)")
    render_tutorial_anchor("standards")
    render_tutorial_callout("standards")
    std_labels = clean_df.loc[(clean_df["Type"].str.lower()=="standard"), "Label"].dropna().unique().tolist()
    std_labels = sorted(std_labels, key=lambda x: (re.sub(r"\D","",x)=="", re.sub(r"\D","",x), x))
    auto_fill = _serial_dilution(std_labels, DEFAULT_TOP_CONC, DEFAULT_DILUTION_FACTOR, highest_first=True)
    default_map = pd.DataFrame({"Label": std_labels, "Concentration": auto_fill})
    st.caption("Default mapping uses a 4-fold serial dilution from the top standard; edit if your plate differs.")

    with st.expander("Auto-fill serial dilution"):
        top = st.number_input("Top concentration", min_value=0.0, value=DEFAULT_TOP_CONC, step=1.0, format="%.6f")
        factor = st.number_input("Dilution factor", min_value=1.0, value=DEFAULT_DILUTION_FACTOR, step=0.5, format="%.3f")
        order = st.selectbox("Order of labels", ["Std1 highest → StdN lowest", "Std1 lowest → StdN highest"])
        if st.button("Fill mapping from settings"):
            default_map["Concentration"] = _serial_dilution(
                std_labels,
                top,
                factor,
                highest_first=order.startswith("Std1 highest")
            )

    map_df = st.data_editor(
        default_map,
        key="stdmap",
        column_config={"Concentration": st.column_config.NumberColumn("Concentration", format="%.6f")},
        hide_index=True,
        width="stretch",
    )
    if map_df["Concentration"].isna().any():
        st.warning("Some standard concentrations are missing. Fill them before fitting curves.")

if quant_mode == "Absolute (std curve)":
    # ------------- fit curves -------------
    st.subheader("4) Fit standard curves")
    render_tutorial_anchor("autoqc")
    render_tutorial_callout("autoqc")
    if plate_scope == "Gene × Plate":
        std_input = collapsed_df[(collapsed_df["Type"].str.lower()=="standard") & (collapsed_df["keep"])].copy()
        std_input["Gene"] = std_input["Gene"].astype(str) + " | " + std_input["Plate"].astype(str)
    else:
        std_input = collapsed_df[(collapsed_df["Type"].str.lower()=="standard") & (collapsed_df["keep"])].copy()

    # Auto-QC: suggest replicate/level exclusions to improve curve fit (no data is deleted unless toggled).
    std_labels_for_qc = (
        std_input.merge(map_df.dropna(subset=["Label", "Concentration"]), on="Label", how="inner")["Label"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    default_min_levels_keep = max(2, min(6, len(std_labels_for_qc))) if std_labels_for_qc else 2

    with st.expander("Auto-QC details (review suggestions + tweak thresholds)", expanded=True):
        st.caption("Defaults are conservative: prefer fewer exclusions, keep at least 6 standard levels when available.")
        c1, c2, c3 = st.columns(3)
        qc_min_levels = c1.number_input("Min levels to keep", min_value=2, max_value=20, value=int(default_min_levels_keep), step=1)
        qc_max_drop_levels = c2.number_input("Max levels to drop", min_value=0, max_value=5, value=1, step=1)
        qc_r2_target = c3.number_input("R² target", min_value=0.0, max_value=0.999, value=0.98, step=0.001, format="%.3f")
        m1, m2 = st.columns(2)
        qc_objective_ui = m1.selectbox(
            "Auto-QC objective",
            [
                "Balanced (targets + minimal drops)",
                "Maximize R² (aggressive)",
            ],
            index=0,
            help="Balanced prioritizes meeting R²/efficiency targets with minimal exclusions. Maximize R² searches for the strongest linearity, then breaks ties by fewer drops.",
        )
        qc_rep_mode_ui = m2.selectbox(
            "Replicate drop candidates",
            [
                "Only outlier-like replicates",
                "Any single replicate per level",
            ],
            index=0,
            help="Outlier-like uses the ΔCq threshold below. Any single replicate per level is more aggressive and can improve R² when differences are subtle.",
        )
        qc_std_rep_thresh = st.slider(
            "Std replicate disagreement ΔCq threshold (Auto-QC)",
            min_value=0.1,
            max_value=3.0,
            value=float(outlier_thresh),
            step=0.05,
            help="Used only for Auto-QC replicate dropping within standards. This is independent of the sidebar Outlier ΔCq threshold (which only flags wells).",
        )
        e1, e2 = st.columns(2)
        qc_eff_min = e1.number_input("Efficiency min (%)", min_value=0.0, max_value=200.0, value=90.0, step=1.0, format="%.1f")
        qc_eff_max = e2.number_input("Efficiency max (%)", min_value=0.0, max_value=200.0, value=110.0, step=1.0, format="%.1f")
        qc_objective = "best_r2" if qc_objective_ui.startswith("Maximize") else "balanced"
        qc_rep_mode = "any_single_replicate" if qc_rep_mode_ui.startswith("Any single") else "outlier_only"

        qc_exclusions, qc_summary = suggest_standard_curve_exclusions(
            std_input,
            map_df,
            outlier_threshold=float(qc_std_rep_thresh),
            min_levels_to_keep=int(qc_min_levels),
            max_levels_to_drop=int(qc_max_drop_levels),
            r2_target=float(qc_r2_target),
            eff_min=float(qc_eff_min),
            eff_max=float(qc_eff_max),
            optimize_for=qc_objective,
            replicate_drop_mode=qc_rep_mode,
        )
        if qc_exclusions.empty:
            st.info("Auto-QC: no exclusions suggested for the current standards.")
        else:
            st.caption("Suggested curve metrics (baseline vs Auto-QC):")
            st.dataframe(qc_summary, width="stretch", height=220)
            st.caption("Suggested excluded wells (review before enabling Auto-QC):")
            st.dataframe(qc_exclusions, width="stretch", height=260)

    # Keep variables defined even if the expander stays collapsed in the UI.
    if "qc_exclusions" not in locals():
        qc_exclusions = pd.DataFrame()

    if qc_exclusions.empty:
        st.caption("Auto-QC: no exclusions suggested for the current standards.")
    else:
        st.caption(f"Auto-QC suggests excluding {qc_exclusions.shape[0]} standard wells. Open Auto-QC details to review.")

    if "auto_qc_apply_toggle" not in st.session_state:
        st.session_state["auto_qc_apply_toggle"] = False

    apply_col, reset_col, toggle_col = st.columns([1.2, 1.2, 2.6])
    apply_clicked = apply_col.button(
        "Apply suggestions",
        type="primary",
        disabled=qc_exclusions.empty,
        help="Apply the current Auto-QC suggestions to curve fitting and downstream calculations.",
    )
    reset_clicked = reset_col.button(
        "Use all kept standards",
        help="Disable Auto-QC exclusions and fit using all currently kept standard wells.",
    )
    if apply_clicked and not qc_exclusions.empty:
        st.session_state["auto_qc_apply_toggle"] = True
    if reset_clicked:
        st.session_state["auto_qc_apply_toggle"] = False

    use_auto_qc = toggle_col.toggle(
        "Apply Auto-QC exclusions to fit + downstream calculations",
        key="auto_qc_apply_toggle",
        help="When enabled, suggested exclusions are removed before fitting standard curves and all downstream quantities.",
    )
    if use_auto_qc and qc_exclusions.empty:
        st.warning("Auto-QC apply is enabled, but there are no suggested exclusions for the current settings.")
    elif use_auto_qc:
        st.success(f"Auto-QC applied: excluding {qc_exclusions.shape[0]} suggested standard well(s).")
    else:
        st.caption("Auto-QC not applied: using all kept standard wells for curve fitting.")

    std_fit_input = std_input
    if use_auto_qc and (not qc_exclusions.empty):
        drop_keys = qc_exclusions[["Gene", "Label", "Plate", "Well"]].drop_duplicates().copy()
        drop_keys["_drop"] = True
        merged = std_fit_input.merge(drop_keys, on=["Gene", "Label", "Plate", "Well"], how="left")
        std_fit_input = merged[merged["_drop"].isna()].drop(columns=["_drop"])

    curves_df, std_points = fit_standard_curve(std_fit_input, map_df)
    st.dataframe(curves_df, width="stretch", height=260)

    if curves_df.empty:
        st.info("No standard curves to fit yet. Add standards (Type=Standard) and fill the standards map.")
    else:
        bad = curves_df[(curves_df["R2"] < 0.98) | (~curves_df["Efficiency_%"].between(90, 110))]
        if bad.shape[0] > 0:
            st.warning("Some curves have R² < 0.98 or efficiency outside 90–110%. Consider revising outliers or concentrations.")

        show_plots = st.toggle("Show curve plots", value=True)
        if show_plots:
            st.caption("Each point: mean Cq at that standard level. Line: linear fit. Efficiency computed from slope.")
            for _, row in curves_df.iterrows():
                gene = row["Gene"]
                sub = std_points[std_points["Gene"] == gene]
                if sub.empty or pd.isna(row["slope"]):
                    continue
                x = sub["log10_conc"].values
                y = sub["meanCq"].values
                fig = make_std_curve_figure(
                    gene=str(gene),
                    x=x,
                    y=y,
                    slope=float(row["slope"]),
                    intercept=float(row["intercept"]),
                    r2=float(row["R2"]),
                    eff=float(row["Efficiency_%"]),
                )
                st.pyplot(fig, width="stretch")

    # ------------- quantify samples -------------
    st.subheader("5) Quantify samples")
    render_tutorial_anchor("quantify")
    render_tutorial_callout("quantify")
    if plate_scope == "Gene × Plate":
        samp_input = collapsed_df[(collapsed_df["Type"].str.lower()=="sample") & (collapsed_df["keep"])].copy()
        samp_input["Gene"] = samp_input["Gene"].astype(str) + " | " + samp_input["Plate"].astype(str)
    else:
        samp_input = collapsed_df[(collapsed_df["Type"].str.lower()=="sample") & (collapsed_df["keep"])].copy()

    quant_df = quantify_samples(samp_input, curves_df)
    st.dataframe(quant_df.head(60), width="stretch", height=440)

    # ------------- normalize -------------
    st.subheader("6) Normalize to reference gene")
    render_tutorial_anchor("normalize")
    render_tutorial_callout("normalize")
    norm_df = normalize_to_ref(quant_df, ref_gene=ref_gene)
    if norm_df["RefQty"].isna().all():
        st.info(f"Reference gene '{ref_gene}' missing in standards/samples; add a ref gene to see normalized values.")
    st.dataframe(norm_df, width="stretch", height=440)

    # Per-well normalized view (attach all metadata for export)
    per_well_norm = quant_df.merge(norm_df[["Label","RefQty"]], on="Label", how="left")
    per_well_norm["Norm_Qty"] = per_well_norm["Quantity"] / per_well_norm["RefQty"]
else:
    # ΔΔCt workflow
    st.subheader("4) ΔΔCt (relative)")
    ref_rows = collapsed_df[collapsed_df["Gene"].str.lower() == ref_gene.lower()].copy()
    ref_mean = ref_rows.groupby("Label", dropna=False)["Cq"].mean().reset_index().rename(columns={"Cq":"RefCq"})
    dd_df = collapsed_df.merge(ref_mean, on="Label", how="left")
    dd_df = dd_df[dd_df["Gene"].str.lower() != ref_gene.lower()].copy()
    dd_df["DeltaCt"] = dd_df["Cq"] - dd_df["RefCq"]
    all_labels = sorted(dd_df["Label"].dropna().unique().tolist())
    calibrators = st.multiselect("Calibrator labels", options=all_labels, default=all_labels[:1])
    exclude = st.multiselect("Exclude from calibrator", options=calibrators, default=[])
    calib_set = [l for l in calibrators if l not in exclude]
    calib_means = (
        dd_df[dd_df["Label"].isin(calib_set)]
        .groupby("Gene")["DeltaCt"]
        .mean()
        .reset_index()
        .rename(columns={"DeltaCt":"CalibDeltaCt"})
    )
    dd_df = dd_df.merge(calib_means, on="Gene", how="left")
    dd_df["DeltaDeltaCt"] = dd_df["DeltaCt"] - dd_df["CalibDeltaCt"]
    dd_df["FoldChange"] = 2 ** (-dd_df["DeltaDeltaCt"])
    st.dataframe(dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]], width="stretch", height=440)
    dd_csv = dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]].to_csv(index=False).encode()
    st.download_button("Download ΔΔCt table (CSV)", dd_csv, file_name="ddct_results.csv", mime="text/csv")

    # placeholders to keep downstream variables defined
    map_df = pd.DataFrame()
    curves_df = pd.DataFrame()
    std_points = pd.DataFrame()
    quant_df = pd.DataFrame()
    norm_df = dd_df.rename(columns={"FoldChange":"Norm_Qty"})
    per_well_norm = dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]].copy()
    per_well_norm["slope"] = np.nan
    per_well_norm["intercept"] = np.nan
    per_well_norm["pred_log10Q"] = np.nan
    per_well_norm["Quantity"] = np.nan
    per_well_norm["RefQty"] = np.nan
    per_well_norm["Norm_Qty"] = per_well_norm["FoldChange"]

# Build a metadata-rich per-well export (one row per well)
meta_cols = ["Plate","Well","Gene","Type","Label","Cq","keep","DeltaCq","Outlier"] + extra_cols
meta_frame = collapsed_df[meta_cols].copy()

# Attach curve + quantity metadata so the Excel "PerWell_Normalized" sheet is fully informative.
merge_cols = ["Plate","Well","Gene","Label"]
per_well_export_cols = merge_cols + ["slope","intercept","pred_log10Q","Quantity","RefQty","Norm_Qty"]
export_norm_df = (
    meta_frame
    .merge(per_well_norm[per_well_export_cols], on=merge_cols, how="left")
    # Ensure only one row per well/gene/label (collapse any accidental duplicate copies)
    .drop_duplicates(subset=["Plate","Well","Gene","Label"])
    # Order by gene first so reference and targets are grouped, not interleaved.
    .sort_values(["Gene","Plate","Well","Label"])
)

# ------------- export -------------
st.subheader("7) Export")
render_tutorial_anchor("export")
render_tutorial_callout("export")
# Pass per-sample normalized table separately from the per-well export table
excel_bytes = download_excel(clean_df, rep_stats, map_df, curves_df, std_points, quant_df, norm_df, export_norm_df)
st.download_button(
    label="Download Excel report",
    data=excel_bytes,
    file_name=f"qpcr_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    type="primary",
)

render_tutorial_focus()
st.success("Done. Paste/upload → clean → averages → fit → quantify → normalize → export.")
