# plot Miami plot + highlight novel hits
# gwaslab v4.1.6

import os
from datetime import datetime

import pandas as pd
import gwaslab as gl
import matplotlib as mpl


# ============================================================
# USER CONFIGURATION
# ============================================================

HIGHLIGHT = True
COLORS = ["#000000", "#ABABAB"]


# -------------------------
# Input GWAS files
# -------------------------

GWAS_TOP = "/path/to/first_gwas.gwaslab.tsv.gz"
GWAS_BOTTOM = "/path/to/second_gwas.gwaslab.tsv.gz"


# -------------------------
# Novel loci files
# -------------------------

NOVEL_TOP = "/path/to/first_gwas_novel_loci.tsv"
NOVEL_BOTTOM = "/path/to/second_gwas_novel_loci.tsv"


# -------------------------
# Plot labels
# -------------------------

TOP_LABEL = "Random effects"
BOTTOM_LABEL = "Fixed effects"

PLOT_TITLE = "SAIGE Phase1 & Phase2 meta-analysis"


# -------------------------
# Output
# -------------------------

OUT_DIR = "./"
OUT_PLOT = os.path.join(OUT_DIR, "gwama_miami.png")

os.makedirs(OUT_DIR, exist_ok=True)


# -------------------------
# Figure configuration
# -------------------------

FIGSIZE = (9, 5)
BASE_W = 9
SCALE = FIGSIZE[0] / BASE_W

DPI = 600
PAD_IN = 0.25 * SCALE


# ============================================================
# FUNCTIONS
# ============================================================

def load_pinpoints(path: str) -> list[str]:
    """Load SNP IDs from a novel-loci TSV file."""

    if not os.path.exists(path):
        print(f"WARNING: Novel loci file not found: {path}")
        return []

    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as e:
        print(f"WARNING: Could not read novel loci file: {path}")
        print(f"         {e}")
        return []

    if df.empty:
        return []

    # Prefer SNPID, otherwise use first column
    col = "SNPID" if "SNPID" in df.columns else df.columns[0]

    snps = (
        df[col]
        .dropna()
        .astype(str)
        .str.strip()
    )

    snps = [
        snp for snp in snps.tolist()
        if snp and snp.lower() != "nan"
    ]

    # Deduplicate while preserving order
    return list(dict.fromkeys(snps))


# ============================================================
# START
# ============================================================

script_start = datetime.now()
print(
    f"Script started: "
    f"{script_start.strftime('%Y-%m-%d %H:%M:%S')}\n"
)

run_start = datetime.now()
print(
    f"[{run_start.strftime('%H:%M:%S')}] "
    f"Processing Miami plot..."
)


# ============================================================
# LOAD NOVEL LOCI
# ============================================================

if HIGHLIGHT:
    pinpoints_top = load_pinpoints(NOVEL_TOP)
    pinpoints_bottom = load_pinpoints(NOVEL_BOTTOM)

    print(f"Top GWAS novel SNPs:    {len(pinpoints_top)}")
    print(f"Bottom GWAS novel SNPs: {len(pinpoints_bottom)}")

else:
    pinpoints_top = []
    pinpoints_bottom = []


# ============================================================
# LOAD GWAS
# ============================================================

gl_top = gl.Sumstats(
    GWAS_TOP,
    build="38",
    fmt="gwaslab"
)

gl_bottom = gl.Sumstats(
    GWAS_BOTTOM,
    build="38",
    fmt="gwaslab"
)


# ============================================================
# MIAMI PLOT
# ============================================================

fig, log = gl.plot_miami2(

    # GWAS datasets
    gl_top,
    gl_bottom,

    id1="SNPID",
    id2="SNPID",

    suffixes=["_TOP", "_BOTTOM"],

    build="38",
    mode="m",

    cut=11,

    # Significance lines
    sig_line=True,
    sig_level=5e-8,

    additional_line=[1e-6],
    additional_line_color=["gray"],

    # Titles
    titles=None,

    # Fonts
    font_family="DejaVu Sans",
    fontsize=8 * SCALE,

    # Annotation
    anno1="GENENAME",
    anno2="GENENAME",

    anno_style="right",
    anno_fontsize=10 * SCALE,

    # Colors
    colors=COLORS,

    # Highlight novel loci
    highlight1=pinpoints_top if pinpoints_top else None,
    highlight2=pinpoints_bottom if pinpoints_bottom else None,

    highlight_color1="#E70B0B",
    highlight_color2="#E70B0B",

    repel_force=0.1,

    # Axes
    xtight=False,
    same_ylim=True,

    # Figure
    fig_kwargs={
        "figsize": FIGSIZE,
        "dpi": DPI,
    },

    save=None
)


# ============================================================
# POST-PROCESS FIGURE
# ============================================================

ax_top, ax_bottom = fig.axes


# -------------------------
# Expand annotation space
# -------------------------

ymin_t, ymax_t = ax_top.get_ylim()
ax_top.set_ylim(
    ymin_t,
    ymax_t * 1.5
)

ymin_b, ymax_b = ax_bottom.get_ylim()
ax_bottom.set_ylim(
    ymin_b,
    ymax_b * 1.25
)


# -------------------------
# Tick labels
# -------------------------

CHR_TICK = 7 * SCALE

ax_top.tick_params(
    axis="both",
    which="both",
    labelsize=CHR_TICK
)

ax_bottom.tick_params(
    axis="both",
    which="both",
    labelsize=CHR_TICK
)


# -------------------------
# Side labels
# -------------------------

labels = [
    TOP_LABEL,
    BOTTOM_LABEL
]

for ax, label in zip(
    (ax_top, ax_bottom),
    labels
):
    ax.set_ylabel("")

    bbox = ax.get_position()
    y_center = (bbox.y0 + bbox.y1) / 2

    fig.text(
        0.97,
        y_center,
        label,
        rotation=90,
        va="center",
        ha="center",
        fontsize=12,
    )


# -------------------------
# Figure margins
# -------------------------

fig.subplots_adjust(
    left=0.08,
    right=0.92,
    top=0.88,
    bottom=0.10
)


# -------------------------
# Main title
# -------------------------

fig.suptitle(
    PLOT_TITLE,
    fontsize=14,
    fontweight="bold",
    x=0.02,
    y=0.97,
    ha="left"
)


# ============================================================
# SAVE
# ============================================================

fig.savefig(
    OUT_PLOT,
    dpi=DPI,
    facecolor="white",
    bbox_inches="tight",
    pad_inches=PAD_IN
)


# ============================================================
# FINISH
# ============================================================

run_end = datetime.now()
elapsed = (run_end - run_start).seconds

print(
    f"[{run_end.strftime('%H:%M:%S')}] "
    f"Finished Miami plot. Elapsed: {elapsed}s"
)

print(f"Output: {OUT_PLOT}")

print(
    f"\nScript finished: "
    f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
)
