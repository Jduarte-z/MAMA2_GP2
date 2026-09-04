#!/usr/bin/env python3

import os
from datetime import datetime

import pandas as pd
import gwaslab as gl


HIGHLIGHT = True
COLORS = ["#000000", "#ABABAB"]

GWAS_TOP = "/home/jupyter/workspace/ws_files/MAMA_2/mama_gp2_non-eur_only/full_mr-mega/mama_gp2Only_non-eur.gwaslab.tsv.gz"
GWAS_BOTTOM = "/home/jupyter/workspace/ws_files/MAMA_2/mama_gp2_non-eur_only/full_random-effects/download/gwama_random_eaf0.01.tsv.gz"

NOVEL_TOP = "../full_mr-mega/out_file_novel"
NOVEL_BOTTOM = "../full_random-effects/out_file_novel"

TOP_LABEL = "MR-MEGA"
BOTTOM_LABEL = "Random effects"
PLOT_TITLE = "Only GP2 non-EUR POPs"

OUT_PLOT = "gwama_miami.png"

FIGSIZE = (9, 5)
BASE_W = 9
SCALE = FIGSIZE[0] / BASE_W

DPI = 600
PAD_IN = 0.25 * SCALE

CUT = 11
REPEL_FORCE = 0.03
YLIM_SCALE = 2


def load_pinpoints(path: str) -> list[str]:
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

    return list(dict.fromkeys(snps))


for path in [GWAS_TOP, GWAS_BOTTOM]:
    if not os.path.exists(path):
        raise SystemExit(f"Input GWAS file not found: {path}")

if HIGHLIGHT:
    for path in [NOVEL_TOP, NOVEL_BOTTOM]:
        if not os.path.exists(path):
            raise SystemExit(f"Novel loci file not found: {path}")

out_dir = os.path.dirname(os.path.abspath(OUT_PLOT))
os.makedirs(out_dir, exist_ok=True)

run_start = datetime.now()

print(
    f"[{run_start.strftime('%H:%M:%S')}] "
    f"Processing Miami plot..."
)

if HIGHLIGHT:
    pinpoints_top = load_pinpoints(NOVEL_TOP)
    pinpoints_bottom = load_pinpoints(NOVEL_BOTTOM)

    print(f"Top GWAS novel SNPs:    {len(pinpoints_top)}")
    print(f"Bottom GWAS novel SNPs: {len(pinpoints_bottom)}")
else:
    pinpoints_top = []
    pinpoints_bottom = []


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


fig, log = gl.plot_miami2(
    gl_top,
    gl_bottom,
    id1="SNPID",
    id2="SNPID",
    suffixes=["_TOP", "_BOTTOM"],
    build="38",
    mode="m",
    cut=CUT,
    sig_line=True,
    sig_level=5e-8,
    additional_line=[1e-6],
    additional_line_color=["gray"],
    titles=None,
    font_family="DejaVu Sans",
    fontsize=7 * SCALE,
    anno1="GENENAME",
    anno2="GENENAME",
    anno_style="expand",
    anno_fontsize=7 * SCALE,
    colors=COLORS,
    highlight1=pinpoints_top if pinpoints_top else None,
    highlight2=pinpoints_bottom if pinpoints_bottom else None,
    highlight_color1="#E70B0B",
    highlight_color2="#E70B0B",
    repel_force=REPEL_FORCE,
    xtight=False,
    same_ylim=True,
    fig_kwargs={
        "figsize": FIGSIZE,
        "dpi": DPI,
    },
    save=None
)


ax_top, ax_bottom = fig.axes

ymin_t, ymax_t = ax_top.get_ylim()
ax_top.set_ylim(
    ymin_t,
    ymax_t * YLIM_SCALE
)

ymin_b, ymax_b = ax_bottom.get_ylim()
ax_bottom.set_ylim(
    ymin_b,
    ymax_b * YLIM_SCALE
)

chr_tick = 7 * SCALE

ax_top.tick_params(
    axis="both",
    which="both",
    labelsize=chr_tick
)

ax_bottom.tick_params(
    axis="both",
    which="both",
    labelsize=chr_tick
)

for ax, label in zip(
    (ax_top, ax_bottom),
    (TOP_LABEL, BOTTOM_LABEL)
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

fig.subplots_adjust(
    left=0.08,
    right=0.92,
    top=0.88,
    bottom=0.10
)

fig.suptitle(
    PLOT_TITLE,
    fontsize=14,
    fontweight="bold",
    x=0.02,
    y=0.97,
    ha="left"
)

fig.savefig(
    OUT_PLOT,
    dpi=DPI,
    facecolor="white",
    bbox_inches="tight",
    pad_inches=PAD_IN
)

run_end = datetime.now()
elapsed = (run_end - run_start).total_seconds()

print(
    f"[{run_end.strftime('%H:%M:%S')}] "
    f"Finished Miami plot. Elapsed: {elapsed:.1f}s"
)

print(f"Output: {OUT_PLOT}")
