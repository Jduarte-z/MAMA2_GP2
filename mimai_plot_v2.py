# Miami plot, MR-MEGA (top) vs random effects (bottom), novel loci in red
# gwaslab 4.1.6

import os
from datetime import datetime

import numpy as np
import pandas as pd
import gwaslab as gl
import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgb
from matplotlib.path import Path
import matplotlib.transforms as mtransforms


HIGHLIGHT = True
COLORS = ["#000000", "#ABABAB"]
HIGHLIGHT_COLOR = "#E70B0B"

# inputs
GWAS_TOP = "../r12_mama_mr-mega_10pops_g2019sAJcond/liftover/r12_mama_mr-mega_10pops_g2019sAJcond_hg19_ndf_fixP_headerFixed.tsv.gz"
GWAS_BOTTOM = "../r12_mama_randomEff_20pops_g2019sAJcond/liftover/r12_mama_randomEff_10pops_g2019sAJcond_hg19_headerfixed.tsv.gz"

NOVEL_TOP = "../r12_mama_mr-mega_10pops_g2019sAJcond/liftover/out_file_novel_v2"
NOVEL_BOTTOM = "../r12_mama_randomEff_20pops_g2019sAJcond/liftover/out_file_novel_v2"

# labels
TOP_LABEL = "MR-MEGA"
BOTTOM_LABEL = "Random effects"
PLOT_TITLE = "GP2's 10 populations (AJ cond. G2019S)"  # None = no title
Y_LABEL = r"$-\log_{10}(P)$"

# output
OUT_DIR = "./"
OUT_STEM = "gwama_miami"
SAVE_PDF = True

os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, f"{OUT_STEM}.png")
OUT_PDF = os.path.join(OUT_DIR, f"{OUT_STEM}.pdf")

# y axis
SIG_LEVEL = 5e-9
CUT = 15
CUT_FACTOR = 8
LOWER_TICK_STEP = 10
UPPER_TICKS = [80, 120, 160, 200, 300]
HEADROOM = 1.12
SKIP = 0

# figure
FIGSIZE = (7.2, 4.4)
DPI = 600

MARKER_SCALE = 0.30
HIGHLIGHT_MARKER_SCALE = 0.60

FS_TICK = 6
FS_AXIS = 7
FS_ANNO = 6
FS_PANEL = 8
FS_TITLE = 9

PREFERRED_FONTS = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]


def load_pinpoints(path):
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
    snps = df[col].dropna().astype(str).str.strip()
    snps = [s for s in snps.tolist() if s and s.lower() != "nan"]
    return list(dict.fromkeys(snps))


def pick_font(preferred):
    available = {f.name for f in font_manager.fontManager.ttflist}
    return next((f for f in preferred if f in available), "DejaVu Sans")


def raw_max_mlog10p(sumstats):
    df = sumstats.data
    if "MLOG10P" in df.columns:
        v = pd.to_numeric(df["MLOG10P"], errors="coerce")
        v = v[np.isfinite(v)]
        if len(v):
            return float(v.max())
    p = pd.to_numeric(df["P"], errors="coerce")
    p = p[np.isfinite(p) & (p > 0)]
    return float(-np.log10(p.min()))


def point_collections(ax):
    return [c for c in ax.collections
            if isinstance(c, PathCollection) and len(c.get_offsets()) > 0]


def plotted_y(ax):
    ys = [np.asarray(c.get_offsets())[:, 1] for c in point_collections(ax)]
    return np.concatenate(ys) if ys else np.array([])


def restyle_points(ax):
    hl_rgb = np.array(to_rgb(HIGHLIGHT_COLOR))
    for c in point_collections(ax):
        fc = c.get_facecolors()
        is_hl = len(fc) > 0 and np.allclose(fc[0][:3], hl_rgb, atol=0.02)
        sizes = c.get_sizes()
        if len(sizes):
            c.set_sizes(sizes * (HIGHLIGHT_MARKER_SCALE if is_hl else MARKER_SCALE))
        c.set_linewidths(0)
        c.set_rasterized(True)
        if is_hl:
            c.set_zorder(max(c.get_zorder(), 3))


def restyle_hlines(ax, cut):
    # drop the cut line, thin the sig line
    for line in list(ax.lines):
        y = np.asarray(line.get_ydata(), dtype=float)
        if y.size == 0 or not np.allclose(y, y[0]):
            continue
        if np.isclose(abs(y[0]), cut):
            line.remove()
        else:
            line.set_linewidth(0.7)


# same compression gwaslab uses: cut + (y - cut) / cutfactor
def to_plot_units(v, cut, cutfactor):
    v = np.asarray(v, dtype=float)
    return np.where(v > cut, cut + (v - cut) / cutfactor, v)


def to_raw_units(y, cut, cutfactor):
    return cut + (y - cut) * cutfactor if y > cut else y


def set_panel_ylim(ax, extent, sign):
    lo, hi = sorted([0.0, sign * extent])
    if ax.yaxis_inverted():
        ax.set_ylim(hi, lo)
    else:
        ax.set_ylim(lo, hi)


def draw_break_mark(ax, y):
    verts = [(-1, -1.3), (1, -0.3), (-1, 0.3), (1, 1.3)]
    codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
    trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.plot([0], [y], marker=Path(verts, codes), markersize=7,
            color="black", markeredgewidth=0.8, linestyle="none",
            transform=trans, clip_on=False, zorder=10)


script_start = datetime.now()
print(f"Script started: {script_start:%Y-%m-%d %H:%M:%S}\n")

FONT = pick_font(PREFERRED_FONTS)
mpl.rcParams["font.family"] = FONT
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
print(f"Font: {FONT}")

if HIGHLIGHT:
    pinpoints_top = load_pinpoints(NOVEL_TOP)
    pinpoints_bottom = load_pinpoints(NOVEL_BOTTOM)
    print(f"Top GWAS novel SNPs:    {len(pinpoints_top)}")
    print(f"Bottom GWAS novel SNPs: {len(pinpoints_bottom)}")
else:
    pinpoints_top, pinpoints_bottom = [], []

gl_top = gl.Sumstats(GWAS_TOP, build="19", fmt="gwaslab")
gl_bottom = gl.Sumstats(GWAS_BOTTOM, build="19", fmt="gwaslab")

raw_max = {
    "top": raw_max_mlog10p(gl_top),
    "bottom": raw_max_mlog10p(gl_bottom),
}
print(f"Max -log10(P): top {raw_max['top']:.1f}, bottom {raw_max['bottom']:.1f}")


fig, log = gl.plot_miami2(
    gl_top,
    gl_bottom,
    id1="SNPID",
    id2="SNPID",
    suffixes=["_TOP", "_BOTTOM"],
    build="19",
    mode="m",
    cut=CUT,
    cutfactor=CUT_FACTOR,
    skip=SKIP,
    sig_line=True,
    sig_level=SIG_LEVEL,
    font_family=FONT,
    fontsize=FS_TICK,
    anno1="GENENAME",
    anno2="GENENAME",
    anno_set1=pinpoints_top or None,
    anno_set2=pinpoints_bottom or None,
    anno_style="right",
    anno_fontsize=FS_ANNO,
    colors=COLORS,
    highlight1=pinpoints_top or None,
    highlight2=pinpoints_bottom or None,
    highlight_color1=HIGHLIGHT_COLOR,
    highlight_color2=HIGHLIGHT_COLOR,
    repel_force=0.03,
    xtight=False,
    same_ylim=True,
    fig_kwargs={"figsize": FIGSIZE, "dpi": DPI},
    save=None,
)

axes = fig.axes
if len(axes) < 2:
    raise RuntimeError(f"Expected 2 axes from plot_miami2, got {len(axes)}")
ax_top, ax_bottom = axes[0], axes[1]
panels = {"top": ax_top, "bottom": ax_bottom}

for ax in panels.values():
    restyle_points(ax)
    restyle_hlines(ax, CUT)

# check the compression gwaslab applied matches CUT_FACTOR
sign, plot_max = {}, {}
for name, ax in panels.items():
    y = plotted_y(ax)
    y = y[np.isfinite(y)]
    sign[name] = 1 if (y.size == 0 or np.median(y) >= 0) else -1
    plot_max[name] = float(np.max(np.abs(y))) if y.size else 0.0

    if raw_max[name] > CUT + 1 and plot_max[name] > CUT:
        inferred = (raw_max[name] - CUT) / (plot_max[name] - CUT)
        if abs(inferred / CUT_FACTOR - 1) > 0.05:
            print(f"WARNING: {name} panel looks compressed by ~{inferred:.2f}, "
                  f"not {CUT_FACTOR}. Tick labels above the break may be off; "
                  f"check the cutfactor argument in your gwaslab version.")

# y range, ticks and break marks
extent = max(plot_max.values()) * HEADROOM
raw_extent = to_raw_units(extent, CUT, CUT_FACTOR)

tick_raw = list(np.arange(0, CUT + 1e-9, LOWER_TICK_STEP))
tick_raw += [t for t in UPPER_TICKS if CUT < t <= raw_extent]
tick_pos = to_plot_units(tick_raw, CUT, CUT_FACTOR)
tick_lab = [f"{int(round(t))}" for t in tick_raw]

compressed = extent > CUT

for name, ax in panels.items():
    set_panel_ylim(ax, extent, sign[name])
    ax.set_yticks(sign[name] * tick_pos)
    ax.set_yticklabels(tick_lab)
    ax.minorticks_off()
    ax.tick_params(axis="both", which="both", labelsize=FS_TICK,
                   width=0.6, length=2.5)
    ax.set_ylabel(Y_LABEL, fontsize=FS_AXIS)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    if compressed:
        draw_break_mark(ax, sign[name] * CUT)

for ax in panels.values():
    ax.xaxis.label.set_size(FS_AXIS)

fig.subplots_adjust(left=0.09, right=0.91, top=0.90, bottom=0.10)

for ax, label in zip((ax_top, ax_bottom), (TOP_LABEL, BOTTOM_LABEL)):
    bbox = ax.get_position()
    fig.text(0.955, (bbox.y0 + bbox.y1) / 2, label, rotation=90,
             va="center", ha="center", fontsize=FS_PANEL)

# title goes above any gene label sticking out of the top panel
if PLOT_TITLE:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_h = fig.bbox.height
    text_tops = [
        t.get_window_extent(renderer).y1 / fig_h
        for t in ax_top.texts
        if t.get_visible() and t.get_text()
    ]
    axes_top = ax_top.get_position().y1
    title_y = max([axes_top + 0.04] + [y + 0.03 for y in text_tops])

    fig.suptitle(PLOT_TITLE, fontsize=FS_TITLE, fontweight="bold",
                 x=0.02, y=title_y, ha="left", va="bottom")

save_kw = dict(dpi=DPI, facecolor="white", bbox_inches="tight", pad_inches=0.05)

fig.savefig(OUT_PNG, **save_kw)
print(f"Output: {OUT_PNG}")

if SAVE_PDF:
    fig.savefig(OUT_PDF, **save_kw)
    print(f"Output: {OUT_PDF}")

elapsed = (datetime.now() - script_start).total_seconds()
print(f"\nScript finished: {datetime.now():%Y-%m-%d %H:%M:%S} (elapsed {elapsed:.0f}s)")
