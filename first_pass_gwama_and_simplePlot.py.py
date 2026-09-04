#!/usr/bin/env python3

import os
import sys
import argparse
import gwaslab as gl


REF_DIR = "/home/jupyter/workspace/ws_files/MAMA_2/gwaslab_referenceFiles"

THREADS = 4
GENOME_BUILD = "38"
INPUT_SEP = "\t"
NA_VALUES = ["NA", ".", "nan", ""]
VERBOSE = True

COLUMN_MAP = {
    "snpid": "rs_number",
    "ea": "reference_allele",
    "eaf": "eaf",
    "nea": "other_allele",
    "n": "n_samples",
    "OR": "OR",
    "OR_95L": "OR_95L",
    "OR_95U": "OR_95U",
    "se": "OR_se",
    "p": "p-value",
    "direction": "effects",
    "i2": "i2",
    "z": "z",
    "mlog10p": "_-log10_p-value"
}

EXTRA_COLS = [
    "OR_se",
    "q_statistic",
    "q_p-value"
]

OUTPUT_EXTRA_COLS = EXTRA_COLS

WINDOW_SIZE_KB = 500
SIG_LEVEL = 5e-8
ANNOTATE_LEADS = True
ANNO_SOURCE = "ensembl"

DROP_BETA_SE = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process GWAMA summary statistics with gwaslab."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="GWAMA result file."
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix/path for the gwaslab output."
    )
    parser.add_argument(
        "--ref-dir",
        default=REF_DIR,
        help="gwaslab reference directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=THREADS,
        help="Number of threads."
    )
    parser.add_argument(
        "--plot",
        default=None,
        help="Output path for the Manhattan/QQ plot."
    )
    parser.add_argument(
        "--title",
        default="GWAMA meta-analysis",
        help="Plot title."
    )
    parser.add_argument(
        "--lead-variants-out",
        default=None,
        help="Output path for lead variants."
    )
    parser.add_argument(
        "--window-size-kb",
        type=int,
        default=WINDOW_SIZE_KB,
        help="Window size in kb for lead variant extraction."
    )
    parser.add_argument(
        "--sig-level",
        type=float,
        default=SIG_LEVEL,
        help="Significance threshold for lead variant extraction."
    )
    parser.add_argument(
        "--drop-beta-se",
        action="store_true",
        help="Drop BETA and SE from the final output."
    )
    parser.add_argument(
        "--no-annotate-leads",
        action="store_true",
        help="Do not annotate lead variants."
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable verbose gwaslab output."
    )

    return parser.parse_args()


args = parse_args()

INPUT_SUMSTATS = args.input
OUTPUT_PREFIX = args.output_prefix
REF_DIR = args.ref_dir
THREADS = args.threads
WINDOW_SIZE_KB = args.window_size_kb
SIG_LEVEL = args.sig_level
DROP_BETA_SE = args.drop_beta_se
ANNOTATE_LEADS = not args.no_annotate_leads
VERBOSE = not args.quiet

LEAD_VARIANTS_OUT = (
    args.lead_variants_out
    or f"{OUTPUT_PREFIX}.lead_variants.tsv"
)

plot = args.plot or f"{OUTPUT_PREFIX}.png"
title = args.title


def check_paths():
    required = {
        "Input summary statistics": INPUT_SUMSTATS,
        "gwaslab data directory": REF_DIR,
    }

    missing = [
        f"  {label}: {path}"
        for label, path in required.items()
        if not os.path.exists(path)
    ]

    if missing:
        sys.exit(
            "ERROR — the following paths do not exist:\n"
            + "\n".join(missing)
        )

    outdir = os.path.dirname(os.path.abspath(OUTPUT_PREFIX))
    if not os.path.isdir(outdir):
        sys.exit(f"ERROR — output directory does not exist: {outdir}")


check_paths()

gl.options.set_option("data_directory", os.path.join(REF_DIR, ""))


ss = gl.Sumstats(
    INPUT_SUMSTATS,
    build=GENOME_BUILD,
    sep=INPUT_SEP,
    other=EXTRA_COLS,
    na_values=NA_VALUES,
    verbose=VERBOSE,
    **COLUMN_MAP,
)


ss.fix_id(fixchrpos=True)
ss.fix_chr(remove=True)
ss.fix_pos(remove=True)
ss.fix_allele(remove=True)
ss.normalize_allele(threads=THREADS)
ss.sort_coordinate()


ss.remove_dup(
    mode="c",
    keep="last",
    keep_col="P",
    keep_ascend=True
)

ss.check_sanity(eaf=(0, 1), p=(0, 1))
ss.check_data_consistency()

ss.fill_data(
    to_fill=["OR", "OR_95L", "OR_95U"],
    overwrite=True
)

if DROP_BETA_SE:
    ss.data.drop(labels=["BETA", "SE"], axis=1, inplace=True)


ss.sort_coordinate()
ss.sort_column()

ss.to_format(
    path=OUTPUT_PREFIX,
    fmt="gwaslab",
    cols=OUTPUT_EXTRA_COLS
)


lead_variants = ss.get_lead(
    windowsizekb=WINDOW_SIZE_KB,
    sig_level=SIG_LEVEL,
    anno=ANNOTATE_LEADS,
    build=GENOME_BUILD,
    source=ANNO_SOURCE,
)

lead_variants.to_csv(
    LEAD_VARIANTS_OUT,
    index=False,
    sep="\t"
)

print(f"Done. Harmonized sumstats: {OUTPUT_PREFIX}")
print(f"Lead variants ({len(lead_variants)}): {LEAD_VARIANTS_OUT}")


ss.plot_mqq(
    mode="mqq",
    title=title,
    sig_level=1e-6,
    anno_sig_level=5e-8,
    anno="GENENAME",
    build="38",
    sig_line=True,
    additional_line=[5e-8],
    additional_line_color=["black"],
    font_family="DejaVu Sans",
    fontsize=11,
    anno_fontsize=12,
    colors=["#000000", "#ABABAB"],
    save=plot,
    save_kwargs={"dpi": 300, "facecolor": "white"}
)
