#!/usr/bin/env python3

import argparse
import gwaslab as gl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify known and novel GWAS loci using gwaslab."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input gwaslab-formatted summary statistics."
    )
    parser.add_argument(
        "--known",
        required=True,
        help="Known loci file."
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for known and novel loci output files."
    )
    parser.add_argument(
        "--sig-level",
        type=float,
        default=5e-8,
        help="Significance threshold."
    )
    parser.add_argument(
        "--window-kb",
        type=int,
        default=500,
        help="Window size in kb."
    )

    return parser.parse_args()


args = parse_args()


ss = gl.Sumstats(
    args.input,
    build="hg38",
    fmt="gwaslab"
)

ss.fix_chr(remove=True)
ss.fix_pos(remove=True)

res = ss.get_novel(
    known=args.known,
    only_novel=False,
    sig_level=args.sig_level,
    windowsizekb=args.window_kb,
    windowsizekb_for_novel=args.window_kb,
    anno=True
)

res.query("NOVEL == True").to_csv(
    f"{args.output_prefix}.novel.tsv",
    sep="\t",
    index=False
)

res.query("NOVEL == False").to_csv(
    f"{args.output_prefix}.known.tsv",
    sep="\t",
    index=False
)
