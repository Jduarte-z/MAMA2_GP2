#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import inspect
import io
import platform
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Iterator

import polars as pl

KEY = "__key__"


def open_maybe_gzip(path: Path) -> BinaryIO:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_header_line(fh: BinaryIO, bufsize: int = 1 << 16) -> tuple[bytes, bytes]:
    buf = b""
    while True:
        block = fh.read(bufsize)
        if not block:
            if buf:
                return buf if buf.endswith(b"\n") else buf + b"\n", b""
            raise ValueError("file is empty - no header line found")
        buf += block
        idx = buf.find(b"\n")
        if idx != -1:
            return buf[: idx + 1], buf[idx + 1 :]


def iter_line_aligned_chunks(
    fh: BinaryIO, leftover: bytes, chunk_bytes: int
) -> Iterator[bytes]:
    buf = leftover
    while True:
        block = fh.read(chunk_bytes)
        if not block:
            break
        buf += block
        if len(buf) < chunk_bytes:
            continue
        idx = buf.rfind(b"\n")
        if idx == -1:
            continue
        yield buf[: idx + 1]
        buf = buf[idx + 1 :]

    if buf.strip():
        yield buf if buf.endswith(b"\n") else buf + b"\n"


def _empty_string_kwargs() -> dict:
    # Argument name changed in Polars 1.43.
    params = inspect.signature(pl.read_csv).parameters
    if "empty_string_is_null" in params:
        return {"empty_string_is_null": False}
    if "missing_utf8_is_empty_string" in params:
        return {"missing_utf8_is_empty_string": True}
    return {}


_EMPTY_STRING_KWARGS = _empty_string_kwargs()


def read_csv_all_strings(source, **kwargs) -> pl.DataFrame:
    return pl.read_csv(
        source,
        infer_schema_length=0,
        quote_char=None,
        **_EMPTY_STRING_KWARGS,
        **kwargs,
    )


def key_expr(chrom_col: str, pos_col: str, strip_chr: bool) -> pl.Expr:
    chrom = pl.col(chrom_col).str.strip_chars()
    if strip_chr:
        chrom = chrom.str.replace(r"(?i)^chr", "")
    pos = pl.col(pos_col).str.strip_chars()
    return pl.concat_str([chrom, pl.lit(":"), pos]).alias(KEY)


def load_variants(path: Path, strip_chr: bool) -> tuple[pl.DataFrame, int]:
    # Read the entire line as one field so arbitrary whitespace can be parsed below.
    raw = read_csv_all_strings(
        path,
        has_header=False,
        separator="\x1f",
        new_columns=["line"],
    )
    n_lines = raw.height

    parsed = (
        raw.select(
            pl.col("line")
            .str.extract_groups(r"^\s*(?<chrom>\S+)[ \t]+(?<pos>\S+)")
            .alias("g")
        )
        .unnest("g")
        .filter(pl.col("chrom").is_not_null() & pl.col("pos").is_not_null())
    )

    if parsed.height == 0:
        raise ValueError(
            f"no usable CHROM POS lines parsed from {path} - check the delimiter"
        )

    parsed = parsed.with_columns(key_expr("chrom", "pos", strip_chr))
    return parsed.unique(subset=[KEY], keep="first"), n_lines


def subset_one_file(
    sumstats_path: Path,
    variants: pl.DataFrame,
    out_handle: BinaryIO,
    write_header: bool,
    args,
    matched_keys: set[str],
) -> dict:
    t0 = time.time()
    rows_in = rows_out = 0
    header_written = not write_header
    header_cols: list[str] | None = None

    with open_maybe_gzip(sumstats_path) as fh:
        header, leftover = read_header_line(fh)

        for chunk in iter_line_aligned_chunks(
            fh, leftover, args.chunk_mb * 1024 * 1024
        ):
            df = read_csv_all_strings(
                io.BytesIO(header + chunk),
                has_header=True,
                separator=args.sep,
            )

            if header_cols is None:
                header_cols = df.columns
                for col in (args.chrom_col, args.pos_col):
                    if col not in header_cols:
                        raise SystemExit(
                            f"{sumstats_path}: column '{col}' not found. "
                            f"Header is: {', '.join(header_cols)}"
                        )

            rows_in += df.height

            hit = df.with_columns(
                key_expr(args.chrom_col, args.pos_col, args.strip_chr)
            ).join(variants.select(KEY), on=KEY, how="semi")

            if hit.height == 0:
                continue

            matched_keys.update(hit.get_column(KEY).unique().to_list())
            rows_out += hit.height

            hit.drop(KEY).write_csv(
                out_handle,
                separator=args.sep,
                include_header=not header_written,
                quote_style="never",
                null_value="",
            )
            header_written = True

    if not header_written and header_cols is not None:
        out_handle.write(args.sep.join(header_cols).encode() + b"\n")

    return {
        "path": str(sumstats_path),
        "columns": header_cols or [],
        "rows_in": rows_in,
        "rows_out": rows_out,
        "seconds": time.time() - t0,
    }


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Subset GWAS summary statistics to a CHROM/POS variant list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "--sumstats",
        nargs="+",
        required=True,
        type=Path,
        help="One or more summary-statistics files (plain or gzipped, tab-delimited).",
    )
    p.add_argument(
        "--variants",
        required=True,
        type=Path,
        help="Headerless space-delimited variant list: CHROM POS POS.",
    )
    p.add_argument(
        "--out",
        required=True,
        type=Path,
        help=(
            "Output .tsv file for a single input or with --merge; "
            "otherwise a directory that will hold one .tsv per input."
        ),
    )
    p.add_argument(
        "--merge",
        action="store_true",
        help="Concatenate all inputs into a single output file.",
    )
    p.add_argument("--chrom-col", default="CHR", help="Chromosome column name.")
    p.add_argument("--pos-col", default="POS", help="Position column name.")
    p.add_argument("--sep", default="\t", help="Sumstats field separator.")
    p.add_argument(
        "--chunk-mb",
        type=int,
        default=64,
        help="Decompressed bytes read per chunk.",
    )
    p.add_argument("--log", type=Path, default=None, help="Run log path.")
    p.add_argument(
        "--unmatched",
        type=Path,
        default=None,
        help="Path for variants with no match in any input.",
    )
    p.add_argument(
        "--no-strip-chr-prefix",
        dest="strip_chr",
        action="store_false",
        help="Do not strip a leading 'chr' before comparing chromosomes.",
    )
    p.set_defaults(strip_chr=True)

    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    started = datetime.now()
    t_start = time.time()

    for path in [*args.sumstats, args.variants]:
        if not path.exists():
            raise SystemExit(f"input not found: {path}")

    per_file_out = len(args.sumstats) > 1 and not args.merge

    if per_file_out:
        args.out.mkdir(parents=True, exist_ok=True)
        log_path = args.log or args.out / "subset_sumstats.log"
        unmatched_path = args.unmatched or args.out / "subset_sumstats.unmatched.tsv"
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        log_path = args.log or args.out.with_suffix(args.out.suffix + ".log")
        unmatched_path = args.unmatched or args.out.with_suffix(
            args.out.suffix + ".unmatched.tsv"
        )

    variants, n_variant_lines = load_variants(args.variants, args.strip_chr)
    n_unique = variants.height

    matched_keys: set[str] = set()
    file_reports: list[dict] = []
    outputs: list[Path] = []

    if per_file_out:
        for path in args.sumstats:
            stem = path.name
            for suffix in (".gz", ".bgz"):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
            stem = Path(stem).stem

            dest = args.out / f"{stem}.subset.tsv"

            with open(dest, "wb") as out_handle:
                file_reports.append(
                    subset_one_file(
                        path,
                        variants,
                        out_handle,
                        True,
                        args,
                        matched_keys,
                    )
                )

            outputs.append(dest)

    else:
        with open(args.out, "wb") as out_handle:
            for i, path in enumerate(args.sumstats):
                rep = subset_one_file(
                    path,
                    variants,
                    out_handle,
                    i == 0,
                    args,
                    matched_keys,
                )

                if (
                    i > 0
                    and file_reports
                    and rep["columns"] != file_reports[0]["columns"]
                ):
                    raise SystemExit(
                        f"header mismatch: {path} does not match "
                        f"{file_reports[0]['path']}"
                    )

                file_reports.append(rep)

        outputs.append(args.out)

    unmatched = variants.filter(~pl.col(KEY).is_in(list(matched_keys)))

    unmatched.select(
        pl.col("chrom").alias("CHROM"),
        pl.col("pos").alias("POS"),
    ).write_csv(
        unmatched_path,
        separator="\t",
        include_header=True,
    )

    total_in = sum(r["rows_in"] for r in file_reports)
    total_out = sum(r["rows_out"] for r in file_reports)
    n_matched = len(matched_keys)

    with open(log_path, "w") as log:
        w = lambda s="": log.write(s + "\n")

        w("# subset_sumstats.py run log")
        w(f"started              : {started.isoformat(timespec='seconds')}")
        w(f"finished             : {datetime.now().isoformat(timespec='seconds')}")
        w(f"elapsed_seconds      : {time.time() - t_start:.1f}")
        w(f"command              : {' '.join(sys.argv)}")
        w(f"working_directory    : {Path.cwd()}")
        w(f"python               : {platform.python_version()} ({sys.executable})")
        w(f"polars               : {pl.__version__}")
        w(f"host                 : {platform.node()} / {platform.platform()}")
        w()

        w("## parameters")
        w(f"chrom_col            : {args.chrom_col}")
        w(f"pos_col              : {args.pos_col}")
        w(f"separator            : {args.sep!r}")
        w(f"strip_chr_prefix     : {args.strip_chr}")
        w(f"chunk_mb             : {args.chunk_mb}")
        w(f"merge                : {args.merge}")
        w("matching             : CHROM + POS only, alleles ignored, no dedup of hits")
        w()

        w("## variant list")
        w(f"path                 : {args.variants.resolve()}")
        w(f"lines_read           : {n_variant_lines}")
        w(f"unique_chrom_pos     : {n_unique}")
        w(f"duplicate_lines      : {n_variant_lines - n_unique}")
        w()

        w("## inputs")
        for r in file_reports:
            size = Path(r["path"]).stat().st_size
            w(f"- file               : {Path(r['path']).resolve()}")
            w(f"  size_bytes         : {size}")
            w(f"  columns            : {len(r['columns'])} ({', '.join(r['columns'])})")
            w(f"  rows_read          : {r['rows_in']}")
            w(f"  rows_written       : {r['rows_out']}")
            w(f"  seconds            : {r['seconds']:.1f}")
        w()

        w("## outputs")
        for path in outputs:
            w(f"- {path.resolve()} ({path.stat().st_size} bytes)")
        w(f"- {unmatched_path.resolve()} (unmatched variants)")
        w()

        w("## summary")
        w(f"rows_read_total      : {total_in}")
        w(f"rows_written_total   : {total_out}")
        w(f"variants_requested   : {n_unique}")
        w(f"variants_matched     : {n_matched}")
        w(f"variants_unmatched   : {unmatched.height}")

        pct = 100.0 * n_matched / n_unique if n_unique else 0.0
        w(f"percent_matched      : {pct:.2f}")

        extra = total_out - n_matched
        w(f"multi_row_positions  : {extra} extra rows beyond one per matched position")

    print(
        f"wrote {total_out} rows from {total_in} read; "
        f"{n_matched}/{n_unique} variants matched "
        f"({unmatched.height} unmatched). Log: {log_path}",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
