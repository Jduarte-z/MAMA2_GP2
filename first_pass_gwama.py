#!/usr/bin/env python3

# GWAS summary statistics harmonization with gwaslab.
# Designed and tested on gwaslab v4.1.6.
# Edit the CONFIG section below.

import os
import sys
import gwaslab as gl

# CONFIG
#Reference resources
#they are usually downloaded from gwaslab
"""
in the past I used these commands in the folder I chose as my reference directory.
by default gwaslab puts them in ~/.gwaslab/


in the desired directory run this in an interactive python session:

python
import gwaslab as gl
#for the 1KGP references:
gl.download_ref('1kg_eur_hg38', directory='./')
gl.download_ref('1kg_amr_hg38', directory='./')
gl.download_ref('1kg_sas_hg38', directory='./')
gl.download_ref('1kg_eas_hg38', directory='./')
gl.download_ref('1kg_afr_hg38', directory='./')
gl.download_ref('1kg_pan_hg38', directory='./')
#the pan is the one for everyone merged

#for the fasta file:
gl.download_ref('ucsc_genome_hg38', directory='./')

"""
REF_DIR = "/home/jupyter/workspace/ws_files/MAMA_2/gwaslab_referenceFiles"


REF_FASTA = f"{REF_DIR}/hg38.fa"


#Change accordignly to the population in question, the only options available are
#AMR, EUR, AFR, EAS and SAS. From the 1KGP, so pick the closest one,
#at the end it is just for infering the strand of palindromics and indels
REF_VCF = f"{REF_DIR}/PAN.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"

#Input / output
#CHANGE ACCORDIGNLY
INPUT_SUMSTATS = "gwama_random.out"
OUTPUT_PREFIX = "gwmam_random_gp2Only_non-eur"
LEAD_VARIANTS_OUT = f"{OUTPUT_PREFIX}.lead_variants.tsv"

plot="manhattan.png"
title="gp2 non-eur random effects meta"

#Run parameters
THREADS = 4
GENOME_BUILD = "38"
INPUT_SEP = "\t"
NA_VALUES = ["NA", ".", "nan", ""]
VERBOSE = True


#Input column names
COLUMN_MAP = {
    "snpid": "MarkerName",
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
    "mlog10p": "-log10_p-value"
}
EXTRA_COLS = ["OR_se", "q_statistic", "q_p-value"]
#EA-referenced frequency columns gwaslab does not know about and will not flip
#EA_FREQ_COLS = ["A1FREQ_CASES", "A1FREQ_CONTROLS"]
#Columns to carry into the output file; the gwaslab format drops anything it doesn't recognize
OUTPUT_EXTRA_COLS = EXTRA_COLS

#Strand inference (infer_strand)
# remove_snp="78": drop palindromic SNPs whose strand cannot be inferred from MAF (7) and those absent from the reference VCF (8)
# remove_indel="8": drop indels absent from the reference VCF
# REF_ALT_FREQ_FIELD = "AF"
# REMOVE_SNP = "78"
# REMOVE_INDEL = "8"
# MAF_THRESHOLD = 0.40
# DAF_TOLERANCE = 0.20

#Lead variant extraction
WINDOW_SIZE_KB = 500
SIG_LEVEL = 5e-8
ANNOTATE_LEADS = True
ANNO_SOURCE = "ensembl"

#Optional: drop BETA/SE from the final output
DROP_BETA_SE = False


# END CONFIG — no edits needed below

def check_paths():
    """Fail fast if any required file is missing."""
    required = {
        "Input summary statistics": INPUT_SUMSTATS,
        "Reference FASTA": REF_FASTA,
        "Reference VCF": REF_VCF,
        "gwaslab data directory": REF_DIR,
    }
    missing = [f"  {label}: {path}" for label, path in required.items()
               if not os.path.exists(path)]
    if missing:
        sys.exit("ERROR — the following paths do not exist:\n" + "\n".join(missing))

    outdir = os.path.dirname(os.path.abspath(OUTPUT_PREFIX))
    if not os.path.isdir(outdir):
        sys.exit(f"ERROR — output directory does not exist: {outdir}")


# def post_alignment_qc(ss):
#     """Flip stats, rebuild SNPID from scratch, renormalize alleles.

#     Run after any step that may change allele orientation
#     (check_ref, infer_strand).
#     """
#     # EAF is the only frequency column flip_allele_stats() modifies, so watching it
#     # identifies the flipped rows exactly. 1 - x is invisible at EAF == 0.5 and
#     # undefined at NaN, so park those on a sentinel for the duration of the call.
#     extra_freq = [c for c in EA_FREQ_COLS if c in ss.data.columns]
#     parked = ss.data["EAF"].isna() | ss.data["EAF"].eq(0.5)
#     parked_values = ss.data.loc[parked, "EAF"].copy()
#     #just an random number place holder to substitute the sentinel values (impossible ones of 0.5 and NAs) for the duration of the call
#     ss.data.loc[parked, "EAF"] = 0.25

#     eaf_before = ss.data["EAF"].to_numpy(copy=True)
#     ss.flip_allele_stats()
#     flipped = ss.data["EAF"].to_numpy() != eaf_before

#     ss.data.loc[parked, "EAF"] = parked_values
#     for col in extra_freq:
#         ss.data.loc[flipped, col] = 1 - ss.data.loc[flipped, col]

#     ss.data.drop(labels=["SNPID"], axis=1, inplace=True)
#     ss.fix_id(fixchrpos=False, fixid=True, fixsep=False,
#               forcefixid=True, overwrite=True)
#     ss.normalize_allele(threads=THREADS)


check_paths()

# gwaslab expects a trailing separator on data_directory
gl.options.set_option("data_directory", os.path.join(REF_DIR, ""))


# Load

ss = gl.Sumstats(
    INPUT_SUMSTATS,
    build=GENOME_BUILD,
    sep=INPUT_SEP,
    other=EXTRA_COLS,
    na_values=NA_VALUES,
    verbose=VERBOSE,
    **COLUMN_MAP,
)


# Strict data cleaning

ss.fix_chr(remove=True)
ss.fix_pos(remove=True)
ss.fix_allele(remove=True)
ss.fix_id(fixchrpos=False, fixid=True, fixsep=False,
          forcefixid=True, overwrite=True)
ss.normalize_allele(threads=THREADS)
ss.sort_coordinate()


# Harmonization: align NEA to REF in the reference FASTA

# ss.check_ref(ref_seq=REF_FASTA)
# post_alignment_qc(ss)

# Infer strand for palindromic SNPs and indels against the reference VCF
# ss.infer_strand(
#     ref_infer=REF_VCF,
#     ref_alt_freq=REF_ALT_FREQ_FIELD,
#     remove_snp=REMOVE_SNP,
#     remove_indel=REMOVE_INDEL,
#     maf_threshold=MAF_THRESHOLD,
#     daf_tolerance=DAF_TOLERANCE,
#     threads=THREADS,
# )
#re do the SNPID making since some variants will have the EA and NEA flipped and the SNPID will be different
# post_alignment_qc(ss)


# QC before filling columns and writing output

ss.remove_dup(mode="c", keep="last", keep_col="P", keep_ascend=True)
ss.check_sanity(eaf=(0, 1), p=(0, 1))
ss.check_data_consistency()

ss.fill_data(to_fill=["OR", "OR_95L", "OR_95U"], overwrite=True)
if DROP_BETA_SE:
    ss.data.drop(labels=["BETA", "SE"], axis=1, inplace=True)
#get the column for the effective sample size filled as well (they use what is already in place for the number of cases and controls)
# ss.get_ess(method="metal")

# Tidy and write

ss.sort_coordinate()
ss.sort_column()
ss.to_format(path=OUTPUT_PREFIX, fmt="gwaslab", cols=OUTPUT_EXTRA_COLS)

# just a raw visualization of the significant variants
lead_variants = ss.get_lead(
    windowsizekb=WINDOW_SIZE_KB,
    sig_level=SIG_LEVEL,
    anno=ANNOTATE_LEADS,
    build=GENOME_BUILD,
    source=ANNO_SOURCE,
)
lead_variants.to_csv(LEAD_VARIANTS_OUT, index=False, sep="\t")

print(f"Done. Harmonized sumstats: {OUTPUT_PREFIX}")
print(f"Lead variants ({len(lead_variants)}): {LEAD_VARIANTS_OUT}")


ss.plot_mqq(
    mode='mqq',
    #cut=14,
    #skip=5,
    title=title, 
    sig_level=1e-6,
    anno_sig_level=5e-8,
    anno="GENENAME",
    build="38",
    sig_line=True,
    additional_line=[5e-8],
    additional_line_color=['black'],
    font_family="DejaVu Sans",
    fontsize=11,
    anno_fontsize=12,
    colors=["#000000","#ABABAB"],
    save=plot, save_kwargs={"dpi":300, "facecolor":"white"}
)
