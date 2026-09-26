#!/usr/bin/env python3

import pandas as pd


# input files
gwas1_file = "../r12_mama_mr-mega_10pops_g2019sAJcond/fine_map/GenomicRiskLoci.txt"
gwas2_file = "../r12_mama_randomEff_20pops_g2019sAJcond/fuma/FUMA_job774100/GenomicRiskLoci_randomEff.txt"

# settings
distance = 250000
out = "combined_loci"


def chr_order(chrom):
    chrom = str(chrom).replace("chr", "")

gwas1 = pd.read_csv(gwas1_file, sep="\t")
gwas2 = pd.read_csv(gwas2_file, sep="\t")

for df in (gwas1, gwas2):
    df["chr"] = df["chr"].astype(str)
    df["start"] = pd.to_numeric(df["start"])
    df["end"] = pd.to_numeric(df["end"])

gwas1["GWAS"] = "GWAS1"
gwas2["GWAS"] = "GWAS2"

gwas1["original_row"] = gwas1.index + 1
gwas2["original_row"] = gwas2.index + 1

loci = pd.concat([gwas1, gwas2], ignore_index=True)

loci["chr_order"] = loci["chr"].apply(chr_order)

loci = loci.sort_values( ["chr_order", "start", "end"]).reset_index(drop=True)


cluster_ids = []
cluster_id = 0
current_chr = None
current_end = None

for _, row in loci.iterrows():

    chrom = row["chr"]
    start = row["start"]
    end = row["end"]

    if chrom != current_chr:
        cluster_id += 1
        current_chr = chrom
        current_end = end

    elif start <= current_end + distance:
        current_end = max(current_end, end)

    else:
        cluster_id += 1
        current_end = end

    cluster_ids.append(cluster_id)


loci["CombinedRegion"] = cluster_ids

regions = []

for region_id, group in loci.groupby("CombinedRegion"):

    g1 = group[group["GWAS"] == "GWAS1"]
    g2 = group[group["GWAS"] == "GWAS2"]

    if len(g1) and len(g2):
        classification = "SHARED"
    elif len(g1):
        classification = "GWAS1_ONLY"
    else:
        classification = "GWAS2_ONLY"

    regions.append({
        "CombinedRegion": region_id,
        "chr": group["chr"].iloc[0],
        "start": group["start"].min(),
        "end": group["end"].max(),
        "classification": classification,
        "n_GWAS1_loci": len(g1),
        "n_GWAS2_loci": len(g2),
        "GWAS1_GenomicLoci": ";".join(g1["GenomicLocus"].astype(str)),
        "GWAS2_GenomicLoci": ";".join(g2["GenomicLocus"].astype(str)),
        "GWAS1_rsIDs": ";".join(g1["rsID"].astype(str)),
        "GWAS2_rsIDs": ";".join(g2["rsID"].astype(str))
    })


regions = pd.DataFrame(regions)

n_shared = (regions["classification"] == "SHARED").sum()
n_gwas1 = (regions["classification"] == "GWAS1_ONLY").sum()
n_gwas2 = (regions["classification"] == "GWAS2_ONLY").sum()


print(f"GWAS1 loci: {len(gwas1)}")
print(f"GWAS2 loci: {len(gwas2)}")
print(f"Combined regions: {len(regions)}")
print(f"Shared: {n_shared}")
print(f"GWAS1 only: {n_gwas1}")
print(f"GWAS2 only: {n_gwas2}")


regions.to_csv(f"{out}.regions.tsv", sep="\t", index=False)

loci.to_csv(f"{out}.locus_membership.tsv", sep="\t", index=False
)

regions[regions["classification"] == "SHARED"].to_csv(f"{out}.shared.tsv", sep="\t", index=False)

regions[regions["classification"] == "GWAS1_ONLY"].to_csv( f"{out}.gwas1_only.tsv", sep="\t",index=False)

regions[regions["classification"] == "GWAS2_ONLY"].to_csv(f"{out}.gwas2_only.tsv",sep="\t", index=False)
