#!/usr/bin/env python3
import pandas as pd
import sys
import os
import re
from glob import glob

# -------------------------
# Helpers: deterministic ordering (natural sort)
# -------------------------
_natsplit = re.compile(r"(\d+)")


def natural_key(s: str):
    """
    Natural sort key: "sample2" < "sample10"
    Splits digits and compares numerically, case-insensitive.
    """
    return [int(t) if t.isdigit() else t.lower() for t in _natsplit.split(s)]


# Read arguments
bed_file = sys.argv[1]
depth_folder = sys.argv[2]  # Folder containing all depth files (*.depth)
prefix = sys.argv[3]

# Load BED (already filtered)
bed = pd.read_csv(
    bed_file,
    sep=r"\s+",
    header=None,
    names=["chrom", "start", "end", "gene_exon", "score", "strand"],
)

# Ensure numeric coordinates
bed["start"] = bed["start"].astype(int)
bed["end"] = bed["end"].astype(int)

# Extract gene name (everything before "_" or ";")
bed["gene"] = bed["gene_exon"].str.extract(r"^([^\_;]+)")

# Region length (BED: end is exclusive, so end - start)
bed["region_length"] = bed["end"] - bed["start"]

# Precompute total target length per gene (sum of region lengths)
gene_lengths = bed.groupby("gene")["region_length"].sum().to_dict()

# Group BED by gene to iterate regions
grouped_bed = bed.groupby("gene", sort=True)  # deterministic gene iteration

# Prepare aggregated files
agg_30X = {}
agg_20X = {}
agg_10X = {}

# Process each depth file (IMPORTANT: sort files by sample name)
depth_files = glob(os.path.join(depth_folder, "*.depth"))

def sample_from_path(p: str) -> str:
    return os.path.basename(p).replace(".depth", "")

depth_files = sorted(depth_files, key=lambda p: natural_key(sample_from_path(p)))

for depth_file in depth_files:
    sample_name = sample_from_path(depth_file)
    output_file = os.path.join(depth_folder, f"{sample_name}_gene_coverage.tsv")

    depth = pd.read_csv(depth_file, sep="\t", header=None, names=["chrom", "pos", "depth"])

    sample_30X = {}
    sample_20X = {}
    sample_10X = {}

    with open(output_file, "w", encoding="utf-8", newline="\n") as out:
        out.write("Gene\tlength\t%_30X\t%_20X\t%_10X\n")

        for gene, regions in grouped_bed:
            gene_depths = []

            for _, row in regions.iterrows():
                chrom = row["chrom"]
                start = int(row["start"])
                end = int(row["end"])

                region_depth = depth[
                    (depth["chrom"] == chrom)
                    & (depth["pos"] >= start)
                    & (depth["pos"] <= end)
                ]
                gene_depths.append(region_depth)

            if not gene_depths:
                pct_30 = pct_20 = pct_10 = 0.00
            else:
                all_depths = pd.concat(gene_depths, ignore_index=True)
                if all_depths.empty:
                    pct_30 = pct_20 = pct_10 = 0.00
                else:
                    total = len(all_depths)
                    pct_30 = round(100 * (all_depths["depth"] >= 30).sum() / total, 2)
                    pct_20 = round(100 * (all_depths["depth"] >= 20).sum() / total, 2)
                    pct_10 = round(100 * (all_depths["depth"] >= 10).sum() / total, 2)

            length_bp = int(gene_lengths.get(gene, 0))

            out.write(f"{gene}\t{length_bp}\t{pct_30:.2f}\t{pct_20:.2f}\t{pct_10:.2f}\n")

            sample_30X[gene] = pct_30
            sample_20X[gene] = pct_20
            sample_10X[gene] = pct_10

    # dict insertion order now matches sorted depth_files
    agg_30X[sample_name] = sample_30X
    agg_20X[sample_name] = sample_20X
    agg_10X[sample_name] = sample_10X


# Save merged outputs
def save_aggregated(file_path, agg_dict, threshold):
    if not agg_dict:
        raise SystemExit(f"No samples found in {depth_folder}/*.depth")

    all_genes = sorted(set().union(*[d.keys() for d in agg_dict.values()]), key=natural_key)

    label = f"all {prefix} - % coverage @ {threshold}X"
    description = f"Coverage % for all {prefix} @ {threshold}X"
    section_id = f"coverage_{prefix}_{threshold}_all"

    # IMPORTANT: enforce sorted sample output regardless of dict order
    ordered_samples = sorted(agg_dict.keys(), key=natural_key)

    with open(file_path, "w", encoding="utf-8", newline="\n") as f:
        # MultiQC metadata header
        f.write("# parent_id: coverage\n")
        f.write("# parent_name: 'Coverage Analysis'\n")
        f.write("# parent_description: 'Proportion of captured sequence which reach one or more coverage value'\n")
        f.write(f"# id: {section_id}\n")
        f.write(f"# section_name: '{label}'\n")
        f.write(f"# description: '{description}'\n")
        f.write("# format: 'tsv'\n")
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write(f"#    id: 'custom_bargraph_{prefix}_{threshold}_all'\n")
        f.write(f"#    title: '% of all {prefix} >= {threshold}X'\n")
        f.write("#    format: '{:,.3f}'\n")
        f.write("#    scale: 'RdYlGn'\n")

        # Column headers: samples × genes
        f.write("Sample\t" + "\t".join(all_genes) + "\n")

        for sample in ordered_samples:
            gene_dict = agg_dict[sample]
            values = [f"{gene_dict.get(gene, 0.00):.2f}" for gene in all_genes]
            f.write(sample + "\t" + "\t".join(values) + "\n")


save_aggregated(os.path.join(depth_folder, f"{prefix}_30X_mqc.tsv"), agg_30X, threshold=30)
save_aggregated(os.path.join(depth_folder, f"{prefix}_20X_mqc.tsv"), agg_20X, threshold=20)
save_aggregated(os.path.join(depth_folder, f"{prefix}_10X_mqc.tsv"), agg_10X, threshold=10)
