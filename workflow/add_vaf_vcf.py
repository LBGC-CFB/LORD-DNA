import sys
import os
import pandas as pd

# Snakemake input/output files
input_vcf = sys.argv[1]
bed_file = sys.argv[2]  # Excel file containing mapping data
output_vcf = sys.argv[3]
input_spip = sys.argv[-1] if len(sys.argv) == 5 else None

def merge_spip_annotations(df, spip_file):
    """Merge Interpretation and InterConfident columns from SPIP file based on CHROM, POS, REF, ALT."""
    if not os.path.exists(spip_file):
        print(f"⚠️  SPIP file not found: {spip_file}")
        df["Interpretation"] = None
        df["InterConfident"] = None
        return df

    spip_df = pd.read_csv(spip_file, sep="\t", dtype={"CHROM": str, "POS": str, "REF": str, "ALT": str})
    spip_df["ID"] = spip_df["CHROM"] + ":" + spip_df["POS"] + "_" + spip_df["REF"] + "_" + spip_df["ALT"]

    for col in ["Interpretation", "InterConfident"]:
        if col not in spip_df.columns:
            spip_df[col] = None

    return df.merge(spip_df[["ID", "Interpretation", "InterConfident"]], on="ID", how="left")


def is_missing(value):
    if value is None:
        return True
    value = str(value).strip()
    return value == "" or value == "." or value.lower() in {"nan", "none", "na"}


def format_float(value):
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text if text else "0"


def compute_af_from_format(format_keys, sample_values):
    """
    Return AF as a comma-separated string.

    Caller behavior differs:
    - Clair3 writes FORMAT/AF.
    - DeepVariant writes FORMAT/VAF.
    - Otherwise use AD/DP, or DV/(DR+DV) when SV read counts are available.
    """
    fmt = {k: v for k, v in zip(format_keys, sample_values)}

    for key in ("AF", "VAF"):
        value = fmt.get(key)
        if not is_missing(value):
            return value

    ad_value = fmt.get("AD")
    if is_missing(ad_value):
        try:
            ref_depth, alt_depth = float(fmt["DR"]), float(fmt["DV"])
        except (KeyError, ValueError, TypeError):
            return None
        # cuteSV uses -1 for unknown reference coverage without genotyping.
        if ref_depth >= 0 and alt_depth >= 0 and ref_depth + alt_depth > 0:
            return format_float(alt_depth / (ref_depth + alt_depth))
        return None

    try:
        depths = [float(x) for x in str(ad_value).split(",") if not is_missing(x)]
    except ValueError:
        return None

    if len(depths) < 2:
        return None

    dp_value = fmt.get("DP")
    try:
        dp = float(dp_value) if not is_missing(dp_value) else sum(depths)
    except ValueError:
        dp = sum(depths)

    if dp <= 0:
        return None

    return ",".join(format_float(alt_depth / dp) for alt_depth in depths[1:])


# Read BED file and validate first line
with open(bed_file, 'r') as f:
    first_line = f.readline().strip()
    if not first_line.startswith("#"):
        raise ValueError("BED file's first line must start with '#' and contain column names.")
    column_names = first_line[1:].split("\t")

# Load BED file into DataFrame
df = pd.read_csv(bed_file, sep="\t", comment="#", header=None, names=column_names, dtype=str)
region_to_genes = {
    str(r).strip(): ",".join(x.strip() for x in str(g).split(",") if x.strip())
    for r, g in zip(df["regionID"], df["controlled_genes"])
    if pd.notna(r) and pd.notna(g) and str(r).strip() and str(g).strip().lower() not in {"", "nan"}
}

# Parse input VCF and collect variant information
vcf_lines = []
vcf_data = []
with open(input_vcf, 'r') as infile:
    for line in infile:
        vcf_lines.append(line)
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 10:
            continue
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        variant_id = f"{chrom}:{pos}_{ref}_{alt}"
        vcf_data.append({"line": line, "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, "ID": variant_id})

# Convert to DataFrame and merge SPIP annotations
df_vcf = pd.DataFrame(vcf_data, columns=["line", "CHROM", "POS", "REF", "ALT", "ID"])
df_vcf = merge_spip_annotations(df_vcf, input_spip) if input_spip else df_vcf
df_unique = df_vcf.drop_duplicates(subset="ID", keep="first")
spip_annot_dict = df_unique.set_index("ID")[["Interpretation", "InterConfident"]].to_dict(orient="index") if input_spip else {}

# Process VCF and write to output
info_headers = {
    "GT": '##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
    "GQ": '##INFO=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n',
    "DP": '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
    "AD": '##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths (ref,alt)">\n',
    "DR": '##INFO=<ID=DR,Number=1,Type=Integer,Description="Reference-supporting reads">\n',
    "DV": '##INFO=<ID=DV,Number=1,Type=Integer,Description="Variant-supporting reads">\n',
    "AF": '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency from caller or read counts">\n',
    "controlled_genes": '##INFO=<ID=controlled_genes,Number=.,Type=String,Description="Genes controlled by this region">\n',
    "Interpretation": '##INFO=<ID=Interpretation,Number=1,Type=String,Description="Interpretation from SPIP file">\n',
    "InterConfident": '##INFO=<ID=InterConfident,Number=1,Type=String,Description="Confidence level from SPIP file">\n',
}
existing_info_headers = set()
with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        if line.startswith("##INFO=<ID="):
            existing_info_headers.add(line.split("ID=", 1)[1].split(",", 1)[0])
        if line.startswith("#CHROM"):
            for key, header in info_headers.items():
                if key not in existing_info_headers:
                    outfile.write(header)

        if line.startswith("#"):
            outfile.write(line)
            continue

        fields = line.strip().split("\t")
        if len(fields) < 10:
            outfile.write(line)
            continue

        info_col = fields[7]
        info_fields = info_col.split(";")
        controlled_genes = set()
        for field in info_fields:
            if field.startswith("greendb_id="):
                region_ids = field.split("=")[1].split(",")
                for region_id in region_ids:
                    gene = region_to_genes.get(region_id)
                    if gene:
                        controlled_genes.add(gene)

        gene_info = f"controlled_genes={','.join(sorted(controlled_genes))}" if controlled_genes else ""

        vep_fields = info_col.split("|")
        if len(vep_fields) > 10 and any(p in vep_fields[10] for p in [":c.", ":n."]):
            vep_fields[10] = vep_fields[10].split(":")[-1]
        updated_info = "|".join(vep_fields)

        # Extract genotype information
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        extracted_info = {k: v for k, v in zip(format_keys, sample_values) if k in ["GT", "GQ", "DP", "AD", "DR", "DV"]}
        existing_info = dict(field.split("=", 1) for field in info_fields if "=" in field)
        # Preserve an AF already supplied by the caller, including cuteSV.
        af_value = existing_info.get("AF")
        if is_missing(af_value):
            af_value = compute_af_from_format(format_keys, sample_values)
        if af_value is not None:
            extracted_info["AF"] = af_value
        updated_info = ";".join(
            field for field in updated_info.split(";")
            if field not in {"", "."} and field.split("=", 1)[0] not in extracted_info
        )
        additional_info = ";".join(f"{k}={v}" for k, v in extracted_info.items())



        # Build new INFO field
        info_parts = [updated_info, additional_info, gene_info]

        if input_spip:
            variant_id = f"{fields[0]}:{fields[1]}_{fields[3]}_{fields[4]}"
            spip_info = spip_annot_dict.get(variant_id, {})
            spip_fields = [f"{k}={v}" for k, v in spip_info.items() if pd.notnull(v)]
            info_parts.extend(spip_fields)

        fields[7] = ";".join(part for part in info_parts if part)
        outfile.write("\t".join(fields) + "\n")
