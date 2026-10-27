#!/bin/bash

cd /mnt/RetD/workspace/LR/workflow


# ---------------------------
# Help message function
# ---------------------------
show_help() {
    echo "Usage: $0 [--quality sup|hac] [-h|--help]"
    echo ""
    echo "Options:"
    echo "  --quality    Model quality to use for basecalling. Default is 'sup'."
    echo "               Accepted values: 'sup', 'hac'."
    echo "  -h, --help   Show this help message and exit."
    echo ""
    echo "Description:"
    echo "  This script runs a Snakemake pipeline for multiple folders defined in a mapping."
    echo "  It dynamically updates the config.yaml file with appropriate input/output paths,"
    echo "  sample barcodes, and model paths before running the pipeline."
    echo ""
    echo "Note:"
    echo "  The pipeline uses Dorado models at:"
    echo "    /models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
    echo "    /models/dna_r10.4.1_e8.2_400bps_hac@v6.0.0"
}

# ---------------------------
# Folder mappings
# ---------------------------

declare -A folders=(

["dlA"]="/mnt/csf/DNAseq_LR_AS/LR_DNA_LIB1_04_09_2026/no_sample_id/20260904_1301_P2S-00833-B_PBO14682_2ddab0e5/pod5 /mnt/RetD/DNAseq_LR_AS/dlApp barcode07 barcode08 barcode15"
)

declare -A folder_maps=(
  ["dlA"]="barcode07=barcode63 barcode08=barcode32 barcode15=barcode15"
)

# ---------------------------
# Default options
# ---------------------------
QUALITY="sup"
DRY_RUN=0

# ---------------------------
# Argument parsing
# ---------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_help
            exit 0
            ;;
        --quality)
            QUALITY="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=1
            shift 1
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Use -h or --help for usage."
            exit 1
            ;;
    esac
done

# Validate quality
if [[ "$QUALITY" != "sup" && "$QUALITY" != "hac" ]]; then
    echo "Invalid quality: '$QUALITY'. Must be 'sup' or 'hac'."
    exit 1
fi

# Model path
if [[ "$QUALITY" == "hac" ]]; then
    MODEL_PATH="/models/dna_r10.4.1_e8.2_400bps_hac@v6.0.0"
else
    MODEL_PATH="/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
fi

# Config file path
CONFIG_FILE="/mnt/RetD/workspace/LR/config/config.yaml"

# Build snakemake args
SNAKEMAKE_ARGS=(--use-singularity --cores --resources gpu_slot=1)
if [[ "$DRY_RUN" -eq 1 ]]; then
  SNAKEMAKE_ARGS=(-n "${SNAKEMAKE_ARGS[@]}")
fi

# ---------------------------
# Main pipeline loop
# ---------------------------
for folder in "${!folders[@]}"; do
  IFS=' ' read -r indir outdir samples <<< "${folders[$folder]}"

  # Basic fields
  sed -i "s|^    indir: .*|    indir: $indir|"   "$CONFIG_FILE"
  sed -i "s|^    outdir: .*|    outdir: $outdir|" "$CONFIG_FILE"
  sed -i "s|^        model: .*|        model: \"$MODEL_PATH\"|" "$CONFIG_FILE"

  # Build folder-local map (optional)
  declare -A MAP=()
  if [[ -n "${folder_maps[$folder]+set}" && -n "${folder_maps[$folder]}" ]]; then
    for kv in ${folder_maps[$folder]}; do
      b="${kv%%=*}"; n="${kv#*=}"
      MAP["$b"]="$n"
    done
  else
    echo "Note: no folder_maps entry for '$folder' — using identity mapping (barcode == sample)." >&2
  fi

  # Produce YAML fragments (identity fallback for any unmapped barcode)
  TMPDIR="$(mktemp -d)"
  SAMPLES_YAML="$TMPDIR/samples.yaml"
  MAP_YAML="$TMPDIR/map.yaml"

  {
    echo "samples:"
    for bc in $samples; do
      name="${MAP[$bc]:-$bc}"   # fallback to the barcode itself
      echo "  - $name"
    done
  } > "$SAMPLES_YAML"

  {
    echo "barcode_to_sample:"
    for bc in $samples; do
      name="${MAP[$bc]:-$bc}"   # fallback to the barcode itself
      echo "  $bc: $name"
    done
  } > "$MAP_YAML"

  # Replace existing 'samples:' and 'barcode_to_sample:' blocks in CONFIG_FILE,
  # then append our fresh ones at the end.
  sed -i -e '/^samples:/,/^$/d' "$CONFIG_FILE" || true
  sed -i -e '/^barcode_to_sample:/,/^$/d' "$CONFIG_FILE" || true
  {
    cat "$SAMPLES_YAML"
    echo ""
    cat "$MAP_YAML"
    echo ""
  } >> "$CONFIG_FILE"

  rm -rf "$TMPDIR"

  # Ensure outdir exists (harmless if it already does)
  mkdir -p "$outdir" 2>/dev/null || true

  # Run Snakemake (sequential; safe since CONFIG_FILE is shared)
  STATS_FILE="${outdir}/stats.json"
  echo "Running Snakemake for $folder... (quality=$QUALITY, dry_run=$DRY_RUN)"
  snakemake "${SNAKEMAKE_ARGS[@]}" --stats "$STATS_FILE" || {
    if [[ "$DRY_RUN" -eq 1 ]]; then
      echo "Dry-run reported an issue for $folder. Skipping..."
      continue
    fi
    echo "Snakemake failed for $folder. Retrying..."
    snakemake "${SNAKEMAKE_ARGS[@]}" --stats "$STATS_FILE" || {
      echo "Snakemake failed again for $folder. Skipping..."
      continue
    }
  }
  echo "Finished processing $folder."
done

echo "All folders processed !"
