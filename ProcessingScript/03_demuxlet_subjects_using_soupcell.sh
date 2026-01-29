#!/usr/bin/env bash
################################################################################
# souporcell Pipeline (Genotype-free demultiplexing)
# DIG CITE-seq / scRNA-seq Project 2022
#
# Purpose:
#   Run souporcell to infer donor clusters (and estimate ambient RNA)
#   WITHOUT an external genotype VCF.
#
# What this does (per sample):
#   1) Prepare/locate Cell Ranger BAM
#   2) Run souporcell pipeline (container recommended)
#   3) Outputs cluster assignments per barcode + summary files
#
# Requirements:
#   - A BAM file from `cellranger count` (possorted_genome_bam.bam)
#   - Reference genome FASTA (same build as BAM, typically GRCh38)
#   - barcodes.tsv.gz (optional but recommended)
#   - souporcell installed OR docker/singularity image available
#
# Notes:
#   - Choose K = expected number of donors in the pooled sample
#   - Recommended to run via docker/singularity for reproducibility
################################################################################

set -euo pipefail
IFS=$'\n\t'

#######################################
# User configuration
#######################################
RUN_DATE="2024-10-16"
RUN_NAME="2ndBatch"

WORKING_DIR="/ProcessedData/Souporcell/${RUN_NAME}"
mkdir -p "${WORKING_DIR}"

# Expected number of donors in the pool
K=2

# Threads to use
THREADS=20

# Reference FASTA (must match BAM build; example GRCh38)
REFERENCE_FASTA="/References/GRCh38/GRCh38.fa"

# Optional: known barcodes list (recommended). If provided, souporcell restricts to these.
# If you only have barcodes.tsv.gz, you can pass that file directly.
# Some souporcell versions accept --barcodes <file>, others use --cell-barcodes.
# We'll generate a plain text barcodes file to be safe.
MAKE_BARCODES_LIST=1

# Choose runtime mode:
#   MODE="docker"   (recommended)
#   MODE="native"   (if you installed souporcell locally)
MODE="docker"

# Docker image (common public image; change if you use another)
SOUPORCELL_IMAGE="Batch_CITEseq/souporcell:latest"

# Native paths (used only if MODE="native")
SOUPORCELL_PIPELINE="/path/to/souporcell_pipeline.py"  # e.g., ~/souporcell/souporcell_pipeline.py
SOUPORCELL_ENV_PYTHON="python3"

# Samples (index-aligned)
SAMPLE_NAMES=(
  "Run1"
  "Run2"
)

CELLRANGER_BAMS=(
  "/Run1/outs/possorted_genome_bam.bam"
  "/Run2/outs/possorted_genome_bam.bam"
)

BARCODES_TSV_GZ=(
  "/Run1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  "/Run2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
)

#######################################
# Logging helpers
#######################################
RUNTIME_FILE="${WORKING_DIR}/runtimes-Souporcell_${RUN_NAME}_${RUN_DATE}.txt"

log() {
  local msg="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S')  ${msg}" | tee -a "${RUNTIME_FILE}"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || die "File not found: $1"
}

require_exe() {
  command -v "$1" >/dev/null 2>&1 || die "Executable not found in PATH: $1"
}

#######################################
# Checks
#######################################
: > "${RUNTIME_FILE}"

n="${#SAMPLE_NAMES[@]}"
[[ "${#CELLRANGER_BAMS[@]}" -eq "${n}" ]] || die "CELLRANGER_BAMS length != SAMPLE_NAMES length"
[[ "${#BARCODES_TSV_GZ[@]}" -eq "${n}" ]] || die "BARCODES_TSV_GZ length != SAMPLE_NAMES length"

require_file "${REFERENCE_FASTA}"

log "Souporcell pipeline started"
log "RUN_NAME=${RUN_NAME}"
log "RUN_DATE=${RUN_DATE}"
log "WORKING_DIR=${WORKING_DIR}"
log "K=${K}"
log "THREADS=${THREADS}"
log "REFERENCE_FASTA=${REFERENCE_FASTA}"
log "MODE=${MODE}"

if [[ "${MODE}" == "docker" ]]; then
  require_exe docker
elif [[ "${MODE}" == "native" ]]; then
  require_file "${SOUPORCELL_PIPELINE}"
  require_exe "${SOUPORCELL_ENV_PYTHON}"
else
  die "MODE must be 'docker' or 'native'"
fi

#######################################
# Run souporcell per sample
#######################################
for i in "${!SAMPLE_NAMES[@]}"; do
  sample="${SAMPLE_NAMES[$i]}"
  bam="${CELLRANGER_BAMS[$i]}"
  bc_gz="${BARCODES_TSV_GZ[$i]}"

  require_file "${bam}"
  require_file "${bc_gz}"

  sample_dir="${WORKING_DIR}/${sample}"
  mkdir -p "${sample_dir}"

  # Prepare barcodes list
  barcodes_txt="${sample_dir}/barcodes_${sample}.txt"
  if [[ "${MAKE_BARCODES_LIST}" -eq 1 ]]; then
    log "Extracting barcodes for ${sample} -> ${barcodes_txt}"
    zcat "${bc_gz}" > "${barcodes_txt}"
  fi

  log "START souporcell: ${sample}"
  log "BAM=${bam}"

  if [[ "${MODE}" == "docker" ]]; then
    # Mount everything we need into the container. Adjust mounts for your environment.
    docker run --rm \
      -u "$(id -u):$(id -g)" \
      -v "${sample_dir}:/out" \
      -v "$(dirname "${bam}"):/bamdir:ro" \
      -v "$(dirname "${REFERENCE_FASTA}"):/ref:ro" \
      "${SOUPORCELL_IMAGE}" \
      souporcell_pipeline.py \
        -i "/bamdir/$(basename "${bam}")" \
        -b "/out/$(basename "${barcodes_txt}")" \
        -f "/ref/$(basename "${REFERENCE_FASTA}")" \
        -t "${THREADS}" \
        -k "${K}" \
        -o "/out" \
        --skip_remap True

  else
    # Native run
    "${SOUPORCELL_ENV_PYTHON}" "${SOUPORCELL_PIPELINE}" \
      -i "${bam}" \
      -b "${barcodes_txt}" \
      -f "${REFERENCE_FASTA}" \
      -t "${THREADS}" \
      -k "${K}" \
      -o "${sample_dir}" \
      --skip_remap True
  fi

  log "END souporcell: ${sample}"
done

log "Souporcell pipeline completed successfully"

