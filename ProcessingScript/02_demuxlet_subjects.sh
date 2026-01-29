#!/usr/bin/env bash
################################################################################
# Demuxlet / Popscle Pipeline (Batch: DIG_2ndBatch)
#
# Original header mentions 05/10/2021; this instance configured for:
#   RUN_DATE = 10162024
#   RUN_NAME = 2ndBatch
#
# Purpose:
#   Run popscle `dsc-pileup` and `demuxlet` on Cell Ranger BAMs to assign
#   barcodes to donors using an external genotype VCF (GRCh38).
#
# Background (VCF provenance):
#   The genotype VCF was produced by applying, in order:
#     1) hg19 Affymetrix annotation
#     2) SNP flipping
#     3) sorting
#     4) imputation
#     5) Picard liftover hg19 -> hg38
#     6) exon filtering of variants
#
# Notes:
#   - BAMs are produced by `cellranger count` (not `aggr`).
#   - If you need a merged BAM across runs, consider samtools merge.
#
# Outputs (per sample):
#   <WORKING_DIR>/<SAMPLE>/
#     Pileup-<SAMPLE>_<RUN_NAME>_<RUN_DATE>.*
#     Demuxlet-<SAMPLE>_<RUN_NAME>_<RUN_DATE>.best
#     Demuxlet-<SAMPLE>_<RUN_NAME>_<RUN_DATE>.singletons
#     Demuxlet-<SAMPLE>_<RUN_NAME>_<RUN_DATE>.doublets
################################################################################

set -euo pipefail
IFS=$'\n\t'

#######################################
# User configuration
#######################################
RUN_DATE="10162024"  # MMDDYYYY (used in output naming)
RUN_NAME="2ndBatch"

WORKING_DIR="ProcessedData/Demuxlet/${RUN_NAME}"

POPSCLE="software/popscle/bin/popscle"

# Samples (must align index-wise with BAM/VCF/BARCODE arrays)
SAMPLE_NAMES=(
  "Run1"
  "Run2"
  "Run3"
  "Run4"
)

CELLRANGER_BAMS=(
  "/Run1/outs/possorted_genome_bam.bam"
  "/Run2/outs/possorted_genome_bam.bam"
  "/Run3/outs/possorted_genome_bam.bam"
  "/Run4/outs/possorted_genome_bam.bam"
)

# IMPORTANT:
# You are using 2 paired-donor VCFs (shared across 2 samples each).
# Keep the duplication explicit to preserve indexing consistency.
REFERENCE_VCFS=(
  "/vcf_files/Run1.vcf.gz"
  "/vcf_files/Run1.vcf.gz"
  "/vcf_files/Run1.vcf.gz"
  "/vcf_files/Run1.vcf.gz"
)

BARCODE_FILES=(
  "/Run1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  "/Run2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  "/Run3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  "/Run4/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
)

# Parallelization:
#   Set to 1 to run samples in parallel (background jobs) then `wait`.
#   Set to 0 to run sequentially (simpler debugging).
RUN_IN_PARALLEL=10

#######################################
# Logging helpers
#######################################
RUNTIME_FILE="runtimes-Demuxlet_${RUN_NAME}_${RUN_DATE}.txt"

log() {
  local msg="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S')  ${msg}" | tee -a "${WORKING_DIR}/${RUNTIME_FILE}"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_exe() {
  [[ -x "$1" ]] || die "Executable not found or not executable: $1"
}

require_file() {
  [[ -f "$1" ]] || die "File not found: $1"
}

#######################################
# Sanity checks
#######################################
mkdir -p "${WORKING_DIR}"
: > "${WORKING_DIR}/${RUNTIME_FILE}"

require_exe "${POPSCLE}"

# Ensure arrays align
n_samples="${#SAMPLE_NAMES[@]}"
[[ "${#CELLRANGER_BAMS[@]}" -eq "${n_samples}" ]] || die "CELLRANGER_BAMS length != SAMPLE_NAMES length"
[[ "${#REFERENCE_VCFS[@]}"   -eq "${n_samples}" ]] || die "REFERENCE_VCFS length != SAMPLE_NAMES length"
[[ "${#BARCODE_FILES[@]}"    -eq "${n_samples}" ]] || die "BARCODE_FILES length != SAMPLE_NAMES length"

log "Demuxlet pipeline started"
log "RUN_NAME=${RUN_NAME}"
log "RUN_DATE=${RUN_DATE}"
log "WORKING_DIR=${WORKING_DIR}"
log "POPSCLE=${POPSCLE}"
log "RUN_IN_PARALLEL=${RUN_IN_PARALLEL}"
log "N_SAMPLES=${n_samples}"

#######################################
# Prepare barcode lists (plain text)
#######################################
BARCODE_LISTS=()

for i in "${!SAMPLE_NAMES[@]}"; do
  sample="${SAMPLE_NAMES[$i]}"
  barcode_gz="${BARCODE_FILES[$i]}"

  require_file "${barcode_gz}"

  barcode_txt="${WORKING_DIR}/ValidBarcodes-${sample}_${RUN_NAME}_${RUN_DATE}.txt"
  log "Extracting barcodes: ${sample} -> ${barcode_txt}"

  zcat "${barcode_gz}" > "${barcode_txt}"
  BARCODE_LISTS+=("${barcode_txt}")
done

#######################################
# Run popscle per sample
#######################################
run_one_sample() {
  local sample="$1"
  local bam="$2"
  local vcf="$3"
  local barcodes="$4"

  require_file "${bam}"
  require_file "${vcf}"
  require_file "${barcodes}"

  local sample_dir="${WORKING_DIR}/${sample}"
  mkdir -p "${sample_dir}"
  cd "${sample_dir}"

  local pileup_prefix="Pileup-${sample}_${RUN_NAME}_${RUN_DATE}"
  local demux_prefix="Demuxlet-${sample}_${RUN_NAME}_${RUN_DATE}"

  log "START dsc-pileup: ${sample}"
  "${POPSCLE}" dsc-pileup \
    --sam "${bam}" \
    --vcf "${vcf}" \
    --out "${sample_dir}/${pileup_prefix}"
  log "END dsc-pileup: ${sample}"

  log "START demuxlet: ${sample}"
  "${POPSCLE}" demuxlet \
    --plp "${sample_dir}/${pileup_prefix}" \
    --vcf "${vcf}" \
    --field GT \
    --group-list "${barcodes}" \
    --out "${sample_dir}/${demux_prefix}"
  log "END demuxlet: ${sample}"

  cd "${WORKING_DIR}"
}

pids=()

for i in "${!SAMPLE_NAMES[@]}"; do
  sample="${SAMPLE_NAMES[$i]}"
  bam="${CELLRANGER_BAMS[$i]}"
  vcf="${REFERENCE_VCFS[$i]}"
  barcodes="${BARCODE_LISTS[$i]}"

  if [[ "${RUN_IN_PARALLEL}" -eq 1 ]]; then
    run_one_sample "${sample}" "${bam}" "${vcf}" "${barcodes}" &
    pids+=("$!")
  else
    run_one_sample "${sample}" "${bam}" "${vcf}" "${barcodes}"
  fi
done

# Wait for parallel jobs to complete
if [[ "${RUN_IN_PARALLEL}" -eq 1 ]]; then
  log "Waiting for ${#pids[@]} background jobs..."
  for pid in "${pids[@]}"; do
    wait "${pid}"
  done
fi

log "Demuxlet pipeline completed successfully"

