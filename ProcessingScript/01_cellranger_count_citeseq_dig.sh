#!/usr/bin/env bash
################################################################################
# Cell Ranger count (Feature Barcoding / CITE-seq)
# DIG CITE-seq Project 2022
#
# Author: Sugandh
# Created: 2024-10-10 (YYYY-MM-DD)
#
# Purpose:
#   Run `cellranger count` for multiple libraries using per-library
#   "Input-CellRangerLibraries-<LIB>-<DATE>.csv" files.
#
# Notes:
#   - Designed for batch/HPC usage. Consider running inside screen/tmux.
#   - Requires: cellranger >= 8, reference transcriptome, library CSV files.
################################################################################

set -euo pipefail
IFS=$'\n\t'

#######################################
# User config
#######################################
RUN_DATE="10102024"                 # MMDDYYYY (used in input CSV filename)
RUN_NAME="240924"
CORES_TO_USE=30

# Where your per-library Input-CellRangerLibraries-*.csv files live
WORKING_DIR="Processing/${RUN_NAME}"

# Cell Ranger binary
CELLRANGER="cellranger-8.0.1/cellranger"

# Transcriptome reference
REFERENCE_10X_DB="refdata-cellranger-GRCh38-1.2.0/"

# Output base directory (Cell Ranger will create per-library folders under this)
PARENT_OUTPUT_DIR="ProcessedData/CellRangerOutput_${RUN_NAME}"
COUNT_OUTPUT_DIR="${PARENT_OUTPUT_DIR}"

# Optional: Feature reference (not used by `--libraries` CSV; included for documentation)
ADT_FEATURE_REFERENCE="TotalSeqA-v01r00a_Reference.csv"

# Libraries to process
INPUT_LIBRARIES=(
  "DIG-13L-14NL"
  "TIG-DIG2_L2"
  "TIG-DIG2_L3"
  "TIG-DIG2_L4"
)

#######################################
# Helpers
#######################################
RUNTIMES_FILE="runtimes-Processing_RNASeq_${RUN_NAME}_${RUN_DATE}.txt"

log() {
  # Timestamped log line to stdout and runtimes file
  local msg="$1"
  echo "$(date '+%Y-%m-%d %H:%M:%S')  ${msg}" | tee -a "${WORKING_DIR}/${RUNTIMES_FILE}"
}

die() {
  echo "ERROR: $*" 1>&2
  exit 1
}

require_exe() {
  local exe="$1"
  [[ -x "$exe" ]] || die "Executable not found or not executable: ${exe}"
}

#######################################
# Setup / checks
#######################################
mkdir -p "${WORKING_DIR}" "${COUNT_OUTPUT_DIR}"

require_exe "${CELLRANGER}"
[[ -d "${REFERENCE_10X_DB}" ]] || die "Reference transcriptome folder not found: ${REFERENCE_10X_DB}"

log "Run started"
log "RUN_NAME=${RUN_NAME}"
log "RUN_DATE=${RUN_DATE}"
log "CORES_TO_USE=${CORES_TO_USE}"
log "WORKING_DIR=${WORKING_DIR}"
log "COUNT_OUTPUT_DIR=${COUNT_OUTPUT_DIR}"
log "CELLRANGER=${CELLRANGER}"
log "REFERENCE_10X_DB=${REFERENCE_10X_DB}"
log "ADT_FEATURE_REFERENCE=${ADT_FEATURE_REFERENCE} (informational)"

#######################################
# Main: run cellranger count per library
#######################################
for current_library in "${INPUT_LIBRARIES[@]}"; do
  # Expect exactly one matching CSV
  csv_glob="${WORKING_DIR}/Input-CellRangerLibraries-${current_library}-${RUN_DATE}.csv"

  # If your naming scheme differs, you can fall back to find:
  # current_library_file="$(find "${WORKING_DIR}" -maxdepth 1 -type f -name "Input-CellRangerLibraries-${current_library}-${RUN_DATE}.csv" | head -n 1)"

  [[ -f "${csv_glob}" ]] || die "Library CSV not found: ${csv_glob}"
  current_library_file="${csv_glob}"

  # Output folder
  mkdir -p "${COUNT_OUTPUT_DIR}/${current_library}"
  cd "${COUNT_OUTPUT_DIR}"

  log "START cellranger count: ${current_library}"
  log "Using libraries CSV: ${current_library_file}"

  "${CELLRANGER}" count \
    --id="${current_library}" \
    --libraries="${current_library_file}" \
    --transcriptome="${REFERENCE_10X_DB}" \
    --create-bam=true \
    --nosecondary \
    --localcores="${CORES_TO_USE}"

  log "END cellranger count: ${current_library}"
done

#######################################
# Cleanup: remove large intermediate folders
#######################################
# WARNING: Removing SC_RNA_COUNTER_CS can break re-runs / troubleshooting.
# Keep if you need detailed logs or want to re-run without recomputing.
for current_library in "${INPUT_LIBRARIES[@]}"; do
  target_dir="${COUNT_OUTPUT_DIR}/${current_library}/SC_RNA_COUNTER_CS"
  if [[ -d "${target_dir}" ]]; then
    log "Cleanup: removing ${target_dir}"
    rm -rf "${target_dir}"
  else
    log "Cleanup: not found (skipping): ${target_dir}"
  fi
done

log "Run finished successfully"

