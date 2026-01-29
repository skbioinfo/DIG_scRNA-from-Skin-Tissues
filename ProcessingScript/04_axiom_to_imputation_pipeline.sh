#!/usr/bin/env bash
################################################################################
# Axiom (.CEL) -> Genotypes/QC -> VCF -> (liftover hg19->GRCh38) -> conform/strand
# -> Michigan Imputation Server submission prep -> post-imputation QC
# -> remove exonic regions (GRCh38)
#
# Author: <your name>
# Created: 2026-01-28
#
# What you must provide:
#   - APT installed (apt-probeset-genotype, apt-geno-qc, etc.)
#   - Axiom array library files (cdf/pgf/clf, annotation, etc.)
#   - Reference FASTA for liftover target build (GRCh38) (and dict)
#   - Chain file hg19ToHg38.over.chain.gz (if liftover needed)
#   - Michigan Imputation Server submission is manual (upload + run)
#
# Outputs:
#   - Genotype calls + QC reports
#   - Pre-imputation VCF (bgzip/tabix) ready for Michigan
#   - Post-imputation filtered VCF (and exonic-removed VCF)
################################################################################

set -euo pipefail
IFS=$'\n\t'

########################################
# 0) CONFIGURATION (EDIT THIS SECTION)
########################################

# Project identifiers
RUN_NAME="Batch2"
OUTDIR="/Genotyping/${RUN_NAME}"
mkdir -p "${OUTDIR}"

# Input CEL directory
CEL_DIR="/Genotyping/CEL/${RUN_NAME}"

# Sample metadata (optional but recommended)
# Can be used later for sex-check, sample IDs, etc.
SAMPLE_SHEET="${OUTDIR}/samples.tsv"   # optional

# APT executables
APT_GENO_QC="/path/to/apt-geno-qc"
APT_PROBESET_GENOTYPE="/path/to/apt-probeset-genotype"
APT_SUMMARY_GENOTYPE="/path/to/apt-summary-genotype"
APT_CHP_TO_VCF="/path/to/apt-format-result"   # array-dependent; see notes below

# Axiom library files (ARRAY-SPECIFIC; ask your core for the correct bundle)
# Commonly includes: *.cdf or *.pgf/*.clf, *.xml, *.arr, annotation files, etc.
AXIOM_LIB_DIR="/path/to/Axiom_LibraryFiles"
# Examples (names vary by array):
#   - *.cdf OR *.pgf + *.clf
#   - AxiomGT1.genotype_model / *.xml / *.arr
#   - *.chrXprobes, *.chrYprobes, etc. (optional)
AXIOM_CDF="${AXIOM_LIB_DIR}/your_array.cdf"      # OR use PGF/CLF
AXIOM_GENO_MODEL="${AXIOM_LIB_DIR}/AxiomGT1.genotype_model"  # if provided

# Genome build your Axiom annotation coordinates are in:
# Set to "hg19" or "GRCh38"
ARRAY_BUILD="hg19"

# Liftover settings (if ARRAY_BUILD=hg19 but you want GRCh38)
DO_LIFTOVER=1
PICARD_JAR="/path/to/picard.jar"
CHAIN_HG19_TO_HG38="/path/to/hg19ToHg38.over.chain.gz"
REF_FASTA_GRCH38="/path/to/GRCh38.fa"          # must have .fai and .dict
# If you also need hg19 reference for sanity checks:
REF_FASTA_HG19="/path/to/hg19.fa"              # optional

# Conform/strand alignment before MIS:
# Michigan recommends: bgzip + tabix, normalized & left-aligned,
# and allele strand alignment to the reference panel.
BCFTOOLS="/usr/bin/bcftools"
TABIX="/usr/bin/tabix"
BGZIP="/usr/bin/bgzip"

# "conform-gt" is commonly used for aligning VCF alleles to panel reference
# https://github.com/BRaVaGenomics/conform-gt
CONFORM_GT_JAR="/path/to/conform-gt.jar"
PANEL_REF_VCF="/path/to/reference_panel_sites.vcf.gz" # match MIS panel (HRC/1KG/TOPMed)

# Exonic region filtering (GRCh38) — provide exon BED.
# To REMOVE exonic variants: we subtract using bcftools view -T ^exons.bed
EXONS_BED_GRCH38="/path/to/gencode.exons.GRCh38.bed"

# Post-imputation filtering thresholds
INFO_MIN="0.8"
MAF_MIN="0.01"

########################################
# 1) UTILITIES
########################################
LOG="${OUTDIR}/pipeline_${RUN_NAME}.log"
: > "${LOG}"

log() { echo "$(date '+%F %T')  $*" | tee -a "${LOG}"; }
die() { echo "ERROR: $*" >&2; exit 1; }

require_file() { [[ -f "$1" ]] || die "Missing file: $1"; }
require_dir()  { [[ -d "$1" ]] || die "Missing dir: $1"; }
require_exe()  { command -v "$1" >/dev/null 2>&1 || die "Missing executable in PATH: $1"; }

########################################
# 2) CHECKS
########################################
require_dir "${CEL_DIR}"
require_dir "${AXIOM_LIB_DIR}"

require_file "${AXIOM_CDF}"
[[ -x "${APT_GENO_QC}" ]] || die "APT_GENO_QC not executable: ${APT_GENO_QC}"
[[ -x "${APT_PROBESET_GENOTYPE}" ]] || die "APT_PROBESET_GENOTYPE not executable: ${APT_PROBESET_GENOTYPE}"

require_exe "${BCFTOOLS}"
require_exe "${BGZIP}"
require_exe "${TABIX}"

if [[ "${DO_LIFTOVER}" -eq 1 ]]; then
  require_file "${PICARD_JAR}"
  require_file "${CHAIN_HG19_TO_HG38}"
  require_file "${REF_FASTA_GRCH38}"
fi

require_file "${CONFORM_GT_JAR}"
require_file "${PANEL_REF_VCF}"

require_file "${EXONS_BED_GRCH38}"

log "Started: ${RUN_NAME}"
log "CEL_DIR=${CEL_DIR}"
log "OUTDIR=${OUTDIR}"
log "ARRAY_BUILD=${ARRAY_BUILD}, DO_LIFTOVER=${DO_LIFTOVER}"

########################################
# 3) STEP A — Axiom QC (DishQC / sample QC)
########################################
QC_DIR="${OUTDIR}/01_axiom_qc"
mkdir -p "${QC_DIR}"

log "STEP A: Running apt-geno-qc"
"${APT_GENO_QC}" \
  --cdf-file "${AXIOM_CDF}" \
  --cel-files "${CEL_DIR}" \
  --out-dir "${QC_DIR}"

log "QC outputs written to: ${QC_DIR}"
log "Review DishQC + call-rate metrics before proceeding."

########################################
# 4) STEP B — Genotype calling (CHP outputs)
########################################
CALL_DIR="${OUTDIR}/02_genotype_calls"
mkdir -p "${CALL_DIR}"

log "STEP B: Running apt-probeset-genotype"
# NOTE: Exact arguments vary by Axiom array bundle.
# You may need:
#   --analysis AxiomGT1  (or array-specific workflow)
#   --read-models-birdseed / --read-models-brlmmp
#   --genotype-model-file / --special-snps / etc.
#
# Use your array provider's APT best-practice command.
"${APT_PROBESET_GENOTYPE}" \
  --cdf-file "${AXIOM_CDF}" \
  --cel-files "${CEL_DIR}" \
  --analysis "AxiomGT1" \
  --out-dir "${CALL_DIR}" \
  ${AXIOM_GENO_MODEL:+--genotype-model-file "${AXIOM_GENO_MODEL}"}

log "Genotype calling complete: ${CALL_DIR}"

########################################
# 5) STEP C — Export to VCF (array-specific)
########################################
# WARNING:
# Exporting Axiom CHP to VCF is the most array-specific part.
# Some workflows export to PLINK first, then to VCF.
#
# Here we assume you will produce an initial VCF named:
#   preimpute.raw.${ARRAY_BUILD}.vcf.gz
#
# Replace this block with your validated export method.
VCF_DIR="${OUTDIR}/03_vcf"
mkdir -p "${VCF_DIR}"

RAW_VCF="${VCF_DIR}/preimpute.raw.${ARRAY_BUILD}.vcf"

log "STEP C: Export genotype calls to VCF (PLACEHOLDER)"
log "You must replace this export step with the correct Axiom export command for your array."
log "Expected output: ${RAW_VCF}"

# ---- PLACEHOLDER: STOP HERE IF YOU HAVEN'T IMPLEMENTED EXPORT ----
# die "Implement Axiom CHP -> VCF export for your array, then re-run."

# If you already created it externally, just check:
require_file "${RAW_VCF}"

log "Compress + index VCF"
"${BGZIP}" -f "${RAW_VCF}"
"${TABIX}" -f -p vcf "${RAW_VCF}.gz"

########################################
# 6) STEP D — Normalize/left-align/sort
########################################
NORM_VCF="${VCF_DIR}/preimpute.norm.${ARRAY_BUILD}.vcf.gz"

log "STEP D: Normalize + left-align + sort"
# For hg19, use hg19 fasta; for GRCh38 use GRCh38 fasta.
REF_FOR_NORM="${REF_FASTA_GRCH38}"
if [[ "${ARRAY_BUILD}" == "hg19" && -f "${REF_FASTA_HG19}" ]]; then
  REF_FOR_NORM="${REF_FASTA_HG19}"
fi

"${BCFTOOLS}" norm -m -both -f "${REF_FOR_NORM}" "${RAW_VCF}.gz" \
  | "${BCFTOOLS}" sort -Oz -o "${NORM_VCF}"

"${TABIX}" -f -p vcf "${NORM_VCF}"

########################################
# 7) STEP E — Liftover to GRCh38 (if needed)
########################################
LIFT_VCF="${VCF_DIR}/preimpute.norm.GRCh38.vcf.gz"

if [[ "${ARRAY_BUILD}" == "hg19" && "${DO_LIFTOVER}" -eq 1 ]]; then
  log "STEP E: Liftover hg19 -> GRCh38 using Picard LiftoverVcf"

  # Picard expects uncompressed VCF input; we'll stream via temp.
  TMP_IN="${VCF_DIR}/tmp.norm.hg19.vcf"
  TMP_OUT="${VCF_DIR}/tmp.lift.GRCh38.vcf"
  TMP_REJ="${VCF_DIR}/tmp.lift.rejects.vcf"

  "${BCFTOOLS}" view "${NORM_VCF}" -Ov -o "${TMP_IN}"

  java -jar "${PICARD_JAR}" LiftoverVcf \
    I="${TMP_IN}" \
    O="${TMP_OUT}" \
    REJECT="${TMP_REJ}" \
    CHAIN="${CHAIN_HG19_TO_HG38}" \
    R="${REF_FASTA_GRCH38}"

  "${BGZIP}" -f "${TMP_OUT}"
  mv "${TMP_OUT}.gz" "${LIFT_VCF}"
  "${TABIX}" -f -p vcf "${LIFT_VCF}"

  log "Liftover complete: ${LIFT_VCF}"
  log "Rejected variants: ${TMP_REJ}"

else
  # Already GRCh38
  if [[ "${ARRAY_BUILD}" == "GRCh38" ]]; then
    LIFT_VCF="${NORM_VCF}"
  else
    LIFT_VCF="${NORM_VCF}"
  fi
  log "STEP E: Liftover skipped. Using: ${LIFT_VCF}"
fi

########################################
# 8) STEP F — Conform/strand-align for Michigan
########################################
# Michigan panels have fixed REF/ALT. Conform-gt aligns your VCF to panel sites.
CONFORM_VCF="${VCF_DIR}/preimpute.conform.GRCh38.vcf.gz"

log "STEP F: Running conform-gt to match reference panel alleles/sites"
java -jar "${CONFORM_GT_JAR}" \
  gt="${LIFT_VCF}" \
  ref="${PANEL_REF_VCF}" \
  out="${VCF_DIR}/preimpute.conform.GRCh38" \
  chrom=chr1-22

require_file "${CONFORM_VCF}"
"${TABIX}" -f -p vcf "${CONFORM_VCF}"

log "Pre-imputation VCF ready for Michigan Imputation Server:"
log "  ${CONFORM_VCF}"

cat <<EOF | tee -a "${LOG}"

===================== MANUAL STEP (Michigan Imputation Server) =====================
1) Upload this file:
   ${CONFORM_VCF}

2) Choose the same panel you used for PANEL_REF_VCF.
3) Select GRCh38 / appropriate genome build.
4) Run phasing + imputation.
5) Download imputed VCFs (usually one per chromosome) into:
   ${OUTDIR}/04_imputation_download/

Then re-run this script with:
   export POST_IMPUTATION=1
====================================================================================
EOF

########################################
# 9) STEP G — Post-imputation filtering (runs only if POST_IMPUTATION=1)
########################################
if [[ "${POST_IMPUTATION:-0}" -ne 1 ]]; then
  log "POST_IMPUTATION not set; stopping after Michigan prep."
  exit 0
fi

IMPUTE_DIR="${OUTDIR}/04_imputation_download"
POST_DIR="${OUTDIR}/05_post_imputation"
mkdir -p "${POST_DIR}"

require_dir "${IMPUTE_DIR}"

log "STEP G: Post-imputation filtering"
log "Looking for imputed VCFs in: ${IMPUTE_DIR}"

# Concatenate per-chromosome VCFs if needed
# Expect: chr*.vcf.gz or similar; adjust pattern to your MIS download naming.
VCF_LIST="${POST_DIR}/imputed_vcfs.list"
ls -1 "${IMPUTE_DIR}"/*.vcf.gz > "${VCF_LIST}" || die "No *.vcf.gz found in ${IMPUTE_DIR}"

IMPUTED_MERGED="${POST_DIR}/imputed.merged.vcf.gz"
log "Merging imputed VCFs"
"${BCFTOOLS}" concat -f "${VCF_LIST}" -Oz -o "${IMPUTED_MERGED}"
"${TABIX}" -f -p vcf "${IMPUTED_MERGED}"

# Filter by INFO and MAF (field names can differ by panel)
# Common: INFO/R2 or INFO/DR2 or INFO/INFO
# We'll try R2 first; edit if your MIS outputs DR2/INFO.
FILTERED_VCF="${POST_DIR}/imputed.filtered.INFO${INFO_MIN}.MAF${MAF_MIN}.vcf.gz"

log "Filtering by INFO(R2) >= ${INFO_MIN} and MAF >= ${MAF_MIN}"
"${BCFTOOLS}" +fill-tags "${IMPUTED_MERGED}" -Oz -o "${POST_DIR}/tmp.tags.vcf.gz" -- -t MAF
"${TABIX}" -f -p vcf "${POST_DIR}/tmp.tags.vcf.gz"

"${BCFTOOLS}" view \
  -i "INFO/R2>=${INFO_MIN} && MAF>=${MAF_MIN}" \
  -Oz -o "${FILTERED_VCF}" \
  "${POST_DIR}/tmp.tags.vcf.gz"

"${TABIX}" -f -p vcf "${FILTERED_VCF}"

########################################
# 10) STEP H — Remove exonic variants (GRCh38)
########################################
NOEXON_VCF="${POST_DIR}/imputed.filtered.noexon.vcf.gz"

log "STEP H: Removing variants that fall in exonic regions (GRCh38)"
# -T ^file.bed means "exclude variants overlapping this BED"
"${BCFTOOLS}" view -T "^${EXONS_BED_GRCH38}" -Oz -o "${NOEXON_VCF}" "${FILTERED_VCF}"
"${TABIX}" -f -p vcf "${NOEXON_VCF}"

log "Done."
log "Filtered VCF:        ${FILTERED_VCF}"
log "Filtered no-exon VCF:${NOEXON_VCF}"

