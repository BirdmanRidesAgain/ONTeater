#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TEST_READS="test_data/barcode42_ONT_salmonella.fq.gz"
PREFIX_DIRECT="smoke_direct"
PREFIX_WRAPPER="smoke_wrapper"
WORKFLOW_MODE="${WORKFLOW_MODE:-run}"
RUN_FULL="${RUN_FULL:-0}"
FLYE_ASM="${FLYE_ASM:-}"
ND_ASM="${ND_ASM:-}"
FINAL_ASM="${FINAL_ASM:-}"

echo "[smoke] ONTeater root: $ROOT_DIR"
echo "[smoke] test reads: $TEST_READS"
echo "[smoke] workflow mode: $WORKFLOW_MODE"

if [[ -n "${CONTAINER_IMAGE:-}" ]]; then
  echo "[smoke] WARNING: CONTAINER_IMAGE is ignored by smoke_test.sh (conda profile). Use scripts/smoke_test_docker.sh for container runs." >&2
fi

build_mode_args() {
  MODE_ARGS=()
  case "$WORKFLOW_MODE" in
    run|trim|assemble)
      MODE_ARGS=(--ONT_rds "$TEST_READS" --genome_size 5m --BUSCO_lineage enterobacterales)
      ;;
    postprocess)
      if [[ -z "$FLYE_ASM" || -z "$ND_ASM" ]]; then
        echo "[smoke] ERROR: WORKFLOW_MODE=postprocess requires FLYE_ASM and ND_ASM." >&2
        exit 1
      fi
      MODE_ARGS=(--flye_asm "$FLYE_ASM" --nd_asm "$ND_ASM")
      ;;
    qc)
      if [[ -z "$FINAL_ASM" ]]; then
        echo "[smoke] ERROR: WORKFLOW_MODE=qc requires FINAL_ASM." >&2
        exit 1
      fi
      MODE_ARGS=(--final_asm "$FINAL_ASM")
      ;;
    *)
      echo "[smoke] ERROR: unsupported WORKFLOW_MODE=$WORKFLOW_MODE" >&2
      exit 1
      ;;
  esac
}

build_mode_args

if [[ ! -f "$TEST_READS" ]]; then
  echo "[smoke] ERROR: test reads not found at $TEST_READS" >&2
  exit 1
fi

echo "[smoke] 1/3 help output"
nextflow run main.nf --help || true

echo "[smoke] 2/3 direct nextflow run"
nextflow run main.nf \
  --workflow "$WORKFLOW_MODE" \
  --prefix "$PREFIX_DIRECT" \
  "${MODE_ARGS[@]}"

echo "[smoke] 3/3 wrapper run"
./ONTeater \
  --workflow "$WORKFLOW_MODE" \
  --prefix "$PREFIX_WRAPPER" \
  "${MODE_ARGS[@]}"

echo "[smoke] validating expected outputs"
SAMPLE_BASE="$(basename "$TEST_READS" .fq.gz)"
RAW_VIZ="${SAMPLE_BASE}_raw_NanoPlot"
FILTERED_VIZ="${SAMPLE_BASE}_trim_filtered_NanoPlot"
TRIM_FASTQ="${SAMPLE_BASE}_trim.fq.gz"

validate_prefix() {
  local prefix="$1"
  local root="results/${prefix}"

  [[ -d "$root" ]] || { echo "[smoke] ERROR: missing $root" >&2; return 1; }
  if [[ "$WORKFLOW_MODE" == "trim" || "$WORKFLOW_MODE" == "run" ]]; then
    [[ -d "$root/reads/$RAW_VIZ" ]] || { echo "[smoke] ERROR: missing $root/reads/$RAW_VIZ (raw read NanoPlot)" >&2; return 1; }
    [[ -f "$root/reads/$RAW_VIZ/${SAMPLE_BASE}_rawNanoStats.txt" ]] || {
      echo "[smoke] ERROR: missing raw NanoStats under $root/reads/$RAW_VIZ" >&2
      return 1
    }
    [[ -d "$root/reads/$FILTERED_VIZ" ]] || {
      echo "[smoke] ERROR: missing $root/reads/$FILTERED_VIZ (filtered read NanoPlot)" >&2
      return 1
    }
    [[ -f "$root/reads/$FILTERED_VIZ/${SAMPLE_BASE}_trim_filteredNanoStats.txt" ]] || {
      echo "[smoke] ERROR: missing filtered NanoStats under $root/reads/$FILTERED_VIZ" >&2
      return 1
    }
    [[ -f "$root/reads/$TRIM_FASTQ" ]] || { echo "[smoke] ERROR: missing trimmed reads $root/reads/$TRIM_FASTQ" >&2; return 1; }
    [[ -s "$root/reads/$TRIM_FASTQ" ]] || { echo "[smoke] ERROR: trimmed reads empty: $root/reads/$TRIM_FASTQ" >&2; return 1; }
  fi
  if [[ "$WORKFLOW_MODE" == "postprocess" || "$WORKFLOW_MODE" == "run" ]]; then
    [[ -d "results/merged_assemblies" ]] || { echo "[smoke] ERROR: missing results/merged_assemblies" >&2; return 1; }
  fi
  if [[ "$WORKFLOW_MODE" == "qc" || "$WORKFLOW_MODE" == "run" ]]; then
    [[ -d "results/QC" ]] || { echo "[smoke] ERROR: missing results/QC" >&2; return 1; }
  fi
}

validate_prefix "$PREFIX_DIRECT"
validate_prefix "$PREFIX_WRAPPER"

if [[ "$RUN_FULL" == "1" ]]; then
  echo "[smoke] optional full run check (workflow=run)"
  nextflow run main.nf \
    --ONT_rds "$TEST_READS" \
    --workflow run \
    --prefix smoke_full \
    --genome_size 5m \
    --BUSCO_lineage enterobacterales
fi

echo "[smoke] PASS"
