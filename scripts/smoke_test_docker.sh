#!/usr/bin/env bash
# Smoke test ONTeater using the Docker profile (single container image, no per-process conda).
# Prereqs: Docker running locally; Nextflow on PATH.
#
# Usage:
#   ./scripts/smoke_test_docker.sh
#
# Optional:
#   CONTAINER_IMAGE=my/oneteater:1.0.0-amd64 ./scripts/smoke_test_docker.sh
#   SKIP_IMAGE_BUILD=1 ./scripts/smoke_test_docker.sh   # assume image already exists
#   WORKFLOW_MODE=trim ./scripts/smoke_test_docker.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TEST_READS="test_data/barcode42_ONT_salmonella.fq.gz"
PREFIX_DIRECT="smoke_direct_docker"
PREFIX_WRAPPER="smoke_wrapper_docker"
WORKFLOW_MODE="${WORKFLOW_MODE:-run}"
RUN_FULL="${RUN_FULL:-0}"
FLYE_ASM="${FLYE_ASM:-}"
ND_ASM="${ND_ASM:-}"
FINAL_ASM="${FINAL_ASM:-}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-oneteater:1.0.0-amd64}"
SKIP_IMAGE_BUILD="${SKIP_IMAGE_BUILD:-0}"
STUB_RUN="${STUB_RUN:-0}"

echo "[smoke-docker] ONTeater root: $ROOT_DIR"
echo "[smoke-docker] container image: $CONTAINER_IMAGE"
echo "[smoke-docker] test reads: $TEST_READS"
echo "[smoke-docker] workflow mode: $WORKFLOW_MODE"

if ! command -v docker >/dev/null 2>&1; then
  echo "[smoke-docker] ERROR: docker not found on PATH." >&2
  exit 1
fi

if [[ "$SKIP_IMAGE_BUILD" != "1" ]]; then
  if ! docker image inspect "$CONTAINER_IMAGE" >/dev/null 2>&1; then
    echo "[smoke-docker] image not found locally; building: $CONTAINER_IMAGE"
    docker buildx build --platform linux/amd64 -t "$CONTAINER_IMAGE" -f docker/Dockerfile --load .
  else
    echo "[smoke-docker] using existing image: $CONTAINER_IMAGE"
  fi
fi

build_mode_args() {
  MODE_ARGS=()
  case "$WORKFLOW_MODE" in
    run|trim|assemble)
      MODE_ARGS=(--ONT_rds "$TEST_READS" --genome_size 5m --BUSCO_lineage enterobacterales)
      ;;
    postprocess)
      if [[ -z "$FLYE_ASM" || -z "$ND_ASM" ]]; then
        echo "[smoke-docker] ERROR: WORKFLOW_MODE=postprocess requires FLYE_ASM and ND_ASM." >&2
        exit 1
      fi
      MODE_ARGS=(--flye_asm "$FLYE_ASM" --nd_asm "$ND_ASM")
      ;;
    qc)
      if [[ -z "$FINAL_ASM" ]]; then
        echo "[smoke-docker] ERROR: WORKFLOW_MODE=qc requires FINAL_ASM." >&2
        exit 1
      fi
      MODE_ARGS=(--final_asm "$FINAL_ASM")
      ;;
    *)
      echo "[smoke-docker] ERROR: unsupported WORKFLOW_MODE=$WORKFLOW_MODE" >&2
      exit 1
      ;;
  esac
}

build_mode_args

if [[ ! -f "$TEST_READS" ]]; then
  echo "[smoke-docker] ERROR: test reads not found at $TEST_READS" >&2
  exit 1
fi

NF_DOCKER=( -profile docker --container_image "$CONTAINER_IMAGE" )
if [[ "$STUB_RUN" == "1" ]]; then
  NF_DOCKER+=( -stub-run )
fi

echo "[smoke-docker] 1/3 help output"
nextflow run main.nf "${NF_DOCKER[@]}" --help || true

echo "[smoke-docker] 2/3 direct nextflow run"
nextflow run main.nf "${NF_DOCKER[@]}" \
  --workflow "$WORKFLOW_MODE" \
  --prefix "$PREFIX_DIRECT" \
  "${MODE_ARGS[@]}"

echo "[smoke-docker] 3/3 wrapper run"
if [[ "$STUB_RUN" == "1" ]]; then
  WRAPPER_STUB=( --stub )
else
  WRAPPER_STUB=()
fi
./ONTeater \
  -pro docker \
  --container_image "$CONTAINER_IMAGE" \
  "${WRAPPER_STUB[@]}" \
  --workflow "$WORKFLOW_MODE" \
  --prefix "$PREFIX_WRAPPER" \
  "${MODE_ARGS[@]}"

echo "[smoke-docker] validating expected outputs"
SAMPLE_BASE="$(basename "$TEST_READS" .fq.gz)"
RAW_VIZ="${SAMPLE_BASE}_raw_NanoPlot"
FILTERED_VIZ="${SAMPLE_BASE}_trim_filtered_NanoPlot"
TRIM_FASTQ="${SAMPLE_BASE}_trim.fq.gz"

validate_prefix() {
  local prefix="$1"
  local root="results/${prefix}"

  [[ -d "$root" ]] || { echo "[smoke-docker] ERROR: missing $root" >&2; return 1; }
  if [[ "$WORKFLOW_MODE" == "trim" || "$WORKFLOW_MODE" == "run" ]]; then
    [[ -d "$root/reads/$RAW_VIZ" ]] || { echo "[smoke-docker] ERROR: missing $root/reads/$RAW_VIZ (raw read NanoPlot)" >&2; return 1; }
    [[ -f "$root/reads/$RAW_VIZ/${SAMPLE_BASE}_rawNanoStats.txt" ]] || {
      echo "[smoke-docker] ERROR: missing raw NanoStats under $root/reads/$RAW_VIZ" >&2
      return 1
    }
    [[ -d "$root/reads/$FILTERED_VIZ" ]] || {
      echo "[smoke-docker] ERROR: missing $root/reads/$FILTERED_VIZ (filtered read NanoPlot)" >&2
      return 1
    }
    [[ -f "$root/reads/$FILTERED_VIZ/${SAMPLE_BASE}_trim_filteredNanoStats.txt" ]] || {
      echo "[smoke-docker] ERROR: missing filtered NanoStats under $root/reads/$FILTERED_VIZ" >&2
      return 1
    }
    [[ -f "$root/reads/$TRIM_FASTQ" ]] || { echo "[smoke-docker] ERROR: missing trimmed reads $root/reads/$TRIM_FASTQ" >&2; return 1; }
    [[ -s "$root/reads/$TRIM_FASTQ" ]] || { echo "[smoke-docker] ERROR: trimmed reads empty: $root/reads/$TRIM_FASTQ" >&2; return 1; }
  fi
  if [[ "$WORKFLOW_MODE" == "postprocess" || "$WORKFLOW_MODE" == "run" ]]; then
    local asm_dir="$root/assemblies"
    compgen -G "$asm_dir/*_major_merged.fa" >/dev/null || { echo "[smoke-docker] ERROR: missing merged assembly under $asm_dir" >&2; return 1; }
    compgen -G "$asm_dir/*_major_merged_purged.fa" >/dev/null || { echo "[smoke-docker] ERROR: missing purged assembly under $asm_dir" >&2; return 1; }
  fi
  if [[ "$WORKFLOW_MODE" == "qc" || "$WORKFLOW_MODE" == "run" ]]; then
    local asm_dir="$root/assemblies"
    compgen -G "$asm_dir/*_final_quast_report.txt" >/dev/null || { echo "[smoke-docker] ERROR: missing QUAST report under $asm_dir" >&2; return 1; }
    compgen -G "$asm_dir/*_final_compleasm_report.txt" >/dev/null || { echo "[smoke-docker] ERROR: missing Compleasm report under $asm_dir" >&2; return 1; }
  fi
}

validate_prefix "$PREFIX_DIRECT"
validate_prefix "$PREFIX_WRAPPER"

if [[ "$RUN_FULL" == "1" ]]; then
  echo "[smoke-docker] optional duplicate full run (workflow=run)"
  nextflow run main.nf "${NF_DOCKER[@]}" \
    --ONT_rds "$TEST_READS" \
    --workflow run \
    --prefix smoke_full_docker \
    --genome_size 5m \
    --BUSCO_lineage enterobacterales
fi

echo "[smoke-docker] PASS"
