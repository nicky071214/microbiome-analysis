#!/usr/bin/env bash
set -euo pipefail

DEMUX="${DEMUX:-results/qiime2/demux-pacbio-ccs.qza}"
METADATA="${METADATA:-config/metadata_template.tsv}"
OUT_DIR="${OUT_DIR:-results/qiime2}"
QZV_DIR="${QZV_DIR:-results/qzv}"

FRONT_PRIMER="${FRONT_PRIMER:-REPLACE_WITH_FORWARD_PRIMER_5_TO_3}"
ADAPTER="${ADAPTER:-}"
MIN_LEN="${MIN_LEN:-1000}"
MAX_LEN="${MAX_LEN:-1800}"
TRUNC_LEN="${TRUNC_LEN:-0}"
TRIM_LEFT="${TRIM_LEFT:-0}"
MAX_EE="${MAX_EE:-2.0}"
TRUNC_Q="${TRUNC_Q:-2}"
MAX_MISMATCH="${MAX_MISMATCH:-2}"
POOLING_METHOD="${POOLING_METHOD:-independent}"
CHIMERA_METHOD="${CHIMERA_METHOD:-consensus}"
MIN_FOLD_PARENT_OVER_ABUNDANCE="${MIN_FOLD_PARENT_OVER_ABUNDANCE:-3.5}"
THREADS="${THREADS:-0}"

TABLE_OUT="${TABLE_OUT:-${OUT_DIR}/feature-table.qza}"
REP_SEQS_OUT="${REP_SEQS_OUT:-${OUT_DIR}/rep-seqs.qza}"
STATS_OUT="${STATS_OUT:-${OUT_DIR}/denoising-stats.qza}"

mkdir -p "$OUT_DIR" "$QZV_DIR"

if [[ ! -f "$DEMUX" ]]; then
  echo "Missing demultiplexed QIIME 2 artifact: $DEMUX" >&2
  exit 1
fi

if [[ "$FRONT_PRIMER" == "REPLACE_WITH_FORWARD_PRIMER_5_TO_3" ]]; then
  echo "Set FRONT_PRIMER to the exact 5-to-3 forward primer used by the sequencing provider." >&2
  echo "Example: FRONT_PRIMER='AGRGTTYGATYMTGGCTCAG' bash scripts/02_denoise_dada2_ccs.sh" >&2
  exit 1
fi

cmd=(
  qiime dada2 denoise-ccs
  --i-demultiplexed-seqs "$DEMUX"
  --p-front "$FRONT_PRIMER"
  --p-min-len "$MIN_LEN"
  --p-max-len "$MAX_LEN"
  --p-trunc-len "$TRUNC_LEN"
  --p-trim-left "$TRIM_LEFT"
  --p-max-ee "$MAX_EE"
  --p-trunc-q "$TRUNC_Q"
  --p-max-mismatch "$MAX_MISMATCH"
  --p-pooling-method "$POOLING_METHOD"
  --p-chimera-method "$CHIMERA_METHOD"
  --p-min-fold-parent-over-abundance "$MIN_FOLD_PARENT_OVER_ABUNDANCE"
  --p-n-threads "$THREADS"
  --o-table "$TABLE_OUT"
  --o-representative-sequences "$REP_SEQS_OUT"
  --o-denoising-stats "$STATS_OUT"
)

if [[ -n "$ADAPTER" ]]; then
  cmd+=(--p-adapter "$ADAPTER")
fi

"${cmd[@]}"

qiime feature-table summarize \
  --i-table "$TABLE_OUT" \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization "${QZV_DIR}/feature-table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "$REP_SEQS_OUT" \
  --o-visualization "${QZV_DIR}/rep-seqs.qzv"

qiime metadata tabulate \
  --m-input-file "$STATS_OUT" \
  --o-visualization "${QZV_DIR}/denoising-stats.qzv"

echo "Created $TABLE_OUT"
echo "Created $REP_SEQS_OUT"
echo "Review ${QZV_DIR}/feature-table.qzv and ${QZV_DIR}/denoising-stats.qzv."
