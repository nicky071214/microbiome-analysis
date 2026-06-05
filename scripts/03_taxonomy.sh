#!/usr/bin/env bash
set -euo pipefail

TABLE="${TABLE:-results/qiime2/feature-table.qza}"
REP_SEQS="${REP_SEQS:-results/qiime2/rep-seqs.qza}"
METADATA="${METADATA:-config/metadata_template.tsv}"
CLASSIFIER="${CLASSIFIER:-references/classifier-full-length-16s.qza}"
OUT_DIR="${OUT_DIR:-results/qiime2}"
QZV_DIR="${QZV_DIR:-results/qzv}"
EXPORTED_DIR="${EXPORTED_DIR:-results/exported}"
TABLES_DIR="${TABLES_DIR:-results/tables}"
THREADS="${THREADS:-1}"

TAXONOMY_OUT="${TAXONOMY_OUT:-${OUT_DIR}/taxonomy.qza}"

mkdir -p "$OUT_DIR" "$QZV_DIR" "$EXPORTED_DIR" "$TABLES_DIR"

for required_file in "$TABLE" "$REP_SEQS" "$METADATA"; do
  if [[ ! -f "$required_file" ]]; then
    echo "Missing required file: $required_file" >&2
    exit 1
  fi
done

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "Missing classifier: $CLASSIFIER" >&2
  echo "Use a classifier trained for PacBio/full-length 16S, or set CLASSIFIER=/path/to/classifier.qza." >&2
  exit 1
fi

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$REP_SEQS" \
  --p-n-jobs "$THREADS" \
  --o-classification "$TAXONOMY_OUT"

qiime metadata tabulate \
  --m-input-file "$TAXONOMY_OUT" \
  --o-visualization "${QZV_DIR}/taxonomy.qzv"

qiime taxa barplot \
  --i-table "$TABLE" \
  --i-taxonomy "$TAXONOMY_OUT" \
  --m-metadata-file "$METADATA" \
  --o-visualization "${QZV_DIR}/taxa-bar-plots.qzv"

qiime tools export \
  --input-path "$TABLE" \
  --output-path "${EXPORTED_DIR}/feature-table"

qiime tools export \
  --input-path "$TAXONOMY_OUT" \
  --output-path "${EXPORTED_DIR}/taxonomy"

qiime tools export \
  --input-path "$REP_SEQS" \
  --output-path "${EXPORTED_DIR}/rep-seqs"

if command -v biom >/dev/null 2>&1; then
  biom convert \
    -i "${EXPORTED_DIR}/feature-table/feature-table.biom" \
    -o "${TABLES_DIR}/feature-table.tsv" \
    --to-tsv
else
  echo "biom command not found; exported BIOM file remains in ${EXPORTED_DIR}/feature-table." >&2
fi

cp "${EXPORTED_DIR}/taxonomy/taxonomy.tsv" "${TABLES_DIR}/taxonomy.tsv"

echo "Created $TAXONOMY_OUT"
echo "Created ${QZV_DIR}/taxa-bar-plots.qzv"
