#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-config/pacbio_manifest.tsv}"
METADATA="${METADATA:-config/metadata_template.tsv}"
OUT_DIR="${OUT_DIR:-results/qiime2}"
QZV_DIR="${QZV_DIR:-results/qzv}"
DEMUX_OUT="${DEMUX_OUT:-${OUT_DIR}/demux-pacbio-ccs.qza}"

mkdir -p "$OUT_DIR" "$QZV_DIR"

if [[ ! -f "$MANIFEST" ]]; then
  echo "Missing manifest: $MANIFEST" >&2
  exit 1
fi

if [[ ! -f "$METADATA" ]]; then
  echo "Missing metadata: $METADATA" >&2
  exit 1
fi

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path "$DEMUX_OUT" \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data "$DEMUX_OUT" \
  --o-visualization "${QZV_DIR}/demux-pacbio-ccs.qzv"

qiime metadata tabulate \
  --m-input-file "$METADATA" \
  --o-visualization "${QZV_DIR}/metadata.qzv"

echo "Created $DEMUX_OUT"
echo "Review ${QZV_DIR}/demux-pacbio-ccs.qzv before denoising."
