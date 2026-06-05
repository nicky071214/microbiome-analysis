# Config Files

## `metadata_template.tsv`

An anonymized sample metadata template. Replace the example rows before running
an actual analysis.

Required starter columns:

- `sample-id`: must match the manifest sample IDs exactly.
- `subject_id`: anonymized subject identifier.
- `group`: experimental or clinical group.
- `timepoint`: sample collection timepoint.
- `sample_type`: sample material, such as stool, saliva, or swab.
- `body_site`: broad body site.
- `extraction_batch`: extraction batch identifier.
- `pcr_batch`: PCR batch identifier.
- `sequencing_run`: sequencing run identifier.
- `notes`: non-sensitive notes.

## `pacbio_manifest.tsv`

QIIME 2 manifest for demultiplexed PacBio CCS/HiFi FASTQ files.

Use absolute paths because QIIME 2 import is usually run from a conda or server
environment where relative paths can become ambiguous.

For demultiplexed PacBio CCS full-length 16S reads, use `forward` as the
direction in the starter manifest.

## `config.yaml`

Project-level parameter record. The starter shell scripts expose the same core
parameters as environment variables so they can be tuned without editing files.
