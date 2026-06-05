# Suggested GitHub Issues

The GitHub connector could not create issues automatically in this session, so
use this backlog to create issues manually if needed.

## 1. Prepare anonymized sample metadata

Fill in `config/metadata_template.tsv` with real anonymized sample information.

Checklist:

- Replace example sample rows.
- Confirm `sample-id` matches the PacBio manifest.
- Remove all identifiable human information.
- Fill in group, timepoint, sample type, body site, extraction batch, PCR batch,
  and sequencing run.
- Run `python scripts/check_metadata.py config/metadata_template.tsv`.

## 2. Prepare PacBio CCS FASTQ manifest

Fill in `config/pacbio_manifest.tsv`.

Checklist:

- Use absolute FASTQ paths.
- Confirm files are demultiplexed.
- Confirm the direction column is `forward`.
- Record checksums for all raw FASTQ files.
- Keep raw FASTQ files outside GitHub.

## 3. Confirm primer sequences and DADA2 CCS parameters

Update DADA2 parameters before running denoising.

Checklist:

- Confirm exact full-length 16S primer sequences with the sequencing provider.
- Set `FRONT_PRIMER`.
- Set `ADAPTER` or reverse primer when appropriate.
- Review `MIN_LEN`, `MAX_LEN`, `MAX_EE`, and chimera settings.
- Document final parameters in `docs/SOP.md`.

## 4. Select or train full-length 16S taxonomy classifier

Choose a classifier suitable for PacBio full-length 16S reads.

Checklist:

- Choose reference database and version.
- Confirm classifier was trained for full-length 16S or the correct target
  region.
- Record source URL, training method, and checksum in `references/README.md`.
- Place the classifier at `references/classifier-full-length-16s.qza` or pass
  `CLASSIFIER=/path/to/classifier.qza`.

## 5. Run first pilot analysis

Run the starter pipeline on a small subset before the full dataset.

Checklist:

- Import data with `scripts/01_import_qiime2.sh`.
- Review demux quality visualization.
- Run `scripts/02_denoise_dada2_ccs.sh`.
- Review denoising stats.
- Run taxonomy and diversity scripts.
- Document any parameter changes.

## 6. Run exploratory TDA analysis

Run persistent homology on the exported feature table.

Checklist:

- Confirm `results/tables/feature-table.tsv` exists.
- Confirm metadata sample IDs match the feature table.
- Run `python scripts/05_tda_persistent_homology.py`.
- Review `results/figures/tda_persistence_diagram.png`.
- Review `results/figures/tda_barcode.png`.
- Check whether TDA structure is driven by group, sequencing run, extraction
  batch, PCR batch, or sample depth.
- Document interpretation as exploratory unless supported by follow-up tests.
