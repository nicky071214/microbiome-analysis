#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

metadata_path = Path(sys.argv[1] if len(sys.argv) > 1 else "config/metadata_template.tsv")
required_columns = [
    "sample-id",
    "subject_id",
    "group",
    "timepoint",
    "sample_type",
    "body_site",
    "extraction_batch",
    "pcr_batch",
    "sequencing_run",
]

if not metadata_path.exists():
    raise SystemExit(f"Missing metadata file: {metadata_path}")

with metadata_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit("Metadata file has no header.")

    missing = [column for column in required_columns if column not in reader.fieldnames]
    if missing:
        raise SystemExit(f"Missing required metadata columns: {', '.join(missing)}")

    seen = set()
    duplicate_ids = set()
    row_count = 0
    for row in reader:
        row_count += 1
        sample_id = row.get("sample-id", "").strip()
        if not sample_id:
            raise SystemExit(f"Row {row_count + 1} has an empty sample-id.")
        if sample_id in seen:
            duplicate_ids.add(sample_id)
        seen.add(sample_id)

    if duplicate_ids:
        raise SystemExit(f"Duplicate sample IDs: {', '.join(sorted(duplicate_ids))}")

    if row_count == 0:
        raise SystemExit("Metadata file contains no sample rows.")

print(f"Metadata check passed: {metadata_path} ({row_count} samples)")
