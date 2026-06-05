# Data Directory

This directory is for data notes and lightweight, non-sensitive files only.

Do not commit raw FASTQ, BAM, sequencing delivery archives, or identifiable
human metadata.

Suggested external data record:

```text
sample-id    file_name    storage_location    sha256    sequencing_run
sample001    sample001.fastq.gz    /secure/storage/path    REPLACE    run01
```

Recommended local-only layout:

```text
data/raw/          raw FASTQ/BAM files, ignored by git
data/external/     externally downloaded reference files, ignored by git
data/processed/    small processed files when safe to track
```

