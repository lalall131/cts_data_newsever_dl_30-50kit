# CTS-Atlas download and input preparation

Download and validate public nascent-RNA sequencing data for eight Kit categories.

- Full manifest: **278 sample/library records, 339 Runs**.
- Added mouse batch: **76 records / 76 Runs / 120 original FASTQ.gz files, 560.856 GB**.
- Default data root: `/srv/projects/CTS_Atlas/liuql/data`.
- This repository downloads/prepares inputs. It does **not** contain eight deployed calculation pipelines or submit SPI/SR jobs.

## Clone and download only the added three Kits

Requires Linux, Python 3.9+ and curl. No Slurm or SRA Toolkit is needed for this added ENA-FASTQ batch.

```bash
cd /srv/projects/CTS_Atlas/liuql
git clone https://github.com/lalall131/cts_data_newsever_dl_30-50kit.git
cd cts_data_newsever_dl_30-50kit
bash start.sh --check --only-added
mkdir -p logs
nohup bash start.sh --stage all --only-added --workers 2 \
  > logs/added_three_kits.out 2>&1 < /dev/null &
echo "PID=$!"
```

Run all eight categories, including the unchanged original five-Kit human batch:

```bash
bash start.sh --install-tools --check
nohup bash start.sh --stage all --workers 2 \
  > logs/all_eight_kits.out 2>&1 < /dev/null &
echo "PID=$!"
```

The project lock prevents concurrent downloaders writing the same data root. Existing verified data are reused. The original manifest includes one GSA Run requiring a separate channel.

## Added batch

| Kit | Species | GEO study | Sample records / Runs | Original FASTQ.gz GB |
|---|---|---|---:|---:|
| SPI-PBLR | Mouse | GSE144205 | 14 / 14 | 0.831 |
| scFLUENT-seq SPI | Mouse | GSE278777 | 12 / 12 | 548.058 |
| scGRO-seq SR | Mouse | GSE242176 | 50 / 50 | 11.967 |

Sample records are Source data ID/GSM sequencing records, not cell counts or necessarily independent biological replicates. PacBio records include 6 HBB-targeted samples requiring a custom reference. The scGRO subset includes 32 paired-end and 18 single-end archives, with barcode/UMI availability still to review.

Added Kits preserve exact ENA FASTQ bytes and headers, check official size/MD5, and read the full gzip/FASTQ. The resulting `raw_fastq_manifest.json` proves raw-file validation, not FLNC or single-cell preparation. Dedicated Kit adapters and reference registration remain necessary before calculation.

## Original ZIP package

The unchanged ZIP is also committed here:

- [CTS_Atlas_download_merge_8kits_20261014.zip](CTS_Atlas_download_merge_8kits_20261014.zip)
- [ZIP checksum](ZIP_SHA256.txt)

After cloning, verify it with `sha256sum -c ZIP_SHA256.txt`. The runnable files are already extracted in the repository, so unpacking the ZIP is optional. If required, unzip it into a separate new directory.

## Documentation and manifests

- [Full package instructions](PACKAGE_README_20261014.md)
- [新增三类Kit数据与下载说明](新增三类Kit数据与下载说明_20261014.md)
- `manifests/selected_runs.tsv`: all 339 Runs.
- `manifests/added_runs.tsv`: the 76 added mouse Runs, official URLs/sizes/MD5 and evidence.
- `manifests/added_libraries.tsv`: added sample records.
- `manifests/excluded_other_branches.tsv`: six same-study comparison samples not selected for the scGRO batch.
- `manifests/extra_candidates_not_selected.tsv`: 26 additional scGRO candidates.
- `evidence/`: official metadata and selection provenance.

The execution logic has Ubuntu simulation tests. Actual downloading on a new server depends on its network, disk and tool environment.
