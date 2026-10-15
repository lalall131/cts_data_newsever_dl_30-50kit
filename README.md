# CTS-Atlas download and input preparation

Download and validate public nascent-RNA sequencing data for eight Kit categories.

## SRA path recovery update — 2026-10-15

The runnable Git sources include a fix for `Expected one SRA payload` when
prefetch 3.4.1 writes an additional accession directory under `-O`.
Both direct and nested `.sra`/`.sralite` files are recognized, existing files
are validated before reuse, and fasterq-dump receives the actual accession
directory. Ambiguous candidates remain a review error.

After updating with `git pull --ff-only`, source your local tool environment
if required, then run the recovery controller in the background:

```bash
source runtime/use_local_tools.sh   # if this server uses the local curl/pigz setup
mkdir -p logs
nohup python3 -u tools/restart_after_sra_path_fix.py \
  > logs/sra_path_fix_restart.out 2>&1 < /dev/null &
tail -n 30 logs/sra_path_fix_restart.out
```

It stops only this user's dispatcher for this exact package, allows its
existing child transfers/conversions to finish, waits for the project lock,
and starts one updated run. Data and `config/site.json` are retained. Watch
for `WAIT_TRANSFER_FINISH` and then `RESTARTED`, which records the new PID/log.
`runtime/sra_path_fix_restart.json` records completion and prevents a repeated
invocation from restarting the same fixed version again.

The original 2026-10-14 ZIP below is retained unchanged as a baseline and does
not include this fix; use the updated Git sources for execution.

If Git's network fetch stalls but ordinary HTTPS file downloads work,
`updates/sra_path_fix_20261015.bundle` carries the same path-fix commit.
Download that small file using raw.githubusercontent.com or GitHub's contents
API, verify `updates/SHA256SUMS`, and run `git pull --ff-only /path/to/file.bundle main`.
This is a local Git update after the file download; it does not contact a Git
server, and it retains the same commit history and local site configuration.
The bundle requires baseline commit `796e9c8` and updates to `a7ea4f3`.

- Full manifest: **278 sample/library records, 339 Runs**.
- Added mouse batch: **76 records / 76 Runs / 120 original FASTQ.gz files, 560.856 GB**.
- Default data root: `/srv/projects/CTS_Atlas/liuql/data`.
- This repository downloads/prepares inputs. It does **not** contain eight deployed calculation pipelines or submit SPI/SR jobs.

## Reviewed singleton input: SRR19908303 (2026-10-15)

This Run is declared PAIRED in the source plan but fasterq-dump reported
33,520,612 spots, 33,520,612 zero-length reads and 33,520,612 written reads.
The first 20 inspected spots contain a zero-length TECHNICAL read followed by
one 100 nt BIOLOGICAL read. This is evidence for a read-layout review, not a
reason to redownload or to permit arbitrary orphaned paired-end inputs.

After updating Git sources, run the targeted review in the background:

```bash
mkdir -p logs
nohup python3 -u tools/review_srr19908303_layout.py \
  > logs/srr19908303_layout_review.out 2>&1 < /dev/null &
tail -n 20 logs/srr19908303_layout_review.out
```

It reads READ_LEN and READ_TYPE for **every** SRA spot and requires precisely
an empty technical read plus one nonempty biological read in every spot.
The spot count must match the existing fully validated FASTQ record count,
and the current compressed FASTQ SHA256 must still match the validation.
If any gate fails, it leaves this Run held. No SRA download or FASTQ conversion
is performed. The currently running dispatcher does not need a restart.

Successful output uses SINGLE with `declared_layout=PAIRED`, retaining original
manifests and an audit under this Run's `stats/read_layout_review.json`.
`ready.json` and the single-Run `library_ready.json` are published, with
calculation submission still disabled. Future dispatches reuse this review
only while its exact source/FASTQ fingerprints and validation evidence match.
Other paired-end inputs continue to require matching R1/R2.

The running dispatcher owns `generated/run_results.json`, so the review does
not rewrite its historical failed entry. Current repaired status is in the
Run's `state.json`, ready manifests and `stats/read_layout_repair_result.json`.
The next full dispatch can reconcile its own result summary.

`updates/single_read_review_20261015.bundle` is an HTTPS-downloadable Git
update for servers where Git network fetches stall. It requires commit
`a7ea4f3` (the installed SRA path fix); verify `updates/SHA256SUMS` and pull
the local bundle with `git pull --ff-only /path/to/file.bundle main`.

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
