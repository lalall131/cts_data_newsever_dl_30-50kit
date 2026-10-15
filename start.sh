#!/usr/bin/env bash
set -Eeuo pipefail
TOOLS=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$TOOLS"
exec python3 -u "$TOOLS/tools/run_download_merge.py" "$@"
