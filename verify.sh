#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

printf '%s\n' 'Checking release file hashes...'
sha256sum --check --strict SHA256SUMS.txt
if command -v python3 >/dev/null && python3 -E -B -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
    python_cmd=python3
elif command -v python >/dev/null && python -E -B -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
    python_cmd=python
else
    printf '%s\n' 'Python 3.11 or later is required.' >&2
    exit 1
fi

# Child interpreters also use assertions and leave the release tree clean.
unset PYTHONOPTIMIZE
export PYTHONDONTWRITEBYTECODE=1
"$python_cmd" -E -B replay.py
printf '%s\n' 'PASS: release hashes, 57 Fourier certificates, three rank bands, 40 orders, and three saved GPU logs.'
