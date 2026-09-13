#!/usr/bin/env bash

# Usage:
#   ./cfinspect <davs-url>

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <davs-url>" >&2
    exit 1
fi

remote_file="$1"
file_name="$(basename "${remote_file%%\?*}")"
local_file="/tmp/${USER}_${file_name}"

cleanup() {
    rm -f "$local_file"
}

trap cleanup EXIT INT TERM

echo "Downloading:"
echo "  $remote_file"
echo "to:"
echo "  $local_file"

gfal-copy -f \
    "$remote_file" \
    "file://${local_file}" || exit 1

if [[ ! -s "$local_file" ]]; then
    echo "ERROR: downloaded file is missing or empty" >&2
    exit 1
fi

cf_inspect "$local_file"
