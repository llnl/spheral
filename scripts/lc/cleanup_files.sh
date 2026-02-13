#!/usr/bin/env bash
set -euo pipefail

# Directory to clean or current directory if not set
TARGET_DIR=${TARGET_DIR:-$PWD}

# Safety check
if [[ ! -d "$TARGET_DIR" ]]; then
  echo "Error: '$TARGET_DIR' is not a directory" >&2
  exit 1
fi

echo "Scanning for files older than 6 months under: $TARGET_DIR"

# Find directories (not the root itself) older than 6 months
# -print0 is used to handle spaces and special characters in names
while IFS= read -r -d '' dir; do
  echo "Deleting file: $file"
  rm -rf -- "$file"
done < <(find "$TARGET_DIR" -mindepth 1 -type f -mtime +200 -print0)

echo "Cleanup complete."
