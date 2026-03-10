#!/usr/bin/env bash

# Determine the path of this script file
SPHERAL_BIN_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Adding -B to prevent byte-compiling
$SPHERAL_BIN_DIR/../.venv/bin/python -B "$@"
exit_result=$?
exit $exit_result
