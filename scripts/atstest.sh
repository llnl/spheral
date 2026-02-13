#!/usr/bin/env bash

# Determine the path of this script file
SPHERAL_BIN_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
$SPHERAL_BIN_DIR/spheral $SPHERAL_BIN_DIR/../scripts/spheral_ats.py "$@"
