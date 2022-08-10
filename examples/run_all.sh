#!/usr/bin/env bash

CURRENT_DIR="$(dirname -- "$0")"
find "$CURRENT_DIR" -type f -name "*.py" -exec python {} \;
