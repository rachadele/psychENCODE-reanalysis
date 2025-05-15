#!/bin/bash
set -e

file=$1

# Check if input is provided
if [ -z "$file" ]; then
  echo "Usage: $0 <GEMMA file identifier>"
  exit 1
fi

# Run and check explicitly
if ! gemma-cli-sc getDataMatrix -f "$file"; then
  echo "Error: gemma-cli-sc failed"
  exit 1
fi
