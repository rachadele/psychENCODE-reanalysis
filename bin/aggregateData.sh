#!/bin/bash
# set to exit on error

study=$1
if [ -z "$study" ]; then
  echo "Usage: $0 <study>"
  exit 1
fi

gemma-cli getSingleCellDataMatrix -e $study --aggregate-by-assay --aggregate-by-preferred-cell-type-assignment 2>> aggregateData.log






