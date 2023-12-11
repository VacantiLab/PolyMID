#!/bin/bash

# retrieve the directory of the bokey script
#   it is in the same directory as this shell script
script_dir="$(dirname "$0")"
parent_dir="$(dirname "$script_dir")"
grandparent_dir="$(dirname "$parent_dir")"
script_file="${script_dir}/PlotCDF.py"

# activate the conda environment
source activate gcms
export PYTHONPATH="${PYTHONPATH}:${grandparent_dir}"

# Find a port that is not in use
port=5006
# increment the port number until we find an available one
while true; do
  if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null ; then
    port=$((port+1))
  else
    break
  fi
done

# run the bokeh server
bokeh serve --show "$script_file" --port $port
