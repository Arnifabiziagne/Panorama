#!/bin/bash
DASH_PORT=${1:-8050}
# Stop on error
set -e

__conda_setup="$($(which conda) 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    echo "Erreur : conda not found."
    exit 1
fi

conda activate panorama_graph

echo "Launching Panorama on port $DASH_PORT..."

python index.py --port $DASH_PORT
