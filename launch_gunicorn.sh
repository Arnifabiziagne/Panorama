#!/bin/bash
# --- Parameters ---
DASH_PORT=${1:-8050}
ENV_NAME="panorama_graph"
ENV_FILE="panorama_graph.yaml"
HASH_FILE="./.${ENV_NAME}_env_hash"

set -e

# Initialize conda
__conda_setup="$($(which conda) 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    echo "Error: conda not found."
    exit 1
fi

# Compute current YAML hash
CURRENT_HASH=$(sha1sum "$ENV_FILE" | awk '{print $1}')
PREV_HASH=""
if [ -f "$HASH_FILE" ]; then
    PREV_HASH=$(cat "$HASH_FILE")
fi

# Check if environment exists
if conda info --envs | grep -q "^$ENV_NAME\s"; then
    if [ "$CURRENT_HASH" != "$PREV_HASH" ]; then
        echo "Environment '$ENV_NAME' exists but YAML changed. Updating..."
        conda env update --name "$ENV_NAME" --file "$ENV_FILE" --prune
        echo "$CURRENT_HASH" > "$HASH_FILE"
    else
        echo "Environment '$ENV_NAME' exists and is up to date."
    fi
else
    echo "Environment '$ENV_NAME' not found. Creating..."
    conda env create --file "$ENV_FILE"
    echo "$CURRENT_HASH" > "$HASH_FILE"
fi

# Activate environment
echo "Activating environment '$ENV_NAME'"
conda activate "$ENV_NAME"
# --- Launch application ---
echo "Launching Panorama on port $DASH_PORT..."


python -m gunicorn wsgi:application -w 4 -b 0.0.0.0:$DASH_PORT --timeout 260000 --config gunicorn_config.py --preload
