#!/bin/bash

set -e  

CONF_FILE="./conf.json"

if [ -f "$CONF_FILE" ]; then
	# Get neo4j container name
	CONTAINER_NAME=$(jq -r '.container_name // empty' "$CONF_FILE")

	# Vérifier si la valeur est non vide
	if [ -z "$CONTAINER_NAME" ]; then
	  :
	  exit 0
	else
		echo "🛑 Stopping Neo4j Docker : $CONTAINER_NAME"
		docker stop "$CONTAINER_NAME" && echo "✅ Neo4j container successfully stopped." || echo "⚠️  Fail to stop Neo4j Container."
	fi

fi



if [ -f gunicorn.pid ]; then
    PID=$(cat gunicorn.pid)
    echo "Stopping Gunicorn server (PID: $PID)..."
    kill -TERM "$PID"
    rm -f gunicorn.pid
    echo "Server stopped."
else
    echo "No file gunicorn.pid found."
fi
