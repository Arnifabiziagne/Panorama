#!/bin/bash

if [ -f gunicorn.pid ]; then
    PID=$(cat gunicorn.pid)
    echo "Stopping Gunicorn server (PID: $PID)..."
    kill -TERM "$PID"
    rm -f gunicorn.pid
    echo "Server stopped."
else
    echo "No file gunicorn.pid found."
fi
