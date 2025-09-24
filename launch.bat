@echo off
REM Stop the script if there is an error
SETLOCAL ENABLEEXTENSIONS
SETLOCAL ENABLEDELAYEDEXPANSION

REM Check if conda is available
where conda >nul 2>&1
IF ERRORLEVEL 1 (
    echo Error : conda not found.
    EXIT /B 1
)

REM Initialize conda (only needed if conda not already in PATH)
CALL conda activate panorama_graph
IF ERRORLEVEL 1 (
    echo Error when activating conda environment.
    EXIT /B 1
)

REM Run the Python script
python index.py
IF ERRORLEVEL 1 (
    echo Error when exceuting index.py
    EXIT /B 1
)

ENDLOCAL
