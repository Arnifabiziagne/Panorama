@echo off
SETLOCAL ENABLEEXTENSIONS
SETLOCAL ENABLEDELAYEDEXPANSION

:: --- Paramètres ---
SET ENV_NAME=panorama_graph
SET ENV_FILE=panorama_graph.yaml
SET HASH_FILE=.\.%ENV_NAME%_env_hash
SET DASH_PORT=8050
IF NOT "%1"=="" SET DASH_PORT=%1

:: --- Vérifie que conda est disponible ---
where conda >nul 2>&1
IF ERRORLEVEL 1 (
    echo Error: conda not found.
    EXIT /B 1
)

:: --- Calculer le hash actuel du YAML ---
for /f "tokens=*" %%H in ('certutil -hashfile "%ENV_FILE%" SHA1 ^| findstr /v "hash" ^| findstr /r /v "^$"') do set CURRENT_HASH=%%H

:: --- Lire le hash précédent s'il existe ---
SET PREV_HASH=
IF EXIST "%HASH_FILE%" (
    set /p PREV_HASH=<"%HASH_FILE%"
)

:: --- Vérifie si l'environnement existe ---
conda env list | findstr /R /C:"^%ENV_NAME%\s" >nul
IF ERRORLEVEL 0 (
    REM L'environnement existe
    IF NOT "!CURRENT_HASH!"=="!PREV_HASH!" (
        echo Environment %ENV_NAME% exists but YAML changed. Updating...
        conda env update --name %ENV_NAME% --file %ENV_FILE% --prune
        echo !CURRENT_HASH! > "%HASH_FILE%"
    ) ELSE (
        echo Environment %ENV_NAME% exists and is up to date.
    )
) ELSE (
    echo Environment %ENV_NAME% not found. Creating...
    conda env create --file %ENV_FILE%
    echo !CURRENT_HASH! > "%HASH_FILE%"
)

:: --- Active l'environnement ---
CALL conda activate %ENV_NAME%
IF ERRORLEVEL 1 (
    echo Error when activating conda environment.
    EXIT /B 1
)

:: --- Lancer l'application Dash ---
echo Launching Panorama on port %DASH_PORT%...
python -m gunicorn wsgi:application -w 4 -b 0.0.0.0:%DASH_PORT% --timeout 260000 --config gunicorn_config.py --preload

IF ERRORLEVEL 1 (
    echo Error when executing index.py
    EXIT /B 1
)

ENDLOCAL
