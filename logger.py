#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 09:57:36 2025

@author: fgraziani
"""

import logging
import logging.handlers
import json
import os


LOGS_DIR = "./logs"
CONF_FILE = os.path.abspath("./conf.json")


def get_log_parameters():
    valid_levels = {"INFO":logging.INFO, "DEBUG":logging.DEBUG, "WARNING":logging.WARNING,
                    "ERROR":logging.ERROR, "CRITICAL":logging.CRITICAL, "NOTSET":logging.NOTSET}
    log_retention = 7
    log_mode = "both"
    log_level=logging.INFO
    log_level_str = "INFO"
    if os.path.exists(CONF_FILE):
        with open(CONF_FILE) as f:
            conf = json.load(f)
            log_retention = int(conf.get("log_retention_days",7))
            log_mode = str(conf.get("server_log_mode","both"))
            log_level_str = str(conf.get("log_level","INFO"))
            log_level = valid_levels.get(str(conf.get("log_level","INFO")).upper(), logging.INFO)

    return log_retention,log_mode,log_level,log_level_str



#Mode log : 'file', 'console' ou 'both'
def setup_logger(name="app"):
    """
    Configure main logger of PanAbyss.

    Args:
    """
    RETENTION_DAYS, SERVER_LOG_MODE, LOG_LEVEL, LOG_LEVEL_STR = get_log_parameters()
    logger = logging.getLogger(name)
    logger.setLevel(LOG_LEVEL)
    

    if not logger.handlers:
    
        formatter = logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        
        if SERVER_LOG_MODE in ("console", "both"):
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        
        if SERVER_LOG_MODE in ("file", "both"):
            os.makedirs(LOGS_DIR, exist_ok=True)
            log_file = os.path.join(LOGS_DIR, "app.log")
        
            file_handler = logging.handlers.TimedRotatingFileHandler(
                filename=log_file,
                when="midnight",       # rotation chaque jour à minuit
                interval=1,
                backupCount=RETENTION_DAYS,  # conserver X fichiers/jours
                encoding="utf-8"
            )
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        
        logger.info(f"Logger initialize in mode {SERVER_LOG_MODE} with level {LOG_LEVEL_STR}")
    return logger