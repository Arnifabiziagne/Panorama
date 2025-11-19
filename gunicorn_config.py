#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:48:11 2025

@author: fgraziani
"""

from gevent import monkey
monkey.patch_all()

import multiprocessing
import os
import json

CONF_FILE = os.path.abspath("./conf.json")

if os.path.exists(CONF_FILE):
    with open(CONF_FILE, "r", encoding="utf-8") as f:
        conf_data = json.load(f)
        log_level = str(conf_data.get("gunicorn_log_level", "")).lower()
        if log_level != "" and log_level in ["info", "error", "debug", "warning", "critical"]:
            accesslog = "logs/gunicorn/gunicorn_access.log"
            errorlog = "logs/gunicorn/gunicorn_error.log"
            loglevel = "error"



def on_starting(server):
    from index import start_container
    start_container()