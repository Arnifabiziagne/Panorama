#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:48:11 2025

@author: fgraziani
"""

from gevent import monkey
monkey.patch_all()

import multiprocessing

accesslog = "logs/gunicorn/gunicorn_access.log"
errorlog = "logs/gunicorn/gunicorn_error.log"
loglevel = "error"



def on_starting(server):
    from index import start_container
    start_container()