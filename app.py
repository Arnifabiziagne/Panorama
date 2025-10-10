#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 11:58:49 2025

@author: fgraziani
"""

from dash import Dash
import os


#Database version. This version is associated with a database model; 
#compatibility with other versions is not guaranteed.
DB_VERSION="1.0.0"

success_style = {"color": "green", "marginTop": "10px"}
warning_style = {"color": "orange", "marginTop": "10px"}
error_style = {"color": "red", "marginTop": "10px"}

app = Dash(__name__, suppress_callback_exceptions=True)
server = app.server


    

