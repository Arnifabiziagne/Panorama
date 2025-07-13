#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""


import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc




def layout():
    return html.Div([
        dcc.Store(id="sequences-page-store", data={'sequences':[]},storage_type="session"),
        html.H1("Sequences"),
        html.Label("click to get the sequences of selected region"),
        html.Button("Get sequences", id='get-sequences-btn', n_clicks=0),
        html.Div(id='sequences-output', style={'marginTop': '20px'})

        ],style={'padding':'20px'})



