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
        html.H2("Sequences"),
        #Help section
        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to display the sequences of the current visualized region. "
                        " It is therefore necessary to select a region to view beforehand (on the home or gwas pages)."
                        " By clicking on the Get sequences button it will get the sequences for each haplotype."),
            ])
        ], style={"marginBottom": "20px"}),
        html.Label("Click to get the sequences of selected region"),
        html.Div([
            html.Button("Get sequences", title='Before using this button, data must be displayed on home page. If data are displayed, then this will computes the sequence of each displayed haplotypes.', id='get-sequences-btn', n_clicks=0),
            dcc.Loading(
                id="loading-sequences-msg",
                #type="circle",
                children=html.Div(id="sequences-message")
            )
        ],style={'padding':'20px'}),
        html.Div(id='sequences-output', style={'marginTop': '20px'})

        ],style={'padding':'20px'})



