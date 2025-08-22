#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 12:08:30 2025

@author: fgraziani
"""

from dash import html, dcc

sidebar = html.Div([
    html.H2("Menu"),
    html.Hr(),
    html.Nav([
        dcc.Link('Home', href='/', className='nav-link'),
        html.Br(),
        dcc.Link('GWAS', href='/gwas', className='nav-link'),
        html.Br(),
        dcc.Link('Phylogenetic', href='/phylogenetic', className='nav-link'),
        html.Br(),
        dcc.Link('Sequences', href='/sequences', className='nav-link'),
        html.Br(),
        dcc.Link('DB management', href='/db_management', className='nav-link'),
        html.Hr()
    ], className='navbar'),
])