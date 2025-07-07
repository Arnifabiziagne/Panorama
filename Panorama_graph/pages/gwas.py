#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 21:55:56 2025

@author: fgraziani
"""

from dash import html, dcc, dash_table


def layout():
    return html.Div([
        dcc.Store(id="gwas-page-store", data={'checkboxes':[],'analyse':[]},storage_type="session"),
        html.H2("GWAS - Shared Region Finder"),

        # Area of genomes selection
        html.Div(id='genome-checkboxes'),
        html.H5("Sélectionnez les génomes :"),
        dcc.Checklist(
            id='genome-list',
            options=[],
            value=[],
            labelStyle={'display': 'inline-block', 'marginRight':'10px'}
        ),
        
        html.Br(), 
        html.Div([
            html.Div([
                html.Label("Min node size to detect a shared region (integer) : ", style={"marginRight": "10px"}),
                dcc.Input(id='min-node-size-int', type='number', step=1, debounce=True),
            ], style={"display": "flex", "align-items": "center", "marginRight": "20px"}
            ),
            dcc.Dropdown(id='gwas_chromosomes_dropdown', placeholder="Limit search to chromosome : ", style={
                "width": "250px",     
                "minWidth": "150px",
                "maxWidth": "100%",   
                "flexShrink": 0
            })
        ], style={"display": "flex", "flexDirection": "row", "align-items": "center", "marginBottom": "20px"}),
        
        html.Br(), 
    
        html.Button("Find shared regions", id='btn-find-shared', n_clicks=0, style={'margin': '15px 0'}),
        dcc.Loading(
            id="gwas_loading-spinner",
            type="circle",  # 'default', 'circle', or 'dot'
            children=html.Div(id="load_spinner_zone")
        ),
    
        html.Div(id='shared-status', style={'marginBottom': '15px'}),
    
        # Analyse array
        dash_table.DataTable(
            id='shared-region-table',
            columns=[
                {"name": "genome", "id": "genome"},
                {"name": "chromosome", "id": "chromosome"},
                {"name": "start", "id": "start"},
                {"name": "stop", "id": "stop"},
                {"name": "annotations", "id": "annotations"},
                {"name": "taille", "id": "taille"}
            ],
            data=[],
            style_table={'overflowX': 'auto'},
            row_selectable='single'
        ),

    html.Div(id='selected-region-output')
    ], style={'padding': '20px'})



