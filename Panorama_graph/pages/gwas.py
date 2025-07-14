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
        html.H5("Select genomes :"),
        dcc.Checklist(
            id='genome-list',
            options=[],
            value=[],
            labelStyle={'display': 'inline-block', 'marginRight':'10px'}
        ),
        
        html.Br(), 
        html.Div([
            html.Div([
                html.Label("Min node size to detect a shared region (integer) : ", style={"marginRight": "10px", "marginLeft":"10px"}),
                dcc.Input(id='gwas-min-node-size-int', type='number', step=1, value=10, debounce=True),
                html.Label("Min percentage of selected genomes to detect shared nodes (set to zéro for one genome min): ", style={"marginRight": "10px", "marginLeft":"10px"}),
                dcc.Input(id='gwas-min-percent_selected', type='number', step=1, value=0, debounce=True),
                html.Label("Max percentage of selected genomes to detect shared nodes : ", style={"marginRight": "10px", "marginLeft":"10px"}),
                dcc.Input(id='gwas-max-percent_selected', type='number', step=1, value=110, debounce=True),
                html.Label("Group detected nodes separate from less than this value into a same region : ", style={"marginRight": "10px", "marginLeft":"10px"}),
                dcc.Input(id='gwas-region-gap', type='number', step=1, value=10000, debounce=True),
                html.Label("Include deletion (takes more time to compute) : ", style={"marginRight": "10px", "marginLeft":"10px"}),
                dcc.Checklist(
                    id='gwas-toggle-deletion',
                    options=[{'label': 'Deletion', 'value': 'show'}],
                    value=['show'],  # Valeur cochée par défaut
                    style={'margin-bottom': '20px'}
                ),
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
                {"name": "size", "id": "size"}
            ],
            data=[],
            style_table={'overflowX': 'auto'},
            row_selectable='single'
        ),
        html.Div([
        html.Button("💾 Export to CSV", id='save-csv-button', n_clicks=0),
        html.Button("📂 Load csv", id='load-csv-button', n_clicks=0),
        html.Div(id='save-feedback'),
        dcc.Upload(
            id='upload-csv',
            children=html.Div(['Glissez un fichier CSV ici ou cliquez pour sélectionner un fichier.']),
            style={
                'display': 'none',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'padding': '10px',
            },
            multiple=False
        )
    ]),

    html.Div(id='selected-region-output')
    ], style={'padding': '20px'})



