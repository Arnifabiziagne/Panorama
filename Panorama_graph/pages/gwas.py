#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 21:55:56 2025

@author: fgraziani
"""

from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

style_help = {
    "cursor": "pointer",
    "color": "black",
    "fontWeight": "bold",
    "marginLeft": "5px"
}

def layout():
    return html.Div([
        dcc.Store(id="gwas-page-store",storage_type="session"),
        html.H2("Shared Region Finder"),

        # Area of genomes selection
        html.Div(id='genome-checkboxes'),
        html.H5(["Select genomes : ", html.Span("?", id="help-select_genomes", n_clicks=0, style=style_help)]),
        dcc.Checklist(
            id='genome-list',
            options=[],
            value=[],
            labelStyle={'display': 'inline-block', 'marginRight':'10px'}
        ),
        
        html.Br(), 
        html.Div([
            html.Div([
                html.Label([
                    "Min node size to detect a shared region (integer) : ",
                    html.Span("?", id="help-min-node-size", n_clicks=0, style=style_help)
                ]),
                dcc.Input(id='gwas-min-node-size-int', type='number', step=1, value=10, debounce=True),
                html.Label(
                    ["Min (%) of selected genomes to detect shared nodes (set to zéro for one genome min): ",
                     html.Span("?", id="help-min-selected", n_clicks=0, style=style_help)
                ]),
                dcc.Input(id='gwas-min-percent_selected', type='number', step=1, value=80, debounce=True),
                html.Label(
                    ["Tolerance (%) for genomes not selected :  ",
                     html.Span("?", id="help-tolerance", n_clicks=0, style=style_help)
                ]),
                dcc.Input(id='tolerance_percentage', type='number', step=1, value=10, debounce=True),
                html.Label(
                    ["Group detected nodes separate from less than this value into a same region : ",
                     html.Span("?", id="help-region-gap", n_clicks=0, style=style_help)
                ]),
                dcc.Input(id='gwas-region-gap', type='number', step=1, value=10000, debounce=True),
                html.Label(
                    ["Include deletion (takes more time to compute) : ",
                     html.Span("?", id="help-deletion", n_clicks=0, style=style_help)
                ]),
                dcc.Checklist(
                    id='gwas-toggle-deletion',
                    options=[{'label': 'Deletion', 'value': 'show'}],
                    value=['show'],  # Valeur cochée par défaut
                    style={'margin-bottom': '20px'}
                ),
            ], style={"display": "flex", "align-items": "center", "marginRight": "20px"}
            ),
                dcc.Dropdown(
                    id='gwas_chromosomes_dropdown',
                    placeholder="Limit search to chromosome : ",
                    style={
                        "width": "250px",     
                        "minWidth": "150px",
                        "maxWidth": "100%",   
                        "flexShrink": 0
                    }
                ),
        ], style={"display": "flex", "flexDirection": "row", "align-items": "center", "marginBottom": "20px"}),
        
        html.Br(), 
        html.Div([
                dcc.Dropdown(id='gwas_ref_genome_dropdown', placeholder="Reference genome : ", style={
                    "width": "250px",     
                    "minWidth": "150px",
                    "maxWidth": "100%",   
                    "flexShrink": 0
                }),
                html.Span("?", id="help-genome_ref-dropdown", n_clicks=0, style=style_help)
            ], style={"display": "flex", "alignItems": "center", "gap": "8px"}),

    
        html.Button("Find shared regions", id='btn-find-shared', n_clicks=0, style={'margin': '15px 0'}),
        dcc.Loading(
            id="gwas_loading-spinner",
            type="circle",  # 'default', 'circle', or 'dot'
            children=html.Div(id="load_spinner_zone")
        ),
    
        html.Div(id='shared-status', style={'marginBottom': '15px'}),
        html.Div(id='sequence-zone', style={"fontSize": "18px", "padding": "10px"}),
        html.Label("Tip : Click on the size value to print the region sequence"), 
        html.Div([
            html.Button("💾 Export to CSV", id='save-csv-button', n_clicks=0),
            html.Button("💾 Export to CSV with sequences", id='save-csv-with_seq-button', n_clicks=0),
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
        dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Help")),
                dbc.ModalBody(id="modal-help-text"),
                dbc.ModalFooter(
                    dbc.Button("Close", id="close-help", className="ms-auto", n_clicks=0)
                ),
            ],
            id="help-modal",
            is_open=False,
        ),
        # Analyse array
        dash_table.DataTable(
            id='shared-region-table',
            columns=[
                {"name": "genome", "id": "genome"},
                {"name": "chromosome", "id": "chromosome"},
                {"name": "start", "id": "start"},
                {"name": "stop", "id": "stop"},
                {"name": "annotation before", "id": "annotation_before"},
                {"name": "annotations", "id": "annotations"},
                {"name": "annotation after", "id": "annotation_after"},
                {"name": "size", "id": "size"}
            ],
            data=[],
            style_cell={
                'whiteSpace':'normal',
                'height':'auto',
                'textAlign':'left'},
            style_table={'overflowX': 'auto'},
            row_selectable='single',
            markdown_options={"html": True},
        ),
        

    html.Div(id='selected-region-output')
    ], style={'padding': '20px'})



