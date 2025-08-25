#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21/08/2025

@author: fgraziani
"""

import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_cytoscape as cyto

from neo4j_requests import *
from app import app, DB_VERSION

from neo4j_available_docker_images_conf import AVAILABLE_DOCKER_IMAGES

PREFIX_CONTAINER_NAME = "DB_"+ DB_VERSION + "_"

# Layout de la page
def layout():
    genomes = get_genomes()
    if genomes is None :
        genomes = []
    return html.Div([
        dcc.Store(id='db-management-page-store'),
        html.H2("DB Management"),
    
        html.H3(id='container-name-label'), 
        html.Hr(),
        html.H3("Create new DB"),
        html.Div([
            html.H4("Important notes:"),
                html.Ul([
                    html.Li("This operation will delete all the data."),
                    html.Li("If dump files or CSV files exist in ./data/import, they will be used to initialize data."),
                ])
            ]),
        html.Label("Select docker image to use :"),
        html.Div([
            dcc.Dropdown(
            id='docker-image-dropdown',
            options=[{'label': img, 'value': img} for img in AVAILABLE_DOCKER_IMAGES],
            value=AVAILABLE_DOCKER_IMAGES[0],
            clearable=False,
            style={'width': '400px'}
            ),
        ], style={'margin-bottom': '20px'}),
        html.Div([
           
        html.Label("📦 Neo4j container name :  "),
        html.Span(PREFIX_CONTAINER_NAME, style={
        'fontWeight': 'bold',
        'paddingRight': '5px'
    }),
            dcc.Input(
                id='container-name-input',
                type='text',
                debounce=True,
                style={'width': '400px'}
            )
        ], style={'margin-bottom': '20px'}),
        html.Button("Create new DB", id="btn-create-db", n_clicks=0),
        html.Div(id="create-db-confirmation", style={"marginTop": "10px"}),
        html.H3("Dumping DB"),
        html.Div([
            html.H4("Important notes:"),
                html.Ul([
                    html.Li("Dumping DB will generate ./data/import/neo4.dump file that can be used to load in databse creation."),
                ])
            ]),
        html.Button("Dump DB", id="btn-dump-db", n_clicks=0),
        dcc.Loading(
            type="circle",
            children=html.Div(id="create-db-message", style={"marginTop": "10px"})
        ),

        html.Hr(),
        # ---- Data Section ----

        html.H3("GFA data loading"),
        html.Br(),
        html.Div([
            html.H4("Instructions for loading a GFA files :"),
            html.Ul([
                html.Li("GFA files must be in the ./data/gfa directory."),
                html.Li("The database must be empty before loading."),
                html.Li("Use the load button to generate data directly into existant database, or the import CSV button to generate CSV data files."),
                html.Li("The CSV files generated will be used in the create new db procedure to generate data."),
                html.Li("For multi GFA processing : each GFA must refer to a chromosom and must be named as follows : either filename_chrxx.gfa, filename_xx.gfa, or xx.gfa")
            ])
        ]),
        html.Br(),
        html.Div([
            html.H4("If GFA concern only one chromosome, or if no chromosome, specify the chromosome value here (0 if no chromosome):  "),
            dcc.Input(id='db-chromosome-input', style={'width': '100px', 'marginRight': '10px'}),
        ], style={'display': 'flex', 'alignItems': 'center'}),
        html.Div([
            html.H4("Batch size. According to your ram available : bigger batch size will go faster but will consume more memory. Recommended 2 000 000):  "),
            dcc.Input(id='db-batch-size-input', type='number', value = 2000000, style={'width': '100px', 'marginRight': '10px'}),
        ], style={'display': 'flex', 'alignItems': 'center'}),
        html.Br(),
        dcc.Upload(
            id="upload-gfa",
            children=html.Div(id="upload-gfa-text", children=[
                "Drag and Drop or ",
                html.A("Select a GFA File")
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px 0'
            },
            multiple=True
        ),
    
        html.Div([
            html.Button("Load", id="btn-load-gfa", n_clicks=0),
            html.Button("Generate CSV Import file", id="btn-csv-import", n_clicks=0),
            html.Button("Create indexes", id="btn-create-index", n_clicks=0),
        ], style={'marginBottom': '10px'}),
        html.Br(),
        html.Label("Create index in database (only after csv import procedure or if this step has failed"),
        html.Div([
            html.Button("Create indexes", id="btn-create-index", n_clicks=0)
        ], style={'marginBottom': '10px'}),
        dcc.Loading(
            id="loading-gfa-msg",
            type="circle",
            children=html.Div(id="gfa-message")
        ),
    

        html.Hr(style={"margin": "30px 0"}),
    
        # ---- Annotation Upload Section ----
        html.H3("Annotations"),
        html.Br(),
        html.Div([
            html.H4("Instructions for loading annotations:"),
            html.Ul([
                html.Li("📄 Load annotation files – they must be located in the ./data/annotations directory."),
                html.Li("✅ The database must already be loaded."),
                html.Li("📌 Required indexes must be created and available in the database before loading annotations.")
            ])
        ]),
        dcc.Dropdown(
            id="dropdown-genome",
            options=[{"label": genome, "value": genome} for genome in genomes],
            placeholder="Select a reference genome"
        ),
        dcc.Upload(
            id="upload-annotation",
            children=html.Div(id="upload-annotation-text", children=[
                "Drag and Drop or ",
                html.A("Select an Annotation File")
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px 0'
            },
            multiple=False
        ),
    
        html.Div([
            html.Button("Load", id="btn-load-annotation", n_clicks=0),
            dcc.Loading(
                id="loading-annotation-msg",
                type="circle",
                children=html.Div(id="annotation-message")
            )
        ]),
        
        
        html.Hr(style={"margin": "30px 0"}),
        # --- Deleting data section ---
        html.H3("Deleting data"),
        html.Br(),
        html.Label("This operation is possible only if data are stored in the ./data directory. Databse will be stopped before."),
        html.Br(),
        # Delete data button
        html.Button("Delete", id="btn-delete", n_clicks=0),
        # Confirm deletion
        html.Div(id="delete-confirmation", style={"marginTop": "10px"}),
        dcc.Loading(
            type="circle",
            children=html.Div(id="delete-message", style={"marginTop": "10px"})
        )
    

    ])



