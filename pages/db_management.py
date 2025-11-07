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
ANNOTATION_DIR = './data/annotations'
GFA_DIR = './data/gfa'

def list_annotation_files():
    """Return all .gff / .gtf / .gff3 files in the annotation directory."""
    return sorted([
        f for f in os.listdir(ANNOTATION_DIR)
        if f.lower().endswith(('.gff', '.gtf', '.gff3'))
    ])

def list_gfa_files():
    """Return all .gfa files in the annotation directory."""
    return sorted([
        f for f in os.listdir(GFA_DIR)
        if f.lower().endswith(('.gfa'))
    ])


# Layout de la page
def layout():
    genomes = get_genomes()
    if genomes is None :
        genomes = []
    return html.Div([
        dcc.Store(id='db-management-page-store'),
        html.H2("DB Management"),
        #Help section
        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to manage database. Once the data have been created, this page is normally only used to reset the database. The creation of the database depends on the data :"),
                    html.Ul([
                        html.Li("Procedure for an intermediate data volume (graph with less than 10 millions nodes): "),  
                        html.Ul([
                            html.Li("First create the database : select the neo4j docker image and give a name to the database. Then click on 'create a new DB' button."),
                            html.Li("Load the GFA : select the GFA files. If the GFA files concern a unique chromosome, enter the name of the chromosome in the appropriate field. Then click on 'Load' button."),                            
                        ]),
                        html.Li("Procedure recommended for big data (gfa with more than 10 millions nodes): "),  
                        html.Ul([
                            html.Li("First convert the GFA into csv file : select GFA files. If the GFA files concern a unique chromosome, enter the name of the chromosome in the appropriate field. Then click on 'Generate CSV import file' button."),
                            html.Li("Then create the database : select the neo4j docker image and give a name to the database. Then click on 'create a new DB' button. The csv will automatically be used to generate databse."),
                            
                        ]),
                        html.Li("Dump procedure (for intermediate data volume): "),  
                        html.Ul([
                            html.Li("If a dump file has been generated (neo4j.dump file) just check if this file is into the /data/import directory, enter the database name, select the neo4j image associated to this dump and click on 'create new db' button."),
                        ]),
                        html.Li("Load annotations files : "),  
                        html.Ul([
                            html.Li("For this step it is necessary that indexes are created. It can take a few time after the GFA loading."),
                            html.Li("For each reference genome : select the reference genome and the drag and drop or select annotations file associated to this genome. Then click on 'Load' button. Depending on the volume of data, this operation may take some time."),
                            
                        ]),
                    ])
            ])
            ], style={"marginBottom": "20px"}),
    
        html.H3(id='container-name-label'), 


        html.Hr(),
        # ---- Data Section ----

        html.H3("GFA data loading"),
        html.Br(),
        html.Div([
            html.H4("Instructions for loading a GFA files :"),
            html.Ul([
                html.Li("This step allow to prepare data before creating database. This will load a gfa file and transform it into csv files that will be imported by neo4j when creatinf the database."),
                html.Li("GFA files must be in the ./data/gfa directory."),
                html.Li("The database must be empty before loading."),
                html.Li("For multi GFA processing : each GFA must refer to a chromosom and must be named as follows : either filename_chrxx.gfa, filename_xx.gfa, or xx.gfa or a chromosome must be associated with the gfa file."),
                html.Li("To add data after database creation, use the the add data button."),
            ])
        ]),
        html.Br(),
        html.H3("GFA files loading (only for small gfa, big must be placed directly into the /data/gfa directory)"),
    
        dcc.Upload(
            id='upload-gfa-data',
            children=html.Div([
                'Dropdown or select ',
                html.A('GFA files')
            ]),
            style={
                'width': '50%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '2px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=True
        ),
        dcc.Loading(
            id="loading-gfa-upload",
            type="default", 
            children=html.Div(id='upload-gfa-output'),
        ),
    
        html.Br(),
        # html.Div([
        #     html.H4("If GFA concern only one chromosome, or if no chromosome, specify the chromosome value here (0 if no chromosome):  ", style={'margin-right': '20px'}),
        #     dcc.Input(id='db-chromosome-input', style={'width': '100px', 'marginRight': '10px'}),
        # ], style={'display': 'flex', 'alignItems': 'center'}),
        html.Div([
            html.H4("Batch size. According to your ram available : bigger batch size will go faster but will consume more memory. Recommended 2 000 000):  ", style={'margin-right': '20px'}),
            dcc.Input(id='db-batch-size-input', type='number', value = 2000000, style={'width': '100px', 'marginRight': '10px'}),
        ], style={'display': 'flex', 'alignItems': 'center'}),
        html.Br(),
        #Checklist of annotations files

        html.Div([
            html.Div("GFA file", style={'fontWeight': 'bold', 'wordBreak': 'break-word'}),
            html.Div("Chromosome (if applicable)", title="If GFA concern only one chromosome, or if no chromosome, specify the chromosome value here (0 if no chromosome) or let this value unset if multiple chromosomes.", style={'fontWeight': 'bold'})
        ],
            style={
                'display': 'grid',
                'gridTemplateColumns': '350px 230px',
                'alignItems': 'center',
                'columnGap': '10px',
                'marginBottom': '8px'
            }),

        html.Div(
            id='gfa-files-container',
            children=[
                html.Div([
                    html.Div(
                        dcc.Checklist(
                            id={'type': 'gfa-checkbox', 'index': f},
                            options=[{'label': f, 'value': f}],
                            value=[],
                            labelStyle={
                                'display': 'inline-block',
                                'whiteSpace': 'normal',
                                'wordBreak': 'break-word',
                                'width': '100%'
                            }
                        )
                    ),

                    html.Div(
                        dcc.Input(
                            id={'type': 'gfa-input', 'index': f},
                            type='text',
                            placeholder='Enter chromosome (optional)',

                            style={'width': '100%'}
                        )
                    )
                ],
                    style={
                        'display': 'grid',
                        'gridTemplateColumns': '350px 230px',
                        'alignItems': 'center',
                        'columnGap': '10px',
                        'marginBottom': '6px'
                    })
                for f in list_gfa_files()
            ]
        ),

        html.Br(),

    
        html.Div([

            html.Button("Generate CSV Import file", title="Recommended procedure for big data : This will generates data into csv files. This files will then be used when creating a new databse.", id="btn-csv-import", n_clicks=0,style={'marginRight': '10px'}),
            html.Button("Add data (DB must be created)",
                        title="Only if the database has been created.",
                        id="btn-load-gfa", n_clicks=0),
        ], style={'marginBottom': '10px'}),
        dcc.Loading(
            id="loading-gfa-msg",
            type="circle",
            children=html.Div(id="gfa-message")
        ),

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

            html.Label("📦 Neo4j container name :  ", title='Set a name for your docker container.', ),
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
        html.Button("Create new DB",
                    title="This will create a new database. If data already exists they will be deleted.",
                    id="btn-create-db", n_clicks=0),
        html.Div(id="create-db-confirmation", style={"marginTop": "10px"}),

        html.H3("Operation on DB"),
        html.Label("Create index or stats in database (if these steps have failed).",
                   style={'display': 'block', 'marginBottom': '8px'}),
        html.Div([
            html.Button("Create indexes",
                        title="This will generate / regenerate indexes. If the already exists there will be no action.",
                        id="btn-create-index", n_clicks=0,
                        style={'marginRight': '10px'}
    ),
            html.Button("Create stats", title="Recuperation procedure for stats node if not exist after loading data.",
                        id="btn-create-stats", n_clicks=0)
        ], style={'marginBottom': '10px'}),

        html.H3("Dumping DB"),
        html.Div([
            html.H4("Important notes:"),
            html.Ul([
                html.Li(
                    "Dumping DB will generate ./data/import/neo4.dump file that can be used to load in databse creation."),
            ])
        ]),
        html.Br(),

        html.Button("Dump DB", id="btn-dump-db", n_clicks=0),
        dcc.Loading(
            type="circle",
            children=html.Div(id="create-db-message", style={"marginTop": "10px"})
        ),
    

        html.Hr(style={"margin": "30px 0"}),
    
        # ---- Annotation Upload Section ----
        html.H3("Annotations"),
        html.Br(),
        html.Div([
            html.H4("Instructions for loading annotations:"),
            html.Ul([
                html.Li("Load annotation files – files must be located in the ./data/annotations directory."),
                html.Li("The database must already be loaded."),
                html.Li("Required indexes must be created and available in the database before loading annotations."),          ])
        ]),
        html.Br(),
        html.H3("Annotations files loading (it is recommended to put big annotations files directly into the /data/annotations directory)"),
    
        dcc.Upload(
            id='upload-annotations-data',
            children=html.Div([
                'Dropdown or select ',
                html.A('GFF or GTF files')
            ]),
            style={
                'width': '50%',
                'height': '100px',
                'lineHeight': '100px',
                'borderWidth': '2px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=True  # Permet l’upload de plusieurs fichiers
        ),
        dcc.Loading(
            id="loading-annotations-upload",
            type="default",
            children=html.Div(id='upload-annotations-output'),
        ),
    
        html.Br(),
        html.Div([
            html.Div("Annotation file", style={'fontWeight': 'bold', 'wordBreak': 'break-word'}),
            html.Div("Associated genome", style={'fontWeight': 'bold'})
        ],
            style={
                'display': 'grid',
                'gridTemplateColumns': '350px 230px',  # largeur fixe pour les colonnes
                'alignItems': 'center',
                'columnGap': '10px',
                'marginBottom': '8px'
            }),


        html.Div(
            id='annotations-files-container',
            children=[
                html.Div([

                    html.Div(
                        dcc.Checklist(
                            id={'type': 'annotation-checkbox', 'index': f},
                            options=[{'label': f, 'value': f}],
                            value=[],
                            labelStyle={
                                'display': 'inline-block',
                                'whiteSpace': 'normal',
                                'wordBreak': 'break-word',
                                'width': '100%'
                            }
                        )
                    ),


                    html.Div(
                        dcc.Dropdown(
                            id={'type': 'annotation-dropdown', 'index': f},
                            options=[{"label": genome, "value": genome} for genome in genomes],
                            placeholder="Select genome",
                            style={'width': '100%'}
                        )
                    )
                ],
                    style={
                        'display': 'grid',
                        'gridTemplateColumns': '350px 230px',
                        'alignItems': 'center',
                        'columnGap': '10px',
                        'marginBottom': '6px'
                    })
                for f in list_annotation_files()
            ]
        ),

        html.Br(),
    
        html.Div([
            html.Button("Load", title="This will load annotations into database. Indexes must be created before (can take some time for big data).", id="btn-load-annotations-with-link", n_clicks=0),
            #html.Button("Load only annotations", id="btn-load-only-annotations", n_clicks=0),
            #html.Button("Link annotations", id="btn-link-annotations", n_clicks=0),
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
        html.Label("This operation is possible only if data are stored in the ./data directory. Database will be stopped before."),
        html.Br(),
        # Delete data button
        html.Button("Delete", title="This will delete all data and indexes in the database.", id="btn-delete", n_clicks=0),
        # Confirm deletion
        html.Div(id="delete-confirmation", style={"marginTop": "10px"}),
        dcc.Loading(
            type="circle",
            children=html.Div(id="delete-message", style={"marginTop": "10px"})
        )
    

    ])



