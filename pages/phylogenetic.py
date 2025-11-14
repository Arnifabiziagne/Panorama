#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""


import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc
import logging


logger = logging.getLogger("panorama_logger")

stylesheet = [
    {
        'selector': '.nonterminal',
        'style': {
            'label': 'data(confidence)',
            'background-opacity': 0,
            "text-halign": "left",
            "text-valign": "top",
        }
    },
    {
        'selector': 'node',
        'style': {
            'label': 'data(name)',
            'font-size':'50px',
            "text-halign": "left",
            "text-valign": "center"
        }
    },
    {
        'selector': '.support',
        'style': {'background-opacity': 0}
    },
    {
        'selector': 'edge',
        'style': {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        }
    },
    {
        'selector': '.terminal',
        'style': {
            'label': 'data(name)',
            'width': 10,
            'height': 10,
            "text-valign": "center",
            "text-halign": "right",
            'background-color': '#222222'
        }
    }
]


def layout():
    return html.Div([
        dcc.Store(id="phylogenetic-page-store",storage_type="session"),
        html.H2("Phylogenetics tree"),
        
        #Help section
        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to display 2 phylogenetic trees :"),
                    html.Ul([
                        html.Li("Load a newick file : this allows you to load a file and display a reference tree, for example."
                                " For that, juste drag / drop or select the newick file."),   
                        html.Li("Plot global tree : this will compute a global tree with a either a RAxML or a distance matrix and neighbor joining (faster than raxml but less accurate)."
                                "These two methodes use a sample matrix of presence / absence of a subset of sampled nodes. Nodes traversed directly or in reverse are considered different nodes."
                                " It is possible to select a chromosome to limit the tree to this chromosome. If no chromosome selected then the tree is computed on all chromosomes."),   
                        html.Li("Plot tree of selected region : This allows you to calculate a tree for the region currently being viewed on the home page."
                                " It is therefore necessary to select a region to view beforehand (on the home or gwas pages)."
                                " The tree is constructed based on a distance matrix. This matrix is calculated using the Jaccard index, taking into account the strand and repetition of each node."
                                " The tree is then calculated using the neighbor joining algorithm."), 
                    ])
            ])
            ], style={"marginBottom": "20px"}),
        
        #First column for reference tree
        html.Div([
            dcc.Upload(
                id='upload-newick',
                children=html.Div([
                    'Drag/drop or ',
                    html.A('selct a Newick file')
                ]),
                style={
                    'width': '60%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
            html.Div([
                html.Div([
                    html.Button(
                        "Plot global tree",
                        title="This will compute the tree for the whole pangenome...",
                        style={'margin-right': '15px'},
                        id="btn-plot-global-tree"
                    ),
                    dcc.Dropdown(
                        id='method-dropdown',
                        options=[
                            {'label': 'RAxML (Maximum Likelihood)', 'value': 'raxml'},
                            {'label': 'Neighbor Joining', 'value': 'nj'}
                        ],
                        value='raxml',
                        clearable=False,
                        style={'width': '300px', 'margin-right':'20px'}
                    ),
                    html.Label(
                        "Chromosome",
                        title="The regions / annotations will be related only to this chromosome.",
                        style={'display': 'block', 'marginRight': "10px"}
                    ),
                    dcc.Dropdown(
                        id='phylogenetic_chromosomes_dropdown',
                        placeholder="Limit tree to chromosome : ",
                        style={
                            "width": "250px",
                            "minWidth": "150px",
                            "maxWidth": "100%",
                            "flexShrink": 0
                        }
                    ),
                ], style={"display": "flex", "flexDirection": "row", "alignItems": "center"}),

                html.Div([
                    html.Button(
                        "Load last computed tree",
                        title="This will load the last computed tree.",
                        id="btn-load-last-tree",
                        style={"display": "none", "marginTop": "10px"}  # un peu d'espace au-dessus
                    ),
                    dcc.Loading(
                        id="phylogenetic_loading-spinner",
                        type="circle",
                        children=html.Div(id="load_spinner_zone")
                    ),
                ], style={"marginTop": "10px"})
            ], style={"marginBottom": "20px"}),
            html.H4("Reference tree"),
            html.Div(id='upload-status'),
            cyto.Cytoscape(
                id='cytoscape-phylo',
                elements=[],
                stylesheet=stylesheet,
                layout={'name':'preset'},
                style={'width': '100%', 'height': '1000px'}
            )
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top'}),
    
        #Second column for specific region tree
        html.Div([
            html.Button("Plot tree of selected region", title="Before using this button, data must be displayed on home page. If data are displayed, then this will computes the phylogenetic tree of displayed haplotypes (see help for details).",style={'margin-right': '15px'}, id="btn-plot-region"),
            html.Button("Save tree", title="Tree will be saved into newick format.", id="btn-save-tree"),
            dcc.Download(id="download-tree"),
            dcc.Loading(
                id="loading-phylogenetic-msg",
                type="circle",
                children=html.Div(id="phylogenetic-message")
            ),
            html.Div(id='region-status', style={'margin': '10px 0'}),
            html.H4("Tree for selected region"),
            cyto.Cytoscape(
                id='cytoscape-phylo-region',
                elements=[],
                stylesheet=stylesheet,
                layout={'name':'preset'},
                style={'width': '100%', 'height': '1000px'}
            )
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top', 'marginLeft':'4%'})
    ],style={'padding':'20px'})


