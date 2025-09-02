#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""


import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc


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
        html.H2("Phylogenetics tree"),
        
        #Help section
        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("This page allows to display 2 phylogenetic trees :"),
                    html.Ul([
                        html.Li("Load a newick file : this allows you to load a file and display a reference tree, for example."
                                " For that, juste drag / drop or select the newick file."),   
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
            html.H4("Reference Phylogenetic tree"),
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
            html.Button("Plot tree of selected region", id="btn-plot-region"),
            dcc.Loading(
                id="loading-phylogenetic-msg",
                type="circle",
                children=html.Div(id="phylogenetic-message")
            ),
            html.Div(id='region-status', style={'margin': '10px 0'}),
            html.H4("Phylogenetic tree for selected region"),
            cyto.Cytoscape(
                id='cytoscape-phylo-region',
                elements=[],
                stylesheet=stylesheet,
                layout={'name':'preset'},
                style={'width': '100%', 'height': '1000px'}
            )
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top', 'marginLeft':'4%'})
    ],style={'padding':'20px'})


