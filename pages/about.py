#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 13:52:04 2025

@author: fgraziani
"""


import dash_cytoscape as cyto
from dash import Dash, html,callback, dcc

from app import DB_VERSION




def layout():
    return html.Div([
        html.H2("About"),
        html.P("Panorama is a tool for managing pangenome graph based on a local graph database. It offers a development framework (back office functions) and an IHM based on this data modelisation."),
        #Help section
        html.Div([
            html.Ul([
                html.Li("Information about this tool :"),
                    html.Ul([
                        html.Li(f"Database version : {DB_VERSION}"),   
                        html.Li("Panorama is developed and maintained at INRAE in the MIA-T laboratory, Genotoul-Bioinfo team, Toulouse, France." ), 
                        html.Li([
                            "Project page: ",
                            html.A("https://forge.inrae.fr/genotoul-bioinfo/panorama",
                                   href="https://forge.inrae.fr/genotoul-bioinfo/panorama",
                                   target="_blank")
                        ]),
                    ]),
                    html.Li("Main functionnalities :"),
                        html.Ul([
                            html.Li("Generate a new database and load GFA / annotations files through the 'DB management' page."),   
                            html.Li("Search and visualize a region of the pangenome graph through the 'Home' page." ), 
                            html.Li("Get the sequences of the displayed region for each haplotype through the 'Sequences' page." ),
                            html.Li("Plot the phylogenetic of the displayed region for each haplotype through the 'Phylogenetic' page." ),
                            html.Li("Find shared regions between selected haplotypes in the pangenome graph through the 'Shared regions finder' page." ),
                            html.Li("The database can be used directly from the 'http://localhost:7474' URL and data can be manipulated with cypher langages." ),
                            html.Li("The database can be used from back office python functions." ),
                        ])
            ])
            ], style={"marginBottom": "100px"}),
        html.Hr(),
        html.Div([
            html.Img(src='/assets/images/logo_genotoul_bioinfo.jpg',
                     style={'height': '150px', 'marginRight': '100px'}),
            html.Img(src='/assets/images/logo_miat.png',
                     style={'height': '100px', 'marginRight': '100px'}),
            html.Img(src='/assets/images/logo_INRAE.png',
                     style={'height': '75px', 'marginRight': '100px'})
        ],
        style={
            'display': 'flex',
            'alignItems': 'center',
            'justifyContent': 'start',  # ou 'center', selon le besoin
            'marginBottom': '20px'
        })
        
        
        ], style={'width': '48%', 'display': 'inline-block', 'verticalAlign': 'top', 'marginLeft':'4%'})


