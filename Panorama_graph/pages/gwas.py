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

        html.Details([
            html.Summary("ℹ️ Click here to display help"),
            html.Ul([
                html.Li("Algorithm details : The principle is to determine the regions shared between a selection of haplotypes."
                        " There are two types of detection :"),
                    html.Ul([
                        html.Li("Shared nodes : Here, the objective is to detect the nodes shared by the selected haplotypes."
                                " Several parameters are involved:"),      
                            html.Ul([
                                html.Li("Min node size : a node will be detected only if it's size is superior to this value"),
                                html.Li("Min (%) of selected haplotypes (= p): a node will be detected only if (p/100) * number of selected haplotypes are present on the node."
                                        "If set to zero it wil require at least one of the selected haplotypes."),
                                html.Li("Tolerance (%) (= t): a node with more than (t/100) * number of haplotypes present on this node will not be detected. If set to zero then nodes with a non selected haplotype will not be detected.")
                                
                            ]),
                            
                     html.Li("Deleted nodes : Here, the objective is to detect deletions shared by the selected haplotypes."
                             " This mode is activated only if 'include deletion' is checked."
                             " It detects nodes with the minimum of selected haplotypes and minimal size (see parameters) and at least one of the unselected haplotypes."
                             " If this node is following by a node with all the unselected haplotypes and only these haplotypes, then a deletion node will be detected."
                             " The following parameters are used :"),      
                         html.Ul([
                             html.Li("Min node size : a node will be detected only if it's size is superior to this value"),
                             html.Li("Min (%) of selected haplotypes (= p): a node will be detected only if (p/100) * number of selected haplotypes are present on the node."
                                     "If set to zero it wil require at least one of the selected haplotypes."),
                             
                         ])
                            
                    ]),
                 html.Ul([
                     html.Li("general settings :"),
                         html.Ul([
                             html.Li("Size of region : this size is used to group nodes separated by less than this value (in bp)."),
                             html.Li("Limit search to chromosom : If a chromosom is selected, it will look for shared region only on this chromosom."),
                             html.Li("Reference haplotype : results will be displayed only for this haplotype, including annotations. If no one is selected then the first annotated haplotype will be displayed."),
                             html.Li("Export to csv / export to csv with sequences : it will save the result into a csv file (without or with the sequences associated to each region). The file is located in the '/gwas' directory."),
                             html.Li("Load csv : it allows to load the saved csv (it must be located in the '/gwas' directory)."),
                             html.Li("First column : by clicking on the first columns it will display the region in the home page."),
                             html.Li("Size column : by clicking on the size columns it will display the sequence associated to the region."),
                             ])
                     ])
            ])
        ], style={"marginBottom": "20px"}),


        # Area of genomes selection
        html.Div(id='genome-checkboxes'),
        html.Div(html.H4("Select genomes : ", title="Select haplotypes for which you want to find shared regions."),
        style={
            'padding': '10px',
            'display': 'inline-block',
            'cursor': 'help'
        }),
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
                    ["Min (%) of selected haplotypes to detect shared nodes (set to zéro for one genome min): ",
                     html.Span("?", id="help-min-selected", n_clicks=0, style=style_help)
                ]),
                dcc.Input(id='gwas-min-percent_selected', type='number', step=1, value=80, debounce=True),
                html.Label(
                    ["Tolerance (%) for haplotypes not selected :  ",
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
                dcc.Dropdown(id='gwas_ref_genome_dropdown', placeholder="Reference haplotype : ", style={
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



