#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc


import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)

from app import app
from neo4j_requests import *

#populates genomes checkboxes
@app.callback(
    Output('genome-list', 'options'),
    Input('shared_storage', 'data')
)
def update_genome_checkboxes(data):
    
    if not data or "genomes" not in data:
        return []
    #print("update data : " + str(data))
    return [{'label': genome, 'value': genome} for genome in data['genomes']]


#Populate chromosome droplist
@app.callback(
    Output('gwas_chromosomes_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data or "chromosomes" not in data:
        return []
    
    chromosomes = data["chromosomes"]
    options = [{"label": "All", "value": "All"}] + [
        {"label": str(chrom), "value": str(chrom)} for chrom in chromosomes
    ]
    return options


@app.callback(
    Output('shared-status', 'children'),
    Output("gwas-page-store", "data"),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('btn-find-shared', 'n_clicks'),
    State('genome-list', 'value'),
    State("gwas-page-store", "data"),
    State("min-node-size-int", 'value'),
    State("gwas_chromosomes_dropdown", 'value'),
    prevent_initial_call=True
)
def handle_shared_region_search(n_clicks, selected_genomes, data, min_node_size, chromosome):
    if min_node_size is not None and min_node_size != "" and isinstance(min_node_size, int):
        min_size = min_node_size
    else:
        min_size = 10
    if chromosome == None or chromosome == "All" :
        c = None
    else:
        c = [chromosome]
    data["checkboxes"]= selected_genomes
    if not selected_genomes:
        return "Choose at least one genome.",data, ""

    try:       
        dic_region, analyse = find_shared_regions(selected_genomes, chromosomes = c,node_min_size = min_size)
        analyse_to_plot = analyse[list(analyse.keys())[0]]
        data["analyse"] = analyse_to_plot
        return f"{len(analyse_to_plot)} shared regions found.",data, ""
    
    except Exception as e:
        return f"Erreur : {e}",data, ""
    
@app.callback(
    Output('selected-region-output', 'children', allow_duplicate=True),
    Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Output("url", "pathname"),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('shared-region-table', 'selected_rows'),
    State('shared-region-table', 'data'),
    State('shared_storage_nodes', 'data'),
    prevent_initial_call=True
)
def handle_row_selection(selected_rows, table_data, data):
    redirect = "/gwas"
    if not selected_rows:
        return no_update, data, redirect, ""

    row = table_data[selected_rows[0]]
    start = row['start']
    stop = row['stop']
    chromosome = row['chromosome']
    genome = row['genome']
    try:
        nodes = get_nodes_by_region(genome, chromosome, start, stop)
        redirect = "/"
        return html.Div([
            html.P(f"Found nodes into the region : {len(nodes)}")
        ]), nodes,redirect,""
    except Exception as e:
        return f"Erreur : {e}", data,redirect,""
    
#Restore checklist
@app.callback(
    Output('shared-region-table', 'data'),
    Output("genome-list", "value"),
    Input("gwas-page-store", "modified_timestamp"),
    Input("gwas-page-store", "data"),
    
    prevent_initial_call=True
)
def restore_checklist_state(ts, data):
    return data["analyse"], data["checkboxes"]


