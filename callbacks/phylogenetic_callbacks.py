#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:08:19 2025

@author: fgraziani
"""

import base64
from dash import html, Input, Output, callback, State, dcc



import os
import sys
import io
import math

from Bio import Phylo
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)
from app import *
from neo4j_requests import *

EXPORT_DIR = "./export/phylo/"

def generate_elements(newick_str, xlen=30, ylen=30, grabbable=False):
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    def get_col_positions(tree, column_width=80):
        taxa = tree.get_terminals()

        # Some constants for the drawing calculations
        max_label_width = max(len(str(taxon)) for taxon in taxa)
        drawing_width = column_width - max_label_width - 1
    
        """Create a mapping of each clade to its column position."""
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        fudge_margin = int(math.ceil(math.log(len(taxa), 2)))
        cols_per_branch_unit = ((drawing_width - fudge_margin) /
                                float(max(depths.values())))
        return dict((clade, int(blen * cols_per_branch_unit + 1.0))
                    for clade, blen in depths.items())

    def get_row_positions(tree):
        taxa = tree.get_terminals()
        positions = dict((taxon, 2 * idx) for idx, taxon in enumerate(taxa))
    
        def calc_row(clade):
            for subclade in clade:
                if subclade not in positions:
                    calc_row(subclade)
            positions[clade] = ((positions[clade.clades[0]] +
                                 positions[clade.clades[-1]]) // 2)
    
        calc_row(tree.root)
        return positions
    
    def add_to_elements(clade, clade_id):
        children = clade.clades
    
        pos_x = col_positions[clade] * xlen
        pos_y = row_positions[clade] * ylen
    
        cy_source = {
            "data": {"id": clade_id},
            'position': {'x': pos_x, 'y': pos_y},
            'classes': 'nonterminal',
            'grabbable': grabbable
        }
        nodes.append(cy_source)
    
        if clade.is_terminal():
            cy_source['data']['name'] = clade.name
            cy_source['classes'] = 'terminal'
    
        for n, child in enumerate(children):
            support_id = clade_id + 's' + str(n)
            child_id = clade_id + 'c' + str(n)
            pos_y_child = row_positions[child] * ylen
    
            cy_support_node = {
                'data': {'id': support_id},
                'position': {'x': pos_x, 'y': pos_y_child},
                'grabbable': grabbable,
                'classes': 'support'
            }
    
            cy_support_edge = {
                'data': {
                    'source': clade_id,
                    'target': support_id,
                    'sourceCladeId': clade_id
                },
            }
    
            cy_edge = {
                'data': {
                    'source': support_id,
                    'target': child_id,
                    'length': clade.branch_length,
                    'sourceCladeId': clade_id
                },
            }
    
            if clade.confidence and clade.confidence.value:
                cy_source['data']['confidence'] = clade.confidence.value
    
            nodes.append(cy_support_node)
            edges.extend([cy_support_edge, cy_edge])
    
            add_to_elements(child, child_id)

    col_positions = get_col_positions(tree)
    row_positions = get_row_positions(tree)
    
    nodes = []
    edges = []
    
    add_to_elements(tree.clade, 'r')
    
    return nodes+edges




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

@app.callback(
    Output('cytoscape-phylo', 'elements'),
    Output('upload-status', 'children'),
    Input('upload-newick', 'contents'),
    State('upload-newick', 'filename'),
    prevent_initial_call=True,
    allow_duplicate=True
)
def update_phylo_graph(contents, filename):
    if contents is None:
        return [], "Aucun fichier chargé."

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    
    try:
        
        newick_str = decoded.decode('utf-8')
        elements = generate_elements(newick_str)
        status = f"File '{filename}' successfully load."
        return elements, status
    except Exception as e:
        return [], f"Parsing error : {str(e)}"

@app.callback(
    Output('cytoscape-phylo', 'stylesheet'),
    Input('cytoscape-phylo', 'mouseoverEdgeData')
    )
def color_children(edgeData):
    if edgeData is None:
        return stylesheet

    if 's' in edgeData['source']:
        val = edgeData['source'].split('s')[0]
    else:
        val = edgeData['source']
    
    children_style = [{
        'selector': 'edge[source *= "{}"]'.format(val),
        'style': {
            'line-color': 'blue'
        }
    }]
    
    return stylesheet + children_style


@app.callback(
    Output('cytoscape-phylo-region', 'elements'),
    Output("phylogenetic-message", "children"),
    Output("phylogenetic-page-store", "data"),
    Input('btn-plot-region', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def plot_region(n_clicks, stored_data, phylo_data):
    if not stored_data:
        return [], html.Div(html.P([
        "❌ No data to compute tree. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("gwas page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ], style=error_style)), phylo_data

    try:
        # Step 1 : compute tree of the region
        newick_str = compute_phylo_tree_from_nodes(stored_data)
        phylo_data = {"newick":newick_str}
    
        # Step2 : draw tree
        elements = generate_elements(newick_str)
        return elements, "", phylo_data
    
    except Exception as e:
        print(f"Error while computing tree : {e}")
        return [], html.Div(f"❌ Error while computing tree : {e}", style=error_style), phylo_data
    
    


@app.callback(
    Output("phylogenetic-message", "children", allow_duplicate=True),
    Input('btn-save-tree', 'n_clicks'),
    State("phylogenetic-page-store", "data"),
    prevent_initial_call=True
)
def save_tree(n_clicks, phylo_data):
    if len(phylo_data):
        save_path = os.path.join(os.getcwd(), EXPORT_DIR, "tree.nwk")
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(phylo_data["newick"])
    
        return f"File saved : {save_path}"
    else:
        return f"No data to save"