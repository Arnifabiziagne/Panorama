#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 12:36:54 2025

@author: fgraziani
"""
import dash
from dash import html, Input, Output, callback, State, dcc, callback_context, ctx
from dash.exceptions import PreventUpdate
from Bio.Seq import Seq
from app import *
from neo4j_requests import *


def generate_sequence_table(sequences_dic):
    return html.Table([
    html.Thead(html.Tr([
        html.Th("Name", style={
            'border': '1px solid black',
            'padding': '6px',
            'width': '20%'
        }),
        html.Th("Sequence", style={
            'border': '1px solid black',
            'padding': '6px',
            'width': '80%',
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-word'
        })
    ])),
    html.Tbody([
        html.Tr([
            html.Td(name, style={
                'border': '1px solid black',
                'padding': '6px',
                'width': '20%'
            }),
            html.Td(seq, style={
                'border': '1px solid black',
                'padding': '6px',
                'width': '80%',
                'whiteSpace': 'pre-wrap',
                'wordBreak': 'break-word'
            })
        ]) for name, seq in sequences_dic.items()
        ])
    ], style={
    'border': '1px solid black',
    'borderCollapse': 'collapse',
    'width': '100%',
    'tableLayout': 'fixed'
    })


@app.callback(
    Output('sequences-output', 'children'),
    Input('sequences-page-store', 'data'),
    prevent_initial_call=False
)
def load_sequences_on_page_load(sequences_dic):
    if sequences_dic and isinstance(sequences_dic, dict) and len(sequences_dic) > 0:
        return generate_sequence_table(sequences_dic)
    return html.Div("No sequences to display.")


@app.callback(
    Output('sequences-page-store', 'data', allow_duplicate=True),
    Output("sequences-message", "children"),
    Input('get-sequences-btn', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    prevent_initial_call=True
)
def display_sequences(n_clicks, nodes_data):
    ctx = dash.callback_context
    if ctx.triggered_id == "get-sequences-btn" and n_clicks == 0:
        raise PreventUpdate

    if not nodes_data:
        return {}, html.Div(html.P([
        "❌ No data to compute sequences. Select a region to visualise on the ",
        dcc.Link("home page", href="/", style={'color': 'blue', 'textDecoration': 'underline'}),
        " or on the ",
        dcc.Link("Shared regions discovery page", href="/gwas", style={'color': 'blue', 'textDecoration': 'underline'})
        ], style=error_style))
    else:
        sequences = []
        genomes_nodes_dic = {}
        set_names = set()
        for n in nodes_data:
            node = nodes_data[n]
            set_names.add(node["ref_node"])
            for g in node["genomes"]:
                if g not in genomes_nodes_dic:
                    genomes_nodes_dic[g] = []
                strand = "P"
                if "strandM" in node and g in node["strandM"]:
                    strand = "M"
                genomes_nodes_dic[g].append({"start":node[g+"_position"], "node_name":node["ref_node"], "strand":strand})
        
        sorted_names_by_genome = {
            genome: {
                "names": [item["node_name"] for item in sorted(nodes, key=lambda x: x["start"])],
                "strands": [item["strand"] for item in sorted(nodes, key=lambda x: x["start"])]
            }
            for genome, nodes in genomes_nodes_dic.items()
        }
        sequences_list = get_sequence_from_names(list(set_names))
        sequences_dic = {}
        for g in sorted_names_by_genome:
            sequence = ""
            for i in range(len(sorted_names_by_genome[g]["names"])):
                if sorted_names_by_genome[g]["strands"][i] == "M":
                    sequence += Seq(sequences_list[sorted_names_by_genome[g]["names"][i]]).reverse_complement()
                else:
                    sequence += sequences_list[sorted_names_by_genome[g]["names"][i]]
            sequences_dic[g] = str(sequence)
        return sequences_dic, ""

@app.callback(
    Output('sequences-page-store', 'data', allow_duplicate=True),
    Input('url', 'pathname'),
    State('sequences-page-store', 'data'),
    prevent_initial_call=True

)
def update_sequences_on_page_load(pathname, sequences_data):
    elements_region = []
    elements_global = []
    return sequences_data
    # if sequences_data is None:
    #     sequences_data = {}
    # else:
    #     if len(sequences_data) > 0:
    #
    # if "newick_region" in pĥylo_data and pĥylo_data["newick_region"] is not None:
    #     elements_region = generate_elements(pĥylo_data["newick_region"])
    # if "newick_global" in pĥylo_data:
    #     elements_global = generate_elements(pĥylo_data["newick_global"])
    # # Display last tree button if a tree has already been computed
    # if os.path.exists(last_tree):
    #     last_tree_btn = {"display": "block"}
    # else:
    #     last_tree_btn = {"display": "none"}
    # return elements_global, elements_region, last_tree_btn