#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 12:36:54 2025

@author: fgraziani
"""
from dash import html, Input, Output, callback, State, dcc
from Bio.Seq import Seq
from app import app
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
    Output('sequences-page-store', 'data'),
    Input('get-sequences-btn', 'n_clicks'),
    State('shared_storage_nodes', 'data'),
    prevent_initial_call=True
)
def display_sequences(n_clicks, nodes_data):
    if not nodes_data:
        return {}
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
                    sequence += sequences_list[sorted_names_by_genome[g]["names"][i]].reverse_complement()
                else:
                    sequence += sequences_list[sorted_names_by_genome[g]["names"][i]]
            sequences_dic[g] = sequence

        return sequences_dic