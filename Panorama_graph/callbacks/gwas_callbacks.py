#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc, ctx, callback_context


import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)

from app import app
from neo4j_requests import *
import base64
import io

EXPORT_DIR = "./export/gwas/"

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

#Populate genome droplist
@app.callback(
    Output('gwas_ref_genome_dropdown', 'options'),
    Input('shared_storage', 'data')
)
def update_dropdown(data):
    if not data or "genomes" not in data:
        return []
    
    genomes = data["genomes"]
    options =  [{"label": str(g), "value": str(g)} for g in genomes]
    return options


@app.callback(
    Output('shared-status', 'children'),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('btn-find-shared', 'n_clicks'),
    State('genome-list', 'value'),
    State("gwas-page-store", "data"),
    State("gwas-min-node-size-int", 'value'),
    State("gwas-min-percent_selected", 'value'),
    State("tolerance_percentage", 'value'),
    State("gwas-region-gap", 'value'),
    State('gwas-toggle-deletion', 'value'),
    State("gwas_chromosomes_dropdown", 'value'),
    State("gwas_ref_genome_dropdown", 'value'),
    prevent_initial_call=True
)
def handle_shared_region_search(n_clicks, selected_genomes, data, min_node_size, min_percent_selected, tolerance_percentage, region_gap, deletion_checkbox, chromosome, ref_genome):
    if min_node_size is not None and min_node_size != "" and isinstance(min_node_size, int):
        min_size = min_node_size
    else:
        min_size = 10
    if chromosome == None or chromosome == "All" :
        c = None
    else:
        c = [chromosome]
    if data is None:
        data = {}
    data["checkboxes"]= selected_genomes
    if not selected_genomes:
        return "Choose at least one genome.",data, ""
    deletion = False
    if 'show' in deletion_checkbox : 
        deletion = True
    
    try:       
        dic_region, analyse = find_shared_regions(selected_genomes, genome_ref = ref_genome, chromosomes = c,node_min_size = min_size, nodes_max_gap=region_gap, deletion = deletion, min_percent_selected_genomes=min_percent_selected, tolerance_percentage = tolerance_percentage)
           
        #take an annotated genome if no reference genome selected
        if ref_genome is None or ref_genome == "":
            no_annotations = True
            i = 0
            keys = list(analyse.keys())
            while no_annotations and i < len(keys):
                current_key = keys[i]
                r = 0
                while r < len(analyse[current_key]) and no_annotations :
                    if len(analyse[current_key][r]["annotations"]) > 0:
                        no_annotations = False
                        analyse_to_plot = analyse[current_key]
                    r += 1
                i += 1
        else:
            analyse_to_plot = analyse[ref_genome]
                
        #print("analyse to plot : " + str(analyse_to_plot))
        
        for r in range(len(analyse_to_plot)):
            annotation = ""
            if len(analyse_to_plot[r]["annotations"]) > 0:
                for annot in analyse_to_plot[r]["annotations"]:
                    annotation += annot["gene_name"] + "\n"
            analyse_to_plot[r]["annotations"] = annotation
            if "annotation_before" in analyse_to_plot[r]:
                annot_before = "gene_name : " + analyse_to_plot[r]["annotation_before"]["gene_name"] \
                                +"\nDistance : " + str(analyse_to_plot[r]["annotation_before"]["distance"])
                analyse_to_plot[r]["annotation_before"] = annot_before
            
            if "annotation_after" in analyse_to_plot[r]:
                annot_after = "gene_name : " + analyse_to_plot[r]["annotation_after"]["gene_name"] \
                                +"\nDistance : " + str(analyse_to_plot[r]["annotation_after"]["distance"])
                analyse_to_plot[r]["annotation_after"] = annot_after

        data["analyse"] = analyse_to_plot
        return f"{len(analyse_to_plot)} shared regions found.",data, ""
    
    except Exception as e:
        return f"Error : {e}",data, ""
    
@app.callback(
    Output('selected-region-output', 'children', allow_duplicate=True),
    Output('shared_storage_nodes', 'data',allow_duplicate=True),
    Output("url", "pathname"),
    Output("home-page-store", "data", allow_duplicate=True),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('shared-region-table', 'selected_rows'),
    State('shared-region-table', 'data'),
    State('shared_storage_nodes', 'data'),
    State('home-page-store', 'data'),
    prevent_initial_call=True
)
def handle_row_selection(selected_rows, table_data, data, home_page_data):
    redirect = "/gwas"
    if home_page_data is None:
        home_page_data = {}
    if not selected_rows:
        return no_update, data, redirect, home_page_data, ""
    print(table_data[selected_rows[0]])
    row = table_data[selected_rows[0]]
    print("selected row to plot : " +str(row))
    start = int(row['start'])
    stop = int(row['stop'])
    chromosome = row['chromosome']
    genome = row['genome']
    print("search region genome " +str(genome) + " chromosome " + str(chromosome) + " start " + str(start) + " stop " + str(stop))
    try:
        nodes = get_nodes_by_region(genome, str(chromosome), start, stop)
        home_page_data["selected_genome"]=genome
        home_page_data["selected_chromosome"]=chromosome
        home_page_data["start"]=start
        home_page_data["end"]=stop
        redirect = "/"
        return html.Div([
            html.P(f"Found nodes into the region : {len(nodes)}")
        ]), nodes,redirect,home_page_data,""
    except Exception as e:
        return f"Erreur : {e}", data,redirect,home_page_data,""
    
#Restore checklist
@app.callback(
    Output('shared-status', 'children',allow_duplicate=True),
    Output('shared-region-table', 'data',allow_duplicate=True),
    Output("genome-list", "value"),
    Input("gwas-page-store", "modified_timestamp"),
    Input("gwas-page-store", "data"),
    State('shared-region-table', 'data'),
    
    prevent_initial_call=True
)
def restore_checklist_state(ts, data, table_data):
    analyse = table_data
    checkbox = []
    if data is not None: 
        if "analyse" in data:
            analyse = data["analyse"]            
        if "checkboxes" in data:
            checkbox = data["checkboxes"]
    return f"{len(analyse)} shared regions found.",analyse, checkbox

#Callback to save the gwas data table into csv file
@app.callback(
    Output('save-feedback', 'children'),
    Output("load_spinner_zone", "children", allow_duplicate=True),
    Input('save-csv-button', 'n_clicks'),
    Input('save-csv-with_seq-button', 'n_clicks'),
    State('shared-region-table', 'data'),
    prevent_initial_call=True
)
def save_csv(n_clicks, n_clicks_seq, table_data):
    #print(f"Callback triggered: n_clicks={n_clicks}, table_data={table_data}")
    triggered_id = ctx.triggered_id
    export_sequences = False
    if triggered_id == 'save-csv-with_seq-button':
        export_sequences = True
    if not table_data:
        return "No data.",""
    df = pd.DataFrame(table_data)
    if export_sequences :
        sequences = []
        for row in tqdm(table_data):
            sequences.append(get_sequence_from_position(row['genome'], row['chromosome'], row['start'], row['stop']))
        df["sequence"] = sequences
    save_path = os.path.join(os.getcwd(), EXPORT_DIR, "shared_regions.csv")
    print("save path : " + str(save_path))
    df.to_csv(save_path, index=False)
    
    return f"File saved : {save_path}",""

#Callback to loads csv file into data table
@app.callback(
    Output('shared-status', 'children',allow_duplicate=True),
    Output('shared-region-table', 'data',allow_duplicate=True),
    Output("gwas-page-store", "data", allow_duplicate=True),
    Input('upload-csv', 'contents'),
    State('upload-csv', 'filename'),
    State("gwas-page-store", "data"),
    prevent_initial_call=True
)
def load_csv(contents, filename, gwas_page_store):
    print("load csv file")
    if contents is None:
        return None, None, gwas_page_store
    if gwas_page_store is None:
        gwas_page_store = {}
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        analyse = df[[c for c in df.columns if c!= "sequence"]].to_dict('records')
        gwas_page_store["analyse"] = analyse
        print("csv file loaded")
        return f"{len(analyse)} shared regions found.", analyse, gwas_page_store
    except Exception as e:
        print(f"Error while loading file : {e}")
        return None, None, gwas_page_store
    

@app.callback(
    Output('upload-csv', 'style'),
    Input('load-csv-button', 'n_clicks'),
    prevent_initial_call=True
)
def show_upload_area(n_clicks):
    return {
        'display': 'block',
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'padding': '10px',
        'marginTop': '10px'
    }


@app.callback(
    Output('sequence-zone', 'children'),
    Input("shared-region-table", "active_cell"),
    State('shared-region-table', 'data'),
    prevent_initial_call=True,
)
def display_sequence_on_button_click(active_cell, table_data):
    if active_cell and active_cell['column_id'] == 'size':
        row_index = active_cell["row"]
        row = table_data[row_index]
        sequence = get_sequence_from_position(row['genome'], row['chromosome'], row['start'], row['stop'])
        return html.Div([
            html.P(f"Sequence : {sequence}")
        ])
    return None


#Help callback
@app.callback(
    Output("help-modal", "is_open"),
    Output("modal-help-text", "children"),
    Input("help-min-node-size", "n_clicks"),
    Input("help-min-selected", "n_clicks"),
    Input("help-tolerance", "n_clicks"),
    Input("help-region-gap", "n_clicks"),
    Input("help-deletion", "n_clicks"),
    Input("help-genome_ref-dropdown", "n_clicks"),
    Input("close-help", "n_clicks"),
    State("help-modal", "is_open"),
)
def toggle_help_modal(n1, n2, n3, n4, n5, n6, n_close, is_open):
    ctx = callback_context

    if not ctx.triggered:
        return is_open, ""

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    help_messages = {
        "help-min-node-size": "The nodes with a size below this value won't be detect by the process.",
        "help-min-selected": "To detect a node as shared, the process requires that a node contains at least [this value * selected haplotypes number / 100] (rounded at tge bottom) haplotypes in the selected haplotypes list. If set to zéro, it will require at least one of the selected haplotypes.",
        "help-tolerance": "Nodes that contains more than [this value * number of haplotypes sharing the node / 100] (rounded up) haplotypes not in the selected list won't be detected.",
        "help-region-gap": "All the nodes detected will be grouped in larger regions. If two nodes are separated by less thant this value (in bp) they will be grouped in the same region.",
        "help-deletion": "If checked, the process will look for deletion : it looks for nodes with minimal size and with the minimum haplotypes (defined by min percent) of the selected list. For these nodes it will look for a following node containing all the haplotypes of the previous node excepted the nodes selected and another following node with all the haplotypes of the previous node.",
        "help-genome_ref-dropdown": "Select the genome for which you want to view the results and obtain annotations. If no genome is selected, the result will be the first genome found with annotations. If there are no annotations, it will be the first genome found."
    }

    if trigger_id == "close-help":
        return False, ""

    help_text = help_messages.get(trigger_id, "No help for this element.")
    return True, help_text




