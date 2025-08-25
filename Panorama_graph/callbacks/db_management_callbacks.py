#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 22:03:10 2025

@author: fgraziani
"""

from dash import html, Output, Input, State, no_update, dcc, ctx, no_update, exceptions


import os
import sys
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if root_path not in sys.path:
    sys.path.append(root_path)

from app import app, DB_VERSION
from neo4j_requests import *
from neo4j_DB_construction import *
from neo4j_container_management import *
from config import *
import base64
import io
import shutil


PREFIX_CONTAINER_NAME = "DB_"+ DB_VERSION + "_"

success_style = {"color": "green", "marginTop": "10px"}
warning_style = {"color": "orange", "marginTop": "10px"}
error_style = {"color": "red", "marginTop": "10px"}


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GFA_FOLDER = os.path.join(PROJECT_ROOT, "data", "gfa")
DATA_FOLDER = os.path.join(PROJECT_ROOT, "data", "data")
ANNOTATIONS_FOLDER = os.path.join(PROJECT_ROOT, "data", "annotations")

BATCH_SIZE = 2000000


def get_container_name_no_prefix(container_name):
    return re.sub(r'^DB_.[^_]+_', '',container_name)


############# Data callbacks#################


@app.callback(
    Output("upload-gfa-text", "children", allow_duplicate=True),
    Input("upload-gfa", "filename"),
    prevent_initial_call=True
)
def update_gfa_upload_text(filenames):
    if filenames:
        if isinstance(filenames, str):
            filenames = [filenames]
        return f"✅ Selected files: {', '.join(filenames)}"
    return ["Drag and Drop or ", html.A("Select one or more GFA Files")]


@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Output("upload-gfa-text", "children", allow_duplicate=True),
    Input("btn-load-gfa", "n_clicks"),
    State("upload-gfa", "filename"),
    State("db-chromosome-input", "value"),
    prevent_initial_call=True
)
def on_click_load_gfa(n_clicks, gfa_file_names, chromosome_file):
    if not gfa_file_names:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), [
            "Drag and Drop or ", html.A("Select one or more GFA Files")
        ]

    # Force to list
    if isinstance(gfa_file_names, str):
        gfa_file_names = [gfa_file_names]

    # ✅ Check all files have .gfa extension
    invalid_files = [f for f in gfa_file_names if not f.lower().endswith(".gfa")]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.", style=error_style), [
            "Drag and Drop or ", html.A("Select one or more GFA Files")
        ]
    if len(gfa_file_names) == 1 and not chromosome_file:
        list_chromosome_file = [None]
    else:
        if len(gfa_file_names) == 1:
            list_chromosome_file = [chromosome_file]
        else:
            list_chromosome_file = []
            for f in gfa_file_names:
                if "_" in gfa_file_name:
                    chromosome = gfa_file_name[:-4].split("_")[-1]
                else:
                    chromosome = gfa_file_name[:-4]

                chromosome = chromosome.lower().removeprefix("chr")
                chromosome = chromosome.lstrip("0")
                list_chromosome_file.append(chromosome)

    

    for i in range(len(gfa_file_names)):
        start_time = time.time()
        file_name = gfa_file_names[i]
        chromosome_file = list_chromosome_file[i]
        if chromosome_file != "" :
            file_path = os.path.join(GFA_FOLDER, file_name)
            load_sequences(file_path, chromosome_file=chromosome_file, create=True)
            load_gfa_data_to_neo4j(file_path, chromosome_file = chromosome_file, batch_size = BATCH_SIZE, start_chromosome = None, create = True, haplotype = True, create_only_relations = False)
        print(f"Graph from {file_name} loaded in {time.time() - start_time:.2f} s")

    create_indexes(base=False, extend=True, genomes_index=True)
    print("✅ All GFA files loaded.")

    return html.Div(f"✅ GFA files loaded successfully: {', '.join(gfa_file_names)}", style=success_style), [
        "Drag and Drop or ", html.A("Select one or more GFA Files")]


@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Output("upload-gfa-text", "children", allow_duplicate=True),
    Input("btn-csv-import", "n_clicks"),
    State("upload-gfa", "filename"),
    State("db-chromosome-input", "value"),
    prevent_initial_call=True
)
def on_click_csv_import(n_clicks, gfa_file_names, chromosome_file):
    print("generate import csv")
    if not gfa_file_names:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), [
            "Drag and Drop or ", html.A("Select one or more GFA Files")
        ]

    # Force to list
    if isinstance(gfa_file_names, str):
        gfa_file_names = [gfa_file_names]

    # ✅ Check all files have .gfa extension
    invalid_files = [f for f in gfa_file_names if not f.lower().endswith(".gfa")]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.", style=error_style), [
            "Drag and Drop or ", html.A("Select one or more GFA Files")
        ]
    if len(gfa_file_names) == 1 and not chromosome_file:
        list_chromosome_file = [None]
    else:
        if len(gfa_file_names) == 1:
            list_chromosome_file = [chromosome_file]
        else:
            list_chromosome_file = []
            for f in gfa_file_names:
                if "_" in gfa_file_name:
                    chromosome = gfa_file_name[:-4].split("_")[-1]
                else:
                    chromosome = gfa_file_name[:-4]

                chromosome = chromosome.lower().removeprefix("chr")
                chromosome = chromosome.lstrip("0")
                list_chromosome_file.append(chromosome)

    for i in range(len(gfa_file_names)):
        start_time = time.time()
        file_name = gfa_file_names[i]
        chromosome_file = list_chromosome_file[i]
        if chromosome_file != "" :
            file_path = os.path.join(GFA_FOLDER, file_name)
            load_gfa_data_to_csv(file_path, import_dir="./data/import", chromosome_file = chromosome_file, chromosome_prefix = False, batch_size = 2000000, start_chromosome = None, haplotype = True)
        print(f"CSV generation from {file_name} loaded in {time.time() - start_time:.2f} s")

    print("✅ All GFA files loaded.")

    return html.Div(f"✅ GFA files loaded successfully: {', '.join(gfa_file_names)}", style=success_style), [
        "Drag and Drop or ", html.A("Select one or more GFA Files")]

    
@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Input("btn-create-index", "n_clicks"),
    prevent_initial_call=True
)
def on_click_create_index(n_clicks):
    print("create indexes")
    create_indexes(base=True, extend=True, genomes_index=True)
    print("Indexes created")
    return html.Div(f"✅ Indexes creation command successfully done.", style=success_style)




############# Annotations callbacks#################


@app.callback(
    Output("annotation-message", "children", allow_duplicate=True),
    Input("dropdown-genome", "value"),
    prevent_initial_call=True
)
def clear_genome_error_on_selection(genome):
    if genome:
        return ""
    dash.exceptions.PreventUpdate 

@app.callback(
    Output("upload-annotation-text", "children", allow_duplicate=True),
    Input("upload-annotation", "filename"),
    prevent_initial_call=True
)
def update_annotation_upload_text(filename):
    if filename:
        return f"✅ File selected: {filename}"
    return ["Drag and Drop or ", html.A("Select an Annotation File")]

@app.callback(
    Output("annotation-message", "children"),
    Output("upload-annotation-text", "children", allow_duplicate=True),
    Input("btn-load-annotation", "n_clicks"),
    State("upload-annotation", "filename"),
    State("dropdown-genome", "value"),
    prevent_initial_call=True
)
def on_click_load_annotation(n_clicks, annotation_file_name, genome):
    print("load annotations")
    if not annotation_file_name:
        return html.Div("❌ Please select an annotation file before loading.", style=error_style), ["Drag and Drop or ", html.A("Select a GFA File")]
    if not genome or genome == "":
        return html.Div("❌ Please select a reference genome.", style=error_style), no_update
    state_index = check_state_index("NodeIndex"+genome+"_position")
    if state_index is None:
        return html.Div(f"❌ Index {state_index} has not been created, please create index before loading annotations.", style=error_style), no_update
    if int(state_index) != 100:
        return html.Div(f"❌ Index {state_index} is not completly created (creation state : {state_index}%). Please wait until this index has been created.", style=warning_style), no_update
    file = os.path.join(ANNOTATIONS_FOLDER, annotation_file_name)
    index_time = time.time()
    load_annotations_neo4j(file, genome_ref = genome, single_chromosome = None)
    annotation_time = time.time()
    print("Annotations loaded in " + str(annotation_time-index_time) + " s.")
    creer_relations_annotations_neo4j(genome)
    annotation_relation_time = time.time()
    print("annotations loaded in " + str(annotation_relation_time-annotation_time) + " s.")
    return html.Div(f"✅ Annotation '{annotation_file_name}' loaded for genome '{genome}'.", style=success_style), ["Drag and Drop or ", html.A("Select a GFA File")]

############# Delete data callbacks#################

@app.callback(
    Output("delete-confirmation", "children"),
    Input("btn-delete", "n_clicks"),
    prevent_initial_call=True
)
def delete_data_ask_confirmation(n_clicks):
    data_dir =  os.path.join(DATA_FOLDER, "databases/neo4j")
    transactions_dir = os.path.join(DATA_FOLDER,"transactions/neo4j")
    if n_clicks > 0:
        return html.Div([
            html.Div("⚠️ Confirm: this operation will delete all data in " + str(data_dir) + " and " + str(transactions_dir)),
            html.Button("Confirm Delete", id="btn-confirm-delete", n_clicks=0, style={"marginTop": "5px", "color": "white", "backgroundColor": "red"})
        ])
    return ""

@app.callback(
    Output("delete-message", "children"),
    Output("delete-confirmation", "children", allow_duplicate=True),
    Input("btn-confirm-delete", "n_clicks"),
    prevent_initial_call=True
)
def confirm_delete_data(n_clicks):
    if not n_clicks:
        raise exceptions.PreventUpdate
    data_dir =  os.path.join(DATA_FOLDER, "databases/neo4j")
    transactions_dir = os.path.join(DATA_FOLDER,"transactions/neo4j")
    stop_container()
    print("Deleting " + str(data_dir) + " and " + str(transactions_dir) + " directories.")
    print("exists : " +str(os.path.exists(data_dir)))
    try:
        if os.path.exists(data_dir):
            shutil.rmtree(data_dir)
            os.makedirs(data_dir)
            shutil.rmtree(transactions_dir)
            os.makedirs(transactions_dir)
            return html.Div("✅ All data deleted successfully.", style=success_style), ""
        else:
            return html.Div("ℹ️ No data to delete.", style=warning_style), ""
    except Exception as e:
        return html.Div(f"❌ Error while deleting data: {str(e)}", style=error_style), ""

############# Container callbacks#################

@app.callback(
    Output('container-name-label', 'children'),
    Output('container-name-input', 'value'),
    Output('db-management-page-store', 'data'),
    Input('db-management-page-store', 'data')
)
def update_label(data):
    if data is None :
        data = {}
    if "container_name" not in data :
        conf = load_config_from_json()
        if not conf:
            return f'No conf file found. Use "create new DB" procedure to generate it.', "container_name", data
        else:
            container_name = conf.get("container_name")
            data['container_name'] = get_container_name_no_prefix(container_name)
            return f'Container name : {container_name}', get_container_name_no_prefix(container_name), data
    else:
        container_name = PREFIX_CONTAINER_NAME+data['container_name']
        return f"Container name : {container_name}", get_container_name_no_prefix(container_name), data
    
    
@app.callback(
    Output("create-db-confirmation", "children"),
    Input("btn-create-db", "n_clicks"),

    prevent_initial_call=True
)
def create_db_ask_confirmation(n_clicks):
    if n_clicks > 0:
        return html.Div([
            html.Div("⚠️ Confirm: this operation will delete all data in " + str(DATA_FOLDER)),
            html.Button("Confirm creation of new DB", id="btn-confirm-create-db", n_clicks=0, style={"marginTop": "5px", "color": "white", "backgroundColor": "red"})
        ])
    return ""

@app.callback(
    Output("create-db-message", "children"),
    Output("create-db-confirmation", "children", allow_duplicate=True),
    Output('db-management-page-store', 'data', allow_duplicate=True),
    Output('dropdown-genome', 'options'),
    Input("btn-confirm-create-db", "n_clicks"),
    State("container-name-input","value"),
    State("docker-image-dropdown","value"),
    State('db-management-page-store', 'data'),
    State('dropdown-genome', 'options'),
    prevent_initial_call=True
)
def confirm_create_db(n_clicks,container_name, docker_image, data, options):
    print("container name : " + container_name)
    if not n_clicks:
        raise exceptions.PreventUpdate
    if container_name is None or container_name == "":
        if 'container_name' in data and data["container_name"] is not None and data["container_name"] != "":
            container_name = data["container_name"]
        else:
            return html.Div("❌ No container name", style=error_style), "", data, options
    else :
        data['container_name'] = container_name
        container_name=PREFIX_CONTAINER_NAME+container_name
    try:
        create_db(container_name, docker_image)
        genomes = get_genomes()
        options = []
        if genomes is not None :
            options = [{"label": genome, "value": genome} for genome in genomes]
        return html.Div("✅ DB successfully created.", style=success_style), "", data, options
    except Exception as e:
        return html.Div(f"❌ Error while creating database: {str(e)}", style=error_style), "", data, options
    

@app.callback(
    Output("create-db-message", "children", allow_duplicate=True),
    Input("btn-dump-db", "n_clicks"),
    State("docker-image-dropdown","value"),
    State('db-management-page-store', 'data'),
    prevent_initial_call=True
)
def dump_db(n_clicks, docker_image, data):
    if not n_clicks:
        raise exceptions.PreventUpdate
        

    if 'container_name' in data and data["container_name"] is not None and data["container_name"] != "":
        container_name = PREFIX_CONTAINER_NAME+data["container_name"]
    else:
        return html.Div("❌ No container name, you have to create the databse before", style=error_style)

    try:
        dump_db(container_name, docker_image=DOCKER_IMAGE)
        return html.Div("✅ DB successfully dumped.", style=success_style)
    except Exception as e:
        return html.Div(f"❌ Error while dumping database: {str(e)}", style=error_style)