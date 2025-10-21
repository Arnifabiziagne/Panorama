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

from app import *
from neo4j_requests import *
from neo4j_DB_construction import *
from neo4j_container_management import *
from config import *
import base64
import io
import shutil


PREFIX_CONTAINER_NAME = "DB_"+ DB_VERSION + "_"




PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GFA_FOLDER = os.path.join(PROJECT_ROOT, "data", "gfa")
DATA_FOLDER = os.path.join(PROJECT_ROOT, "data", "data")
ANNOTATIONS_FOLDER = os.path.join(PROJECT_ROOT, "data", "annotations")

DEFAULT_BATCH_SIZE = 2000000
MIN_BATCH_SIZE = 1000


def get_container_name_no_prefix(container_name):
    return re.sub(r'^DB_.[^_]+_', '',container_name)


############# Data callbacks#################



@app.callback(
        Output('upload-gfa-output', 'children'),
        Input('upload-gfa-data', 'contents'),
        State('upload-gfa-data', 'filename'),
        State('upload-gfa-data', 'last_modified'),
        prevent_initial_call=True
    )
def save_uploaded_files(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        saved_files = []
        for content, filename in zip(list_of_contents, list_of_names):
            if not filename.endswith(".gfa"):
                continue
            try:
                data = content.encode("utf8").split(b";base64,")[1]
                file_path = os.path.join(GFA_FOLDER, filename)
                with open(file_path, "wb") as f:
                    f.write(base64.b64decode(data))
                saved_files.append(html.Li(f"File saved : {filename}"))
            except Exception as e:
                saved_files.append(html.Li(f"Error for saving file {filename} : {str(e)}"))

        return html.Ul(saved_files)
    return html.Div("No file.")


@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Output("gfa-files-selector", "value"),
    Input("btn-load-gfa", "n_clicks"),
    State('gfa-files-selector', 'value'),
    State("db-chromosome-input", "value"),
    State("db-batch-size-input", "value"),
    prevent_initial_call=True
)
def on_click_load_gfa(n_clicks, gfa_file_names, chromosome_file, batch_size=DEFAULT_BATCH_SIZE):
    if not gfa_file_names:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), no_update
    if batch_size < MIN_BATCH_SIZE:
        batch_size = MIN_BATCH_SIZE
    # Force to list
    if isinstance(gfa_file_names, str):
        gfa_file_names = [gfa_file_names]

    # ✅ Check all files have .gfa extension
    invalid_files = [f for f in gfa_file_names if not f.lower().endswith(".gfa")]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.", style=error_style), no_update
    if len(gfa_file_names) == 1 and not chromosome_file:
        list_chromosome_file = [None]
    else:
        if len(gfa_file_names) == 1:
            list_chromosome_file = [chromosome_file]
        else:
            list_chromosome_file = []
            for f in gfa_file_names:
                if "_" in f:
                    chromosome = f[:-4].split("_")[-1]
                else:
                    chromosome = f[:-4]

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
            load_gfa_data_to_neo4j(file_path, chromosome_file = chromosome_file, batch_size = batch_size, start_chromosome = None, create = True, haplotype = True, create_only_relations = False)
        print(f"Graph from {file_name} loaded in {time.time() - start_time:.2f} s")

    create_indexes(base=True, extend=True, genomes_index=True)
    print("✅ All GFA files loaded.")

    return html.Div(f"✅ GFA files loaded successfully: {', '.join(gfa_file_names)}", style=success_style), []


@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Output("gfa-files-selector", "value", allow_duplicate=True),
    Input("btn-csv-import", "n_clicks"),
    State('gfa-files-selector', 'value'),
    State("db-chromosome-input", "value"),
    State("db-batch-size-input", "value"),
    prevent_initial_call=True
)
def on_click_csv_import(n_clicks, gfa_file_names, chromosome_file, batch_size=DEFAULT_BATCH_SIZE):
    print("generate import csv")
    if not gfa_file_names:
        return html.Div("❌ Please select at least one GFA file before importing.", style=error_style), no_update
    
    if batch_size < MIN_BATCH_SIZE:
        batch_size = MIN_BATCH_SIZE
    # Force to list
    if isinstance(gfa_file_names, str):
        gfa_file_names = [gfa_file_names]

    # ✅ Check all files have .gfa extension
    invalid_files = [f for f in gfa_file_names if not f.lower().endswith(".gfa")]
    if invalid_files:
        return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gfa files.", style=error_style), no_update
    if len(gfa_file_names) == 1 and not chromosome_file:
        list_chromosome_file = [None]
    else:
        if len(gfa_file_names) == 1:
            list_chromosome_file = [chromosome_file]
        else:
            list_chromosome_file = []
            for f in gfa_file_names:
                if "_" in f:
                    chromosome = f[:-4].split("_")[-1]
                else:
                    chromosome = f[:-4]

                chromosome = chromosome.lower().removeprefix("chr")
                chromosome = chromosome.lstrip("0")
                list_chromosome_file.append(chromosome)

    for i in range(len(gfa_file_names)):
        start_time = time.time()
        file_name = gfa_file_names[i]
        chromosome_file = list_chromosome_file[i]
        if chromosome_file != "" :
            file_path = os.path.join(GFA_FOLDER, file_name)
            load_gfa_data_to_csv(file_path, import_dir="./data/import", chromosome_file = chromosome_file, chromosome_prefix = False, batch_size = batch_size, start_chromosome = None, haplotype = True)
        print(f"CSV generation from {file_name} loaded in {time.time() - start_time:.2f} s")
    print("✅ All GFA files loaded.")

    return html.Div(f"✅ GFA files loaded successfully: {', '.join(gfa_file_names)}", style=success_style), []

    
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


@app.callback(
    Output("gfa-message", "children", allow_duplicate=True),
    Input("btn-create-stats", "n_clicks"),
    prevent_initial_call=True
)
def on_click_create_index(n_clicks):
    print("create stats")
    create_stats_from_nodes()
    print("Stats created")
    return html.Div(f"✅ Stats creation command successfully done.", style=success_style)



############# Annotations callbacks#################


@app.callback(
        Output('upload-annotations-output', 'children'),
        Input('upload-annotations-data', 'contents'),
        State('upload-annotations-data', 'filename'),
        State('upload-annotations-data', 'last_modified'),
        prevent_initial_call=True
    )
def save_uploaded_files(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        saved_files = []
        for content, filename in zip(list_of_contents, list_of_names):
            if not filename.endswith(('.gff', '.gff3', '.gtf')):
                continue
            try:
                data = content.encode("utf8").split(b";base64,")[1]
                file_path = os.path.join(ANNOTATIONS_FOLDER, filename)
                with open(file_path, "wb") as f:
                    f.write(base64.b64decode(data))
                saved_files.append(html.Li(f"File saved : {filename}"))
            except Exception as e:
                saved_files.append(html.Li(f"Error for saving file {filename} : {str(e)}"))

        return html.Ul(saved_files)
    return html.Div("No file.")

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
    Output("annotation-message", "children"),
    Output("annotations-files-selector", "value"),
    Input("btn-load-annotations-with-link", "n_clicks"),
    #Input("btn-load-only-annotations", "n_clicks"),
    #Input("btn-link-annotations", "n_clicks"),
    State("dropdown-genome", "value"),
    State('annotations-files-selector', 'value'),
    prevent_initial_call=True
)
def on_click_load_annotation(n_clicks_load_all, genome, annotation_file_names):
    triggered_id = ctx.triggered_id
    annotation_time = time.time()
    print("annotations file names : " + str(annotation_file_names))
    print("triggered id : " + str(triggered_id))
    if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-load-only-annotations":
        if not annotation_file_names:
            return html.Div("❌ Please select an annotation file before loading.", style=error_style), no_update
        if not genome or genome == "":
            return html.Div("❌ Please select a reference haplotype.", style=error_style), no_update
        state_index = check_state_index("NodeIndex"+genome+"_position")
        if state_index is None:
            return html.Div(f"❌ Index {state_index} has not been created, please create index before loading annotations.", style=error_style), no_update
        if int(state_index) != 100:
            return html.Div(f"❌ Index {state_index} is not completly created (creation state : {state_index}%). Please wait until this index has been created.", style=warning_style), no_update

        # Force to list
        if isinstance(annotation_file_names, str):
            annotation_file_names = [annotation_file_names]

        # ✅ Check all files have .gtf / .gff3 / .gff extension
        invalid_files = [f for f in annotation_file_names if not f.lower().endswith(".gff") and not f.lower().endswith(".gff3") and not f.lower().endswith(".gtf")]
        if invalid_files:
            return html.Div(f"❌ Invalid file(s): {', '.join(invalid_files)}. Please select only .gff, .gff3 or gtf files.", style=error_style), no_update

        for f in annotation_file_names:
            print(f"Load annotations for file {f}")
            state_index = check_state_index("NodeIndex"+genome+"_position")
            if state_index is None:
                return html.Div(f"❌ Index {state_index} has not been created, please create index before loading annotations.", style=error_style), no_update
            if int(state_index) != 100:
                return html.Div(f"❌ Index {state_index} is not completly created (creation state : {state_index}%). Please wait until this index has been created.", style=warning_style), no_update
            file = os.path.join(ANNOTATIONS_FOLDER, f)
            annotation_time = time.time()
            load_annotations_neo4j(file, genome_ref = genome, single_chromosome = None)
        print("Annotations loaded in " + str(time.time()-annotation_time) + " s.")
    if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-link-annotations" :
        print("Link annotations")
        annotation_relation_time = time.time()
        creer_relations_annotations_neo4j(genome)
        print("Link annotations in " + str(time.time()-annotation_relation_time) + " s.")
    if triggered_id == "btn-load-annotations-with-link" or triggered_id == "btn-load-only-annotations":
        return html.Div(f"✅ Annotation '{annotation_file_names}' loaded for genome '{genome}'.", style=success_style), []
    else:
        return html.Div(f"✅ Annotation linked.", style=success_style),[]

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
        creation_mode = create_db(container_name, docker_image)
        #If creation by importing csv files it is necessary to create stats and indexes
        if creation_mode == "csv" :
            print("creating base indexes")
            create_indexes(base=True, extend=True, genomes_index=False)
            if check_state_index("NodeIndexChromosome") is not None:
                t = 0
                while int(check_state_index("NodeIndexChromosome")) < 100 and t < MAX_TIME_INDEX:
                    time.sleep(10)
                    t+=10
                print("creating stats")
                create_stats_from_nodes()
            print("creating other indexes")    
            create_indexes(base=False, extend=False, genomes_index=True)
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
def dump_db_callback(n_clicks, docker_image, data):
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