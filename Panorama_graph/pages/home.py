import plotly.graph_objects as go


from neo4j_requests import *

import dash
from dash import Dash, dcc, html, Input, Output, State, callback_context, ctx, MATCH, ALL
import dash_cytoscape as cyto
import dash_bootstrap_components as dbc
import matplotlib.colors as mcolors

import pandas as pd
import numpy as np
import itertools
import json

from app import app

cyto.load_extra_layouts()

def records_to_dataframe(nodes_data):
    rows = []
    for record in nodes_data:
        rows.append(nodes_data[record])
    return pd.DataFrame(rows)


def compute_stylesheet(color_number):
    if color_number > 1:
        stylesheet =  [ {
             'selector': 'node',
             'style': {
                 'label': 'data(label)',
                 'min-zoomed-font-size': 1
                 
                 }
             },
            {
            'selector':'edge',
            'style':{
            'curve-style': 'unbundled-bezier',   
            'control-point-weights': [0.5],
            'target-arrow-shape': 'none',
            'control-point-distances':[1],
            'opacity': 1
            }
            }]
        for i in range(1, color_number+1):
            sign = 1 if i % 2 == 0 else -1
            distance = sign * (20 + 10 * (i //2))
            distance = 1 + 5*i
            stylesheet.append({
                'selector' : f'.offset-{i}',
                'style':{'control-point-distances': [distance], 'opacity': 0.6,}
                })
    else:
        stylesheet =  [
            {
                'selector': 'node',
                'style': {
                    'background-color': '#BFD7B5',
                    'label': 'data(label)',
                    'min-zoomed-font-size': 1
                    
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': '#A3C4BC'
                }
            }
            ]
    return stylesheet


def flow_to_rgb(flow):
    r = int(255 * flow)
    g = int(0)
    b = int(255 * (1 - flow))
    return f'rgb({r},{g},{b})'


def get_color_palette(n):
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab20")
    return [f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})' for r, g, b, _ in cmap(np.linspace(0, 1, n))]


def compute_graph_elements(data, selected_genomes, size_min, specifics_genomes=None, color_genomes = [], x_max=1000, y_max=1000, labels=True):
    all_genomes = get_genomes()
    all_chromosomes = get_chromosomes()
    if data != None and len(data) > 0:
        df = records_to_dataframe(data)
        df = df[df["size"] >= size_min].copy()
        df = df[df["genomes"].apply(lambda g: any(x in selected_genomes for x in g))].copy()
        
    
        def moyenne_position(row):
            positions = [row.get(f"{g}_node") for g in row["genomes"] if f"{g}_node" in row]
            positions = [p for p in positions if p is not None]
            return np.mean(positions) if positions else 0
    
        df["mean_pos"] = df.apply(moyenne_position, axis=1)
        x_min, x_max_data = df["mean_pos"].min(), df["mean_pos"].max()
        df["x"] = ((df["mean_pos"] - x_min) / (x_max_data - x_min + 1e-6)) * x_max
        
        df["genome_key"] = df["genomes"].apply(lambda g: "".join(sorted(g)))
        genome_keys = sorted(df["genome_key"].drop_duplicates(), key=lambda x: (-len(x), x))
        y_positions = {k : 1 for k in enumerate(genome_keys)}
        df["y"] = df["genome_key"].map(y_positions)
        
        color_map = {k: c for k, c in zip(genome_keys, get_color_palette(len(genome_keys)))}
        
        nodes = []
        size_max = df['size'].max()
        size_min = df['size'].min()
        size_max_noeud = 100

        for _, row in df.iterrows():
            nodes.append({
                'data': {'id': row['name'], 'size' : row['size'], 'flow':row['flow'], 'genomes':row['genomes'], 'chromosome':row['chromosome'], 'annotations':row['annotations']},
                'position': {'x': row['x'], 'y': row['y']},
                'style': {
                    'background-color':flow_to_rgb(row['flow']),
                    'shape':'circle',
                    'width': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
                    'height': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
                }
            })
        
        edges = []
        edges_dict = {}
        for genome in selected_genomes:
            nodes_g = df[df["genomes"].apply(lambda g: genome in g)].copy()
            col = f"{genome}_position"
            if col not in nodes_g.columns:
                continue
            nodes_g = nodes_g[nodes_g[col].notnull()].copy()
            if nodes_g.empty:
                continue
            nodes_g = nodes_g.sort_values(by=col, ascending=True).reset_index(drop=True)
            for i in range(len(nodes_g) - 1):
                source = nodes_g.loc[i, 'name']
                target = nodes_g.loc[i + 1, 'name']
                edge_key = tuple(sorted([source, target]))
                if edge_key not in edges_dict:
                    edges_dict[edge_key] = {}
                    edges_dict[edge_key]["flow"] = 1
                    edges_dict[edge_key]["annotations"] = set()
                    edges_dict[edge_key]["genomes"] = [genome]
                else :
                    edges_dict[edge_key]["flow"] += 1
                    if genome not in edges_dict[edge_key]["genomes"] :
                        edges_dict[edge_key]["genomes"].append(genome)
                
                for a in nodes_g.loc[i, 'annotations']+nodes_g.loc[i + 1, 'annotations'] :
                    edges_dict[edge_key]["annotations"].add(a)
        
        
        #Specific colors for a selected set of genomes
        colored_genomes = {}
        for g,c in zip(all_genomes, color_genomes):
            if c != '#000000' :
                colored_genomes[g] = c
        for (source, target), dic in edges_dict.items(): 
            link_color = 'gray'
            flow = dic["flow"]
            if specifics_genomes is not None and len(specifics_genomes) > 0:
                #Uncomment to get the paths shared by the wholes selected genomes and only these genomes
                #if all(g in specifics_genomes for g in dic["genomes"]) and all(g in dic["genomes"] for g in specifics_genomes):
                #Path are colored in red if it contains at leat one of seleted genomes but no unselected genome
                if all(g in specifics_genomes for g in dic["genomes"]):
                    link_color = 'red'
                    flow = len(all_genomes)
            label = ""
            if labels:

                first_label = True
    
                for a in dic["annotations"]:
                    if first_label :
                        label += str(a)
                        first_label = False
                    else :
                        label += ", " + str(a)
            edges.append({
                'data': {
                    'source': source,
                    'target': target,
                    'flow': flow,
                    'genomes' : dic["genomes"]
                },
                'style': {
                    'line-color': link_color,
                    'label':label,
                    'text-rotation':'autorotate',
                    'width': (flow+int(0.2*len(all_genomes)))/len(all_genomes)*10
                    }
            })

            if len(colored_genomes) > 0:
                i = 0
                for g in list(colored_genomes.keys()):
                    if g in dic["genomes"]:
                        sign = 1 if i % 2 == 0 else -1
                        distance = sign * (20 + 10 * (i //2))
                        edges.append({
                            'data': {
                                'source': source,
                                'target': target,
                                'flow':1,
                                'genomes' : [g]
                            },
                            'classes':f'offset-{i}',
                            'style':{
                                'line-color':colored_genomes[g],
                                'width':4

                                }
                            })
                        i+=1
        print("nb nodes : " + str(len(nodes)) + " - nb edges : " + str(len(edges)))                
        return nodes + edges
    else :
        return []


def layout(data=None, initial_size_limit = 10):
    all_genomes = get_genomes()
    all_chromosomes = get_chromosomes()
    if data != None :
        elements = compute_graph_elements(data,all_genomes, initial_size_limit, [], [])
    else :
        elements = []
    print("nombre d'élements : " + str(len(elements)))
    max_size = 500
    
    return html.Div([

        #Upper block : settings
        html.Div([
            #Left block
            html.Div([
                html.Div([
                    html.Div([
                        html.Label("Genome", style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='genomes-dropdown',
                            options=[{'label': genome, 'value': genome} for genome in all_genomes],
                            value='HER',
                            style={'width': '200px'}
                        )
                    ], style={'marginRight': '30px'}),
                    html.Div([
                        html.Label("Ref genome", style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='genomes-ref-dropdown',
                            options=[{'label': genome, 'value': genome} for genome in all_genomes],
                            value='HER',
                            style={'width': '200px'}
                        )
                    ], style={'marginRight': '30px'}),
                    html.Div([
                        html.Label("Chromosome", style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='chromosomes-dropdown',
                            options=[{'label': chrom, 'value': chrom} for chrom in all_chromosomes],
                            placeholder="Choisissez un chromosome",
                            value="1",
                            style={'width': '200px'}
                        )
                    ])
                ], style={'display': 'flex', 'padding': '20px', 'border': '1px solid #ccc', 'minWidth': '300px', 'boxSizing': 'border-box'}),
            
        
                html.Div([
                    html.H5("Search region", style={'textAlign': 'left', 'marginBottom':'15px'}),
                    html.Div([
                        html.Label("Start : "),
                        dcc.Input(id='start-input', type='number', style={'width': '100px', 'marginRight': '10px'}),
                        html.Label("End : "),
                        dcc.Input(id='end-input', type='number', style={'width': '100px', 'marginRight': '20px'})
                    ], style={'marginBottom': '10px'}),
                    html.Div([
                        html.Label("Gene name : "),
                        dcc.Input(
                            id='genename-input',
                            type='text',
                            placeholder='Nom de gene',
                            debounce=True,  # Déclenche le callback uniquement après validation (Entrée ou perte de focus)
                            style={'marginRight':"10px"}
                        ),
                        html.Label("Gene id : "),
                        dcc.Input(
                            id='geneid-input',
                            type='text',
                            placeholder='id de gene',
                            debounce=True  # Déclenche le callback uniquement après validation (Entrée ou perte de focus)
                        )
                    ], style={'marginBottom':'20px'}),
                    html.Button('Search', id='search-button', n_clicks=0, style={'marginTop': '10px'}),
                    html.Div([
                        html.Label("Genomes to visualize :", style={'marginBottom': '5px'}),
                        dcc.Checklist(
                            id="genome_selector",
                            options=[{"label": g, "value": g} for g in all_genomes],
                            value=all_genomes,
                            inline=True
                        )
                        
                    ]),
                    dcc.Loading(id="loading-spinner", type="circle", children=html.Div(id="output-zone"))
                ], 
                style={'marginBottom': '20px'}
                )
            ], style={'flex': '1', 'padding': '20px', 'border': '1px solid #ccc'}),
            
            #Right block
            html.Div([
                
                html.Div([
                    html.Label("Minimal size of nodes :"),
           
                    dcc.Slider(
                        id='size_slider',
                        min=0,
                        max=size_max,
                        step=int(100/size_max),
                        value=10,
                        tooltip={"placement": "bottom", "always_visible": False},
                    ),
                    
                    html.Div(id='size_stats', style={'marginTop': '10px'})
                    
                    ]),
                html.Div([
                    html.Div(id='size-output', children='Min node size : 10', style={'margin':'10px'}),
                    html.H4("Layout", style={'marginLeft':'50px'}),
                    dcc.Dropdown(
                        id='layout-dropdown',
                        options=[
                            {'label': 'fcose', 'value': 'fcose'},
                            {'label': 'dagre', 'value': 'dagre'}
                        ],
                        value='fcose',
                        clearable=False,
                        style={'width': '120px', 'display': 'inline-block','marginRight':'50px'}
                    ),
                    dcc.Checklist(
                        options=[{'label': 'Hide labels', 'value': 'hide'}],
                        id='show-labels'
                    ),
                ], style={'display': 'flex', 'alignItems': 'center', 'gap': '8px'}),
                html.Div(id='nb-noeuds', style={'margin':'10px'}),
                
                dcc.Checklist(
                    options=[{'label': 'Search shared path', 'value': 'shared'}],
                    id='shared-mode',
                    style={'marginBottom': '20px'}
                ),
                html.Div([

                    
                        html.Div([
                            html.Label("Vizualize shared paths :", style={'marginBottom': '20px'}),
                            dcc.Checklist(
                                id="specific-genome_selector",
                                options=[{"label": g, "value": g} for g in all_genomes],
                                value=[],
                                inline=True
                            )
                            
                        ],id='shared-checklist-container'),
                    html.Div([
                            html.Div([
                                html.Label(s),
                                dbc.Input(
                                    id={'type': 'color-picker', 'index': i},
                                    type='color',
                                    value="#000000",
                                    style={'width': '25px','height': '25px', 'marginLeft': '10px'}
                                )
                                
                            ], style={'marginRight': '10px'}) for i, s in enumerate(all_genomes)
                    ], style={'display':'flex',
                              'flexWrap':'wrap',
                              'gap':'5px'}, id='color-picker-container')
                    
                    ], id='sample-controls'),
                html.Button("Update graph", id="update-btn", n_clicks=0, style={'marginTop' : '10px'}),
                
                html.Div(html.H4(id='node-info', style={'margin':'10px'})),
                html.Div(html.Label("Annotations :", style={'marginBottom': '5px'})),
                html.Div(html.H4(id='annotations-info', style={'margin':'10px'}))
            ], style={'flex': '1','padding': '20px','border': '1px solid #ccc','marginLeft': '20px','minWidth': '300px',
            'boxSizing': 'border-box','display': 'flex', 'flex-direction': 'column'})
         ], style={
        'display': 'flex',
        'flex-direction': 'row',
        'flexWrap': 'wrap',
        'justify-content': 'space-between'
        }),   
        
        #Graph block
        cyto.Cytoscape(
            id='graph',
            #layout is important to get good visualization result
            #There are many algorithms : cose, cose-bilkent-fcose, euler, dagre, etc.
            #fcose seems to be the most performant
            #dagre is usefull to get a linear representation
            #Le layout est super important pour une bonne visu du graphe, pour le moment fcose semble le plus adapté
            layout={
                'name': 'fcose',
                'maxIterations':100000, 
                'maxSimulationTime':5000,
                #'nodeRepulsion': 10000,
                #'gravity': 0.1,
                #'gravityRangeCompound': 1.5,
                #'idealEdgeLength': 100,
                #'componentSpacing': 100,
                #'nodeDimensionsIncludeLabels': True,
                #'edgeElasticity': 0.1,
                #'nestingFactor': 0.8,
                #'tile': True,
                'quality': "proof",
                'fit': True
                },
            style={'width': '100%', 'height': '1000px'},
            elements=elements,
            #minZoom=0.1,
            #maxZoom=5,
            zoomingEnabled=True,
            userZoomingEnabled=True,
            #userPanningEnabled=True,
            wheelSensitivity=0.1,
            boxSelectionEnabled= True,
            autoungrabify=False,
            stylesheet=compute_stylesheet(0),
            
        )
        
    ])
    
#Callback to get nodes or link info when clicking on it
@app.callback(
Output('node-info', 'children'),
Input('graph', 'tapNodeData'),
Input('graph', 'tapEdgeData')
)

def display_element_data(node_data, edge_data):
    triggered_id = ctx.triggered_id

    if triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapEdgeData' and edge_data:
        return (
            f"Lien sélectionné : {edge_data.get('source')} → {edge_data.get('target')}\n"
            f"• Flow : {edge_data.get('flow')}\n"
            f"• Génomes : {', '.join(edge_data.get('genomes', []))}"
        )
    elif triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapNodeData' and node_data:
        return (
            f"Selected node : {node_data.get('label', node_data.get('id'))}\n"
            f"• Size : {node_data.get('size')}\n"
            f"• Flow : {node_data.get('flow')}\n"
            f"• Genomes : {', '.join(node_data.get('genomes', []))}"
            f"• Annotations : {', '.join(node_data.get('annotations', []))}"
        )
    return "Click on a node or link to display data."
   
#Main callback to update graph when changing size, or selecting genomes, etc.
@app.callback(
    Output("graph", "elements"),
    Output("nb-noeuds", 'children'),
    Output("size-output", 'children'),
    Output('shared_storage_nodes', 'data', allow_duplicate=True),
    Output('output-zone', 'children'),
    Output('annotations-info', 'children'),
    Output('graph', 'stylesheet'),
    Input('genome_selector', 'value'), 
    State('shared-mode', 'value'), 
    State('specific-genome_selector', 'value'),
    State({'type':'color-picker', 'index': ALL}, 'value'),
    Input('show-labels', 'value'),
    Input('update-btn', 'n_clicks'), 
    Input('size_slider', 'value'),
    Input('search-button', 'n_clicks'),
    State('start-input', 'value'),
    State('end-input', 'value'),
    State('genename-input', 'value'),
    State('geneid-input', 'value'),
    State('genomes-dropdown', 'value'),
    State('genomes-ref-dropdown', 'value'),
    State('chromosomes-dropdown', 'value'),
    State('shared_storage', 'data'),
    State('shared_storage_nodes', 'data'),
    prevent_initial_call=True
)
def update_graph(selected_genomes, shared_mode, specifics_genomes, color_genomes, show_labels, update_n_clicks, size_slider_val, n_clicks, start, end, gene_name, gene_id, genome, genome_ref, chromosome,data_storage, data_storage_nodes):
    ctx = dash.callback_context
    print("update graph : " + str(ctx.triggered[0]['prop_id']))
    stylesheet = []
    if shared_mode and 'shared' in shared_mode:
        specifics_genomes_list = specifics_genomes
        color_genomes_list = []
    else:
        specifics_genomes_list = []
        color_genomes_list = color_genomes
    labels = True
    if show_labels and 'hide' in show_labels:
        labels = False

    if ctx.triggered and len(ctx.triggered[0]['prop_id'].split('.')) > 0:
        input_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if input_id == "search-button" and n_clicks > 0:
            if start is not None :
                new_data = get_nodes_by_region(genome, chromosome=chromosome, start=start, end=end)
                data_storage_nodes = new_data
                print("len new_data : " + str(len(new_data)))
            else : 
                if gene_name is not None :
                    new_data = get_nodes_by_gene(genome_ref, gene_name=gene_name,chromosome=chromosome)
                    data_storage_nodes = new_data
            elements = compute_graph_elements(new_data, selected_genomes, size_slider_val, specifics_genomes_list, color_genomes_list, labels=labels)
        else:
            elements = compute_graph_elements(data_storage_nodes, selected_genomes, size_slider_val, specifics_genomes_list, color_genomes_list, labels=labels)
    else :
        elements = compute_graph_elements(data_storage_nodes, selected_genomes, size_slider_val, specifics_genomes_list, color_genomes_list, labels=labels)

    defined_color = 0
    for c in color_genomes:
        if c != "#000000":
            defined_color += 1
    stylesheet = compute_stylesheet(defined_color)
    count = len(elements)
    annotations = ""
    set_annot = set()
    if data_storage_nodes != None:
        for n in data_storage_nodes: 
            if "annotations" in data_storage_nodes[n]:
                for a in data_storage_nodes[n]["annotations"]:
                    set_annot.add(a)
    for a in set_annot:
        annotations += str(a) + "\n"

    return elements, f"{count} displayed nodes", f"Min node size  : {size_slider_val}",data_storage_nodes,"",annotations, stylesheet


#color picker
@app.callback(
    Output('color-picker-container', 'style'),
    Output('shared-checklist-container', 'style'),
    Input('shared-mode', 'value')
)
def toggle_inputs(shared_mode):
    if shared_mode and 'shared' in shared_mode:
        return {'display':'none'},{'display':'block'}
    else:
        return {'display':'flex','flexWrap':'wrap'},{'display':'none'}


#Algorithm cytoscape choice
@app.callback(
    Output('graph', 'layout'),
    Input('layout-dropdown', 'value')
)
def toggle_layout(layout_choice):
    if layout_choice and 'dagre' in layout_choice:
        return {
            'name': 'dagre',
            'rankDir': "RL",
            'nodeDimensionsIncludeLabels': True
        }
    else:
        return {
            'name': 'fcose',
            'maxIterations':100000, 
            'maxSimulationTime':5000,
            #'nodeRepulsion': 10000,
            #'gravity': 0.1,
            #'gravityRangeCompound': 1.5,
            #'idealEdgeLength': 100,
            #'componentSpacing': 100,
            #'nodeDimensionsIncludeLabels': True,
            #'edgeElasticity': 0.1,
            #'nestingFactor': 0.8,
            #'tile': True,
            'quality': "proof",
            'fit': True
            }
        


