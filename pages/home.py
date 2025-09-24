import plotly.graph_objects as go


from neo4j_requests import *

import dash
from dash import Dash, dcc, html, Input, Output, State, callback_context, ctx, MATCH, ALL
from dash.exceptions import PreventUpdate
import dash_cytoscape as cyto
import dash_bootstrap_components as dbc
import matplotlib.colors as mcolors

import pandas as pd
import numpy as np
import itertools
import json

from app import *

import os
import base64

cyto.load_extra_layouts()

DEFAULT_SIZE_VALUE = 10
DEFAULT_SHARED_REGION_COLOR = "#008000"
DEFAULT_EXONS_COLOR = "#008000"
EXPORT_DIR = './export/graphs'


def records_to_dataframe(nodes_data):
    rows = []
    for record in nodes_data:
        rows.append(nodes_data[record])
    return pd.DataFrame(rows)


def compute_stylesheet(color_number):
    if color_number > 1:
        stylesheet = [
            {
            'selector': 'node',
            'style': {
                'label': 'data(label)',
                'background-color':'data(color)',
                'min-zoomed-font-size': 10,
                'text-opacity':1,
                'opacity':1,
                'width':'data(displayed_node_size)',
                'height':'data(displayed_node_size)',
                'z-index':9999

            }
        },
            {
            'selector': 'edge',
            'style': {
                'curve-style': 'unbundled-bezier',
                'control-point-weights': [0.5],
                'target-arrow-shape': 'none',
                'control-point-distances': [1],
                'opacity':0.9,
                'z-index':0
            },
        },
            {'selector': ':selected', 'style': {
                'background-color': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }}  
            
            ]
        for i in range(1, color_number+1):
            sign = 1 if i % 2 == 0 else -1
            distance = sign * (20 + 10 * (i // 2))
            distance = 1 + 5*i
            stylesheet.append({
                'selector': f'.offset-{i}',
                'style': {'control-point-distances': [distance], 'opacity': 0.6, }
            })
    else:
        stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'background-color':'data(color)',
                    'label': 'data(label)',
                    'min-zoomed-font-size': 10,
                    'opacity':1,
                    'text-opacity':1,
                    'width':'data(displayed_node_size)',
                    'height':'data(displayed_node_size)',
                    'z-index':9999

                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': '#A3C4BC',
                    'opacity':0.9,
                    'z-index':0
                }
            },
            {'selector': ':selected', 'style': {
                'background-color': 'grey',
                'line-color': '#FF4136',
                'border-width': 3,
                'border-color': 'black'
            }}
        ]
    return stylesheet


def flow_to_rgb(flow, node_style="default", exons_color=DEFAULT_EXONS_COLOR):
    if node_style=="exon":
        defined_color = hex_to_rgb_string(exons_color)
    else:
        r = int(255 * flow)
        g = int(0)
        b = int(255 * (1 - flow))
        defined_color = f'rgb({r},{g},{b})'
    return defined_color

def hex_to_rgb_string(hex_color):
    hex_color = hex_color.lstrip('#')
    if len(hex_color) != 6:
        raise ValueError(f"Invalid hex color: {hex_color}")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f'rgb({r}, {g}, {b})'


def get_color_palette(n):
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab20")
    return [f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})' for r, g, b, _ in cmap(np.linspace(0, 1, n))]



def compute_graph_elements(data, selected_genomes, size_min, all_genomes, all_chromosomes, specifics_genomes=None, color_genomes=[], x_max=1000, y_max=1000, labels=True, min_shared_genome=100, tolerance=0, color_shared_regions=DEFAULT_SHARED_REGION_COLOR, exons=False, exons_color=DEFAULT_EXONS_COLOR):
    if data != None and len(data) > 0:
        df = records_to_dataframe(data)
        df = df[df["size"] >= size_min].copy()
        df = df[df["genomes"].apply(lambda g: any(
            x in selected_genomes for x in g))].copy()

        def mean_position(row):
            positions = [row.get(f"{g}_node")
                         for g in row["genomes"] if f"{g}_node" in row]
            positions = [p for p in positions if p is not None]
            return np.mean(positions) if positions else 0

        df["mean_pos"] = df.apply(mean_position, axis=1)
        x_min, x_max_data = df["mean_pos"].min(), df["mean_pos"].max()
        df["x"] = ((df["mean_pos"] - x_min) /
                   (x_max_data - x_min + 1e-6)) * x_max

        df["genome_key"] = df["genomes"].apply(lambda g: "".join(sorted(g)))
        genome_keys = sorted(
            df["genome_key"].drop_duplicates(), key=lambda x: (-len(x), x))
        y_positions = {k: 1 for k in enumerate(genome_keys)}
        df["y"] = df["genome_key"].map(y_positions)

        color_map = {k: c for k, c in zip(
            genome_keys, get_color_palette(len(genome_keys)))}

        nodes = []
        size_max = df['size'].max()
        size_min = df['size'].min()
        size_max_noeud = 100
        
        for _, row in df.iterrows():
            node_style = "default"
            if exons :
                if "features" in row and "exon" in row['features']:
                    node_style="exon"
            node_color = flow_to_rgb(row['flow'],node_style,exons_color)
            displayed_node_size = (10+row['size']-size_min) / size_max*size_max_noeud+size_min
            main_style = {
                #'background-color': flow_to_rgb(row['flow'],node_style),
                'shape': 'circle',
                #'width': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
                #'height': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
                
            }
            degenerate_node_style = {
                #'background-color': flow_to_rgb(row['flow'],node_style),
                'shape': 'square',
                #'width': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
                #'height': (10+row['size']-size_min) / size_max*size_max_noeud+size_min,
            }

                # else:
                #     main_style["background-color"]="#000000"
                #     degenerate_node_style["background-color"]="#000000"
            if row['ref_node'] == row['name']:
                nodes.append({
                    'data': {'id':row['name'] ,'name':row['name'], 'displayed_node_size':displayed_node_size, 'ref_node': row['ref_node'], 'size': row['size'], 'flow': row['flow'], 'genomes': row['genomes'], 'chromosome': row['chromosome'], 'annotations': row['annotations'], 'features': row['features'], 'color':node_color},
                    'position': {'x': row['x'], 'y': row['y']},
                    'style': main_style
                })
            else:
                nodes.append({
                    'data': {'id':row['name'], 'name': row['name'],  'displayed_node_size':displayed_node_size, 'ref_node': row['ref_node'], 'size': row['size'], 'flow': row['flow'], 'genomes': row['genomes'], 'chromosome': row['chromosome'], 'annotations': row['annotations'], 'features': row['features'], 'color':node_color},
                    'position': {'x': row['x'], 'y': row['y']},
                    'style': degenerate_node_style
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
            nodes_g = nodes_g.sort_values(
                by=col, ascending=True).reset_index(drop=True)
            for i in range(len(nodes_g) - 1):
                source = nodes_g.loc[i, 'name']
                target = nodes_g.loc[i + 1, 'name']
                edge_key = tuple(sorted([source, target]))
                if edge_key not in edges_dict:
                    edges_dict[edge_key] = {}
                    edges_dict[edge_key]["flow"] = 1
                    edges_dict[edge_key]["annotations"] = set()
                    edges_dict[edge_key]["genomes"] = [genome]
                else:
                    edges_dict[edge_key]["flow"] += 1
                    if genome not in edges_dict[edge_key]["genomes"]:
                        edges_dict[edge_key]["genomes"].append(genome)

                for a in nodes_g.loc[i, 'annotations']+nodes_g.loc[i + 1, 'annotations']:
                    edges_dict[edge_key]["annotations"].add(a)

        # Specific colors for a selected set of genomes
        colored_genomes = {}
        for g, c in zip(all_genomes, color_genomes):
            if c != '#000000':
                colored_genomes[g] = c
        for (source, target), dic in edges_dict.items():
            link_color = 'gray'
            flow = dic["flow"]
            virtual_flow = flow
            if specifics_genomes is not None and len(specifics_genomes) > 0:
                min_required_shared = max(
                    ceil(min_shared_genome * len(specifics_genomes) / 100), 1)
                max_allowed_extra = ceil(tolerance * len(dic["genomes"]) / 100)
                set_specifics_genomes = set(specifics_genomes)
                set_genomes = set(dic["genomes"])

                intersect = set_specifics_genomes.intersection(set_genomes)
                if len(intersect) >= min_required_shared and len(dic["genomes"]) - len(intersect) <= max_allowed_extra:
                    link_color = color_shared_regions
                    virtual_flow = len(all_genomes)
            label = ""
            if labels:

                first_label = True

                for a in dic["annotations"]:
                    if first_label:
                        label += str(a)
                        first_label = False
                    else:
                        label += ", " + str(a)
            edges.append({
                'data': {
                    'source': source,
                    'target': target,
                    'flow': flow,
                    'genomes': dic["genomes"]
                },
                'style': {
                    'line-color': link_color,
                    'label': label,
                    'text-rotation': 'autorotate',
                    'width': (virtual_flow+int(0.2*len(all_genomes)))/len(all_genomes)*10
                }
            })

            if len(colored_genomes) > 0:
                i = 0
                for g in list(colored_genomes.keys()):
                    if g in dic["genomes"]:
                        sign = 1 if i % 2 == 0 else -1
                        distance = sign * (20 + 10 * (i // 2))
                        edges.append({
                            'data': {
                                'source': source,
                                'target': target,
                                'flow': 1,
                                'genomes': [g]
                            },
                            'classes': f'offset-{i}',
                            'style': {
                                'line-color': colored_genomes[g],
                                'width': 4

                            }
                        })
                        i += 1
        print("nb nodes : " + str(len(nodes)) +
              " - nb edges : " + str(len(edges)))
        return nodes + edges
    else:
        return []


def layout(data=None, initial_size_limit=10):
    all_genomes = get_genomes()
    all_genomes.sort()
    all_chromosomes = get_chromosomes()
    if data != None:
        elements = compute_graph_elements(
            data, all_genomes, initial_size_limit, all_genomes, all_chromosomes, [], [])
    else:
        elements = []
    size_max = 500

    return html.Div([
        dcc.Store(id='zoom_shared_storage_nodes', storage_type='memory'),
        html.Div([
            html.H2("PANORAMA"),

            html.Details([
                html.Summary("ℹ️ Click here to display help"),
                html.P("Panorama is an application that allows you to view and manipulate pan-genomic data. "
                       "This page allows you to view portions of the pan-genome:"),
                html.Ul([
                    html.Li(
                        "Haplotype selection: choose the haplotype that will be used for annotations and coordinates."),
                    html.Li(
                        "Select a chromosome: searches will only be performed on this chromosome."),
                    html.Li("Choose from the following options:"),
                    html.Ul([
                        html.Li("Select the start and end of the region on the selected chromosome. If there is no defined end, then the searched region will be from the defined start to the end of the pangenome. The region should not be too large, otherwise the display will take too long or will not be possible."),
                        html.Li("Search by annotation with a gene name or ID.")
                    ]),
                    html.Li("Select the haplotypes to be viewed: it is possible to exclude some haplotypes. "
                            "In this case, nodes containing only these haplotypes will not be displayed."),
                    html.Li("You can download the current graph by clicking on 'as jpg' or 'as png' button. "
                            "Graphs are saved into ./export/graphs directory.")
                ]),
                html.P("Display settings:"),
                html.Ul([
                    html.Li(
                        "The “minimum node size” slider allows you to hide nodes smaller than this value."),
                    html.Li(
                        "“Hide labels”: allows you to hide labels (they can take up too much space)."),
                    html.Li(
                        "Select layout: allows you to change the representation algorithm if it is not suitable."),
                    html.Li("Search shared paths : if the box is checked, the display is modified to allow the selection of haplotypes for which you want to view shared links. This display can be configured:"),
                    html.Ul([
                        html.Li(
                            "Selection of haplotypes for which common links are sought"),
                        html.Li("Selection of the minimum percentage of selected haplotypes that must be present on the link "
                                "(for example, if the value is set to 50 and 10 haplotypes are selected, at least 5 of the haplotypes must be on the link). "
                                "If the value is zero, then at least one of the selected haplotypes will be required."),
                        html.Li("Tolerance: links containing fewer than [(tolerance / 100) × number of haplotypes passing through this node] "
                                "unselected haplotypes will be reported."),
                        html.Li("Link color: choice of color for reported links")
                    ]),
                    html.Li(
                        "Search shared paths : if the box is unchecked, then it is possible to select a specific color for a each haplotype."),
                    html.Li(
                        "If one of this option is modified, click on the Update graph button."),
                ]),
                html.P("Display description :"),
                html.Ul([
                    html.Li(
                        "By clicking on a node or a link, data of this node will be displayed under the 'update graph' button."),
                    html.Li(
                        "If annotations exists in the visualized region they will be concatenated under the node / link description area."),
                    html.Li("Graph description :"),
                    html.Ul([
                        html.Li(
                            "Node shape : a node is drawn as a circle, unless if it's a repeated node in which case it will be displayed as a square."),
                        html.Li("Node color : The color of the nodes ranges from blue to red. The bluer the color, the less common the node (node associated with only one or a small number of haplotypes). Conversely, red nodes are those of the core genome shared by all haplotypes."),
                        html.Li(
                            "Node size : the size of a node is proportionnal to the size of the sequence associated."),
                        html.Li(
                            "Link size : the size of a link is proportional to the number of haplotypes passing through that link."),
                    ]),
                    html.Li("Zoom :"),
                    html.Ul([
                        html.Li(
                            "To zoom in : first select the nodes you want to zoom in by holding down the left mouse button and Shift key. Then push the 'zoom on selection' button."),
                        html.Li(
                            "To retrieve the initial region just push the 'Reset zoom' button."),
                        html.Li("Important : only the selected nodes will be used to display the new graph, this may hide areas displayed before zooming because nodes are linked by position ascending.")
                    ]),
                ]),
            ], style={"marginBottom": "20px"}),


        ]),
        # Upper block : settings
        html.Div([


            # Left block
            html.Div([
                html.Div([
                    html.Div([
                        html.Label("Reference haplotype", title="Select an haplotype to search / display annotation and to define genomic coordinates to search region.",
                                   style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='genomes-dropdown',
                            options=[{'label': genome, 'value': genome}
                                     for genome in all_genomes],
                            # value=all_genomes[0],
                            clearable=False,
                            style={'width': '300px'}
                        )
                    ], style={'marginRight': '30px'}),
                    html.Div([
                        html.Label("Chromosome", title="The regions / annotations will be related only to this chromosome.",
                                   style={'display': 'block', 'marginBottom': '5px'}),
                        dcc.Dropdown(
                            id='chromosomes-dropdown',
                            options=[{'label': chrom, 'value': chrom}
                                     for chrom in all_chromosomes],
                            clearable=False,
                            placeholder="Choose a chromosome",
                            # value="1",
                            style={'width': '200px'}
                        )
                    ])
                ], style={'display': 'flex', 'padding': '20px', 'border': '1px solid #ccc', 'minWidth': '300px', 'boxSizing': 'border-box'}),


                html.Div([
                    html.H5("Search region", style={'textAlign': 'left', 'marginBottom': '15px'}),
                    html.Div([
                        html.Label("Start : ", title="Start on the selected haplotype / chromosome."),
                        dcc.Input(id='start-input', type='number', style={'width': '100px', 'marginRight': '10px'}),
                        html.Label("End : ", title="End on the selected haplotype / chromosome."),
                        dcc.Input(id='end-input', type='number', style={'width': '100px', 'marginRight': '20px'})
                    ], style={'marginBottom': '10px'}),
                    html.Div([
                        html.Label(
                            "Gene name : ", title="Will be searched on annotation of the selected haplotype / chromosome."),
                        dcc.Input(
                            id='genename-input',
                            type='text',
                            placeholder='Nom de gene',
                            debounce=True,
                            style={'marginRight': "10px"}
                        ),
                        html.Label(
                            "Gene id : ", title="Will be searched on annotation of the selected haplotype / chromosome."),
                        dcc.Input(
                            id='geneid-input',
                            type='text',
                            placeholder='id de gene',
                            debounce=True
                        )
                    ], style={'marginBottom': '20px'}),
                    html.Button('Search', id='search-button',
                                n_clicks=0, style={'marginTop': '10px'}),
                    dcc.Loading(
                        id="loading-search-msg",
                        type="circle",
                        children=html.Div(id="search-message")
                    ),
                    html.Div([
                        html.Label("Haplotypes to visualize :", title="Nodes containing only unselected haplotypes won't be displayed.", style={
                                   'marginBottom': '5px'}),
                        dcc.Checklist(
                            id="genome_selector",
                            options=[{"label": g, "value": g}
                                     for g in all_genomes],
                            value=all_genomes,
                            inline=True
                        )

                    ],
                        style={'marginBottom': '20px'}),
                    # Download graph div
                    html.Div([
                        html.Label(
                            'Download graph:', title="This will download the displayed graph into ./export/graphs directory in jpg or png format.", style={"marginLeft": "10px"}),
                        html.Button("as jpg", id="btn-save-jpg",
                                    style={"marginLeft": "10px"}),
                        html.Button("as png", id="btn-save-png",
                                    style={"marginLeft": "10px"}),
                        # html.Button("as svg", id="btn-save-svg", style={"marginLeft":"10px"}),
                        html.Div(id="dummy-output")
                    ]),
                    # dcc.Loading(id="loading-spinner", type="circle",
                    #             children=html.Div(id="output-zone"))
                ],
                    style={'marginBottom': '20px'}
                )
            ], style={'flex': '1', 'padding': '20px', 'border': '1px solid #ccc'}),

            # Right block
            html.Div([

                html.Div([
                    html.Label("Minimal size of nodes :",
                               title="Nodes with size under this value won't be displayed."),

                    dcc.Slider(
                        id='size_slider',
                        min=0,
                        max=size_max,
                        step=int(100/size_max),
                        value=DEFAULT_SIZE_VALUE,
                        tooltip={"placement": "bottom",
                                 "always_visible": False},
                    ),

                    html.Div(id='size_stats', style={'marginTop': '10px'})

                ]),
                html.Div([
                    html.Div(
                        id='size-output', children='Min node size : 10', style={'margin': '10px'}),
                    html.H4("Layout", title="The layout computes the graph. If the graph is not readable, it is possible to modify the algorithm. Dagre layout will display more linear graphs and fcose will display more compact graphs.", style={
                            'marginLeft': '50px'}),
                    dcc.Dropdown(
                        id='layout-dropdown',
                        options=[
                            {'label': 'fcose', 'value': 'fcose'},
                            {'label': 'dagre', 'value': 'dagre'}
                        ],
                        value='fcose',
                        clearable=False,
                        style={'width': '120px',
                               'display': 'inline-block', 'marginRight': '50px'}
                    ),
                    
                    

                    dcc.Checklist(
                        options=[
                            {'label': 'Hide labels', 
                             'title': 'Uncheck if labels takes too much space on the graph.', 
                             'value': 'hide'
                            }],
                        id='show-labels',
                        style={'marginRight': '30px'}
                    ),
            
                    dcc.Checklist(
                        options=[{
                            'label': 'Show exons',
                            'title': 'Check to display exons.',
                            'value': 'exons'
                        }],
                        id='show-exons',
                        style={'marginRight': '10px'},
                        value=[]
                    ),
            
                    dbc.Input(
                        id='exon-color-picker',
                        type='color',
                        value=DEFAULT_EXONS_COLOR,
                        style={'width': '25px',
                               'height': '25px', 'marginLeft': '10px'}
                    ),

                        

                ], style={'display': 'flex', 'alignItems': 'center', 'gap': '8px'}),
                html.Div(id='nb-noeuds', style={'margin': '10px'}),

                dcc.Checklist(
                    options=[{
                        'label': 'Search shared paths',
                        'title': 'Check to switch to shared path search options.',
                        'value': 'shared'
                    }],
                    id='shared-mode',
                    style={'marginRight': '30px'}
                ),
            
                html.Div([

                    html.Div([
                        html.Label(
                            "Vizualize shared paths :", title="Detect the links through which the selected haplotypes pass - the following input fields allow you to refine the search. ", style={'marginBottom': '20px'}),
                        dcc.Checklist(
                            id="specific-genome_selector",
                            options=[{"label": g, "value": g}
                                     for g in all_genomes],
                            # value=[],
                            inline=True
                        ),
                        html.Label("Min (%) of shared haplotypes : ", title="Min (%) of shared haplotypes = M. Number of selected haplotypes = N. To detect a shared node it must contains almost (M/100) x N of the selected haplotypes. If M = 0 then the minimum number of selected haplotypes will be 1."),
                        dcc.Input(id='min_shared_genomes-input', type='number',
                                  value=0, style={'width': '100px', 'marginRight': '10px'}),
                        html.Label("Tolerance (%) : ", title="Tolerance = T. Number of haplotypes on a node = n. To detect a shared node it must contains less than (T/100) x n of the non selected haplotypes. If T = 0 then detected nodes should contain only selected haplotypes."),
                        dcc.Input(id='tolerance-input', type='number', value=0,
                                  style={'width': '100px', 'marginRight': '20px'}),

                        html.Label(
                            "Link color", title="This color will be used to color links between detected shared nodes."),

                        dbc.Input(
                            id='shared-region-color-picker',
                            type='color',
                            # value=DEFAULT_SHARED_REGION_COLOR,
                            style={'width': '25px',
                                   'height': '25px', 'marginLeft': '10px'}
                        ),

                    ], id='shared-checklist-container'),
                    html.Div([
                        html.Div([
                            html.Label(s),
                            dbc.Input(
                                id={'type': 'color-picker', 'index': i},
                                type='color',
                                value="#000000",
                                style={'width': '25px',
                                       'height': '25px', 'marginLeft': '10px'}
                            )

                        ], style={'marginRight': '10px'}) for i, s in enumerate(all_genomes)
                    ], style={'display': 'flex',
                              'flexWrap': 'wrap',
                              'gap': '5px'}, id='color-picker-container')

                ], id='sample-controls'),
                html.Button("Update graph", id="update-btn",
                            n_clicks=0, style={'marginTop': '10px'}),

                html.Div(html.H4(id='node-info', style={'margin': '10px'})),
                html.Div(html.Label("Annotations :", title="Compiles all annotations for the displayed nodes.", style={
                         'marginBottom': '5px'})),
                html.Div(html.H4(id='annotations-info',
                         style={'margin': '10px'}))
            ], style={'flex': '1', 'padding': '20px', 'border': '1px solid #ccc', 'marginLeft': '20px', 'minWidth': '300px',
                      'boxSizing': 'border-box', 'display': 'flex', 'flex-direction': 'column'})
        ], style={
            'display': 'flex',
            'flex-direction': 'row',
            'flexWrap': 'wrap',
            'justify-content': 'space-between'
        }),
        html.Div([
            # Button on top of the graph
            html.Button(
                "🔍 Zoom on selection",
                id='btn-zoom',
                title='Before using this button, nodes must be selected by holding left mouse button and Shift key.'
            ),
            html.Button(
                "🔄 Reset Zoom",
                id='btn-reset-zoom',
            )
        ]),

        # Graph block
        cyto.Cytoscape(
            id='graph',
            # layout is important to get good visualization result
            # There are many algorithms : cose, cose-bilkent-fcose, euler, dagre, etc.
            # fcose seems to be the most performant
            # dagre is usefull to get a linear representation
            layout={
                'name': 'fcose',
                'maxIterations': 100000,
                'maxSimulationTime': 5000,
                # 'nodeRepulsion': 10000,
                # 'gravity': 0.1,
                # 'gravityRangeCompound': 1.5,
                # 'idealEdgeLength': 100,
                # 'componentSpacing': 100,
                # 'nodeDimensionsIncludeLabels': True,
                # 'edgeElasticity': 0.1,
                # 'nestingFactor': 0.8,
                # 'tile': True,
                'quality': "proof",
                'fit': True
            },
            style={'width': '100%', 'height': '1000px'},
            elements=elements,
            # minZoom=0.1,
            # maxZoom=5,
            zoomingEnabled=True,
            userZoomingEnabled=True,
            userPanningEnabled=True,
            wheelSensitivity=0.1,
            #responsive=True,
            autoRefreshLayout=True,
            boxSelectionEnabled=True,
            autoungrabify=False,
            stylesheet=compute_stylesheet(0),

        )
        ])


# Callback to get nodes or link info when clicking on it
@app.callback(
    Output('node-info', 'children'),
    Input('graph', 'tapNodeData'),
    Input('graph', 'tapEdgeData')
)
def display_element_data(node_data, edge_data):
    triggered_id = ctx.triggered_id

    if triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapEdgeData' and edge_data:
        return (
            f"Selected link : {edge_data.get('source')} → {edge_data.get('target')}\n"
            f"• Flow : {edge_data.get('flow')}\n"
            f"• Haplotypes : {', '.join(edge_data.get('genomes', []))}"
        )
    elif triggered_id == 'graph' and ctx.triggered[0]['prop_id'] == 'graph.tapNodeData' and node_data:
        return (
            f"Selected node : {node_data.get('label', node_data.get('name'))}\n"
            f"• Size : {node_data.get('size')}\n"
            f"• Flow : {node_data.get('flow')}\n"
            f"• Ref node : {node_data.get('ref_node')}\n"
            f"• Haplotypes : {', '.join(node_data.get('genomes', []))}"
            f"• Annotations : {', '.join(node_data.get('annotations', []))}"
            f"• Features : {', '.join(node_data.get('features', []))}"
        )
    return "Click on a node or link to display data."

# Main callback to update graph when changing size, or selecting genomes, etc.


@app.callback(
    Output("graph", "elements"),
    Output("nb-noeuds", 'children'),
    Output('shared_storage_nodes', 'data', allow_duplicate=True),
    Output('search-message', 'children'),
    Output('annotations-info', 'children'),
    Output('graph', 'stylesheet'),
    Output('home-page-store', 'data', allow_duplicate=True),
    Output('graph', 'selectedNodeData'),
    Output('graph', 'selectedEdgeData'),
    Output('zoom_shared_storage_nodes', 'data', allow_duplicate=True),
    State('genome_selector', 'value'),
    State('shared-mode', 'value'),
    State('specific-genome_selector', 'value'),
    State({'type': 'color-picker', 'index': ALL}, 'value'),
    Input('show-labels', 'value'),
    Input('update-btn', 'n_clicks'),
    Input('btn-zoom', 'n_clicks'),
    Input('btn-reset-zoom', 'n_clicks'),
    State('graph', 'selectedNodeData'),
    # Input('size_slider', 'value'),
    State('home-page-store', 'data'),
    Input('search-button', 'n_clicks'),
    State('start-input', 'value'),
    State('end-input', 'value'),
    State('genename-input', 'value'),
    State('geneid-input', 'value'),
    State('genomes-dropdown', 'value'),
    State('chromosomes-dropdown', 'value'),
    State('shared_storage', 'data'),
    State('shared_storage_nodes', 'data'),
    State('min_shared_genomes-input', 'value'),
    State('tolerance-input', 'value'),
    State('shared-region-color-picker', 'value'),
    State('zoom_shared_storage_nodes', 'data'),
    State('show-exons', 'value'),
    State('exon-color-picker', 'value'),
    prevent_initial_call=True
)
def update_graph(selected_genomes, shared_mode, specifics_genomes, color_genomes, show_labels, update_n_clicks, zoom_clicks, reset_zoom_bouton_clicks, selected_nodes_data, home_data_storage, n_clicks, start, end, gene_name, gene_id, genome, chromosome, data_storage, data_storage_nodes, min_shared_genome, tolerance, shared_regions_link_color, zoom_shared_storage, show_exons, exons_color):
    ctx = dash.callback_context
    message = ""
    triggered_id = ctx.triggered_id
    if home_data_storage is None:
        home_data_storage = {}
    if home_data_storage is not None and 'slider_value' in home_data_storage:
        size_slider_val = home_data_storage['slider_value']
    else:
        size_slider_val = DEFAULT_SIZE_VALUE
    # save the parameters into store
    if genome is not None:
        home_data_storage["selected_genome"] = genome
    if chromosome is not None:
        home_data_storage["selected_chromosome"] = chromosome
    if shared_regions_link_color is not None:
        home_data_storage["shared_regions_link_color"] = shared_regions_link_color
    if start is not None:
        home_data_storage["start"] = start
    if end is not None:
        home_data_storage["end"] = end
    if gene_name is not None:
        home_data_storage["gene_name"] = gene_name
    if gene_id is not None:
        home_data_storage["gene_id"] = gene_id
    if min_shared_genome is None:
        min_shared_genome = 100
    if tolerance is None:
        tolerance = 0
    if color_genomes is not None:
        home_data_storage["color_genomes"] = color_genomes
    if specifics_genomes is not None:
        home_data_storage["specifics_genomes"] = specifics_genomes
    # zoom on selected nodes
    zoom_shared_storage_out = zoom_shared_storage or {}
    if triggered_id == "btn-zoom":
        if selected_nodes_data is not None and len(selected_nodes_data) > 0:
            data_to_plot = {}
            selected_ids = set([node['name']
                               for node in selected_nodes_data if 'name' in node])
            data_to_plot = {
                k: v for k, v in data_storage_nodes.items() if k in selected_ids}
            zoom_shared_storage_out = data_to_plot
            print("Zoom : datat to plot length : " + str(len(data_to_plot)))
    if triggered_id == "btn-reset-zoom":
        zoom_shared_storage_out = {}
    if len(zoom_shared_storage_out) > 0:
        data_to_plot = zoom_shared_storage_out
    else:
        data_to_plot = data_storage_nodes
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
    exons=False
    if show_exons and 'exons' in show_exons:
        exons = True
        
    all_genomes = data_storage["genomes"]
    all_chromosomes = data_storage["chromosomes"]
    if ctx.triggered and len(ctx.triggered[0]['prop_id'].split('.')) > 0:
        input_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if input_id == "search-button" and n_clicks > 0:
            if start is not None:
                new_data = get_nodes_by_region(
                    genome, chromosome=chromosome, start=start, end=end)
                data_storage_nodes = new_data
                print("len new_data : " + str(len(new_data)))
            else:
                if ((gene_name is not None and gene_name != "") or (gene_id is not None and gene_id != "")) and chromosome is not None:
                    if gene_name is not None and gene_name != "":
                        new_data = get_nodes_by_gene(
                            genome, chromosome=chromosome, gene_name=gene_name)
                    else:
                        new_data = get_nodes_by_gene(
                            genome, chromosome=chromosome, gene_id=gene_id)
                    data_storage_nodes = new_data
                else:
                    new_data = get_nodes_by_region(
                        genome, chromosome=chromosome, start=0, end=end)
            elements = compute_graph_elements(new_data, selected_genomes, size_slider_val, all_genomes, all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome, 
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color, exons=exons, exons_color=exons_color)
            zoom_shared_storage_out = {}
            if len(elements) == 0:
                message=html.Div("❌ No data found.", style=warning_style)
        else:
            print(f"min node size : {size_slider_val}")
            elements = compute_graph_elements(data_to_plot, selected_genomes, size_slider_val, all_genomes, all_chromosomes, specifics_genomes_list,
                                              color_genomes_list, labels=labels, min_shared_genome=min_shared_genome, 
                                              tolerance=tolerance, color_shared_regions=shared_regions_link_color, exons=exons, exons_color=exons_color)
    else:
        print(f"min node size : {size_slider_val}")
        elements = compute_graph_elements(data_to_plot, selected_genomes, size_slider_val, all_genomes, all_chromosomes, specifics_genomes_list,
                                          color_genomes_list, labels=labels, min_shared_genome=min_shared_genome, 
                                          tolerance=tolerance, color_shared_regions=shared_regions_link_color, exons=exons, exons_color=exons_color)

    defined_color = 0
    if color_genomes is not None:
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

    

    return elements, f"{count} displayed nodes", data_storage_nodes, message, annotations, stylesheet, home_data_storage, [], [], zoom_shared_storage_out


# color picker
@app.callback(
    Output('color-picker-container', 'style'),
    Output('shared-checklist-container', 'style'),
    Input('shared-mode', 'value')
)
def toggle_inputs(shared_mode):
    if shared_mode and 'shared' in shared_mode:
        return {'display': 'none'}, {'display': 'block'}
    else:
        return {'display': 'flex', 'flexWrap': 'wrap'}, {'display': 'none'}


@app.callback(
    Output('home-page-store', 'data', allow_duplicate=True),
    Output("size-output", 'children'),
    Input('size_slider', 'value'),
    State('home-page-store', 'data'),
    prevent_initial_call=True
)
def save_slider_value(size_slider_val, data):
    if data is None:
        data = {}
    data['slider_value'] = size_slider_val
    return data, f"Min node size  : {size_slider_val}"

# Restore value after navigation


@app.callback(
    Output('size_slider', 'value'),
    Output('chromosomes-dropdown', 'value'),
    Output('genomes-dropdown', 'value'),
    Output('start-input', 'value'),
    Output('end-input', 'value'),
    Output('shared-region-color-picker', 'value'),
    Output('specific-genome_selector', 'value'),
    Input('url', 'pathname'),
    Input('home-page-store', 'data'),
    State('genomes-dropdown', 'options'),
    State('chromosomes-dropdown', 'options'),
    State('specific-genome_selector', 'value'),
)
def update_parameters_on_page_load(pathname, data, options_genomes, options_chromosomes, specifics_genomes):
    slider_value = DEFAULT_SIZE_VALUE
    selected_genome = None
    selected_chromosome = None
    start_input = None
    end_input = None
    shared_regions_link_color = DEFAULT_SHARED_REGION_COLOR
    if specifics_genomes is not None:
        selected_shared_genomes = specifics_genomes
    else:
        selected_shared_genomes = []
    if data is None:
        data = {}
    if "slider_value" in data:
        slider_value = data["slider_value"]
    if "selected_genome" in data:
        selected_genome = data["selected_genome"]
    else:
        if options_genomes:
            selected_genome = options_genomes[0]["value"]
    if "selected_chromosome" in data:
        selected_chromosome = data["selected_chromosome"]
    else:
        if options_chromosomes:
            selected_chromosome = options_chromosomes[0]["value"]
    if "start" in data:
        start_input = data["start"]
    if "end" in data:
        end_input = data["end"]
    if "shared_regions_link_color" in data:
        shared_regions_link_color = data["shared_regions_link_color"]
    if "specifics_genomes" in data:
        selected_shared_genomes = data["specifics_genomes"]

    return slider_value, selected_chromosome, selected_genome, start_input, end_input, shared_regions_link_color, selected_shared_genomes


# Algorithm cytoscape choice
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
            'maxIterations': 100000,
            'maxSimulationTime': 5000,
            # 'nodeRepulsion': 10000,
            # 'gravity': 0.1,
            # 'gravityRangeCompound': 1.5,
            # 'idealEdgeLength': 100,
            # 'componentSpacing': 100,
            # 'nodeDimensionsIncludeLabels': True,
            # 'edgeElasticity': 0.1,
            # 'nestingFactor': 0.8,
            # 'tile': True,
            'quality': "proof",
            'fit': True
        }


######## download graph callbacks ###############
@app.callback(
    Output('dummy-output', 'children'),
    Input('graph', 'imageData'),
    State('chromosomes-dropdown', 'value'),
    State('genomes-dropdown', 'value'),
    State('start-input', 'value'),
    State('end-input', 'value')
)
def save_image_to_file(image_data, chromosome, genome, start, end):
    if not image_data:
        raise PreventUpdate
    # S'assurer que le dossier existe
    os.makedirs(EXPORT_DIR, exist_ok=True)

    # imageData est typiquement sous la forme : 'data:image/jpeg;base64,...'
    header, base64_data = image_data.split(',', 1)

    # Get image formaty from header
    if 'image/png' in header:
        ext = 'png'
    elif 'image/jpeg' in header:
        ext = 'jpg'
    elif 'image/svg+xml' in header:
        ext = 'svg'
    else:
        raise ValueError("Unrecognized format in imageData")

    file_name = "graph_"+str(genome)+"_chr_"+str(chromosome) + \
        "_start_"+str(start)+"_end_"+str(end)+"."+ext
    save_path = os.path.join(os.getcwd(), EXPORT_DIR, file_name)
    # Décode les données base64 en bytes
    image_bytes = base64.b64decode(base64_data)

    # Enregistre dans un fichier
    with open(save_path, 'wb') as f:
        f.write(image_bytes)

    return f"Image downloaded in {save_path}"


@app.callback(
    Output("graph", "generateImage"),
    Input("btn-save-jpg", "n_clicks"),
    Input("btn-save-png", "n_clicks"),
    # Input("btn-save-svg", "n_clicks")
)
def trigger_image_save(n_clicks_jpg, n_clicks_png):
    if not ctx.triggered_id:
        raise PreventUpdate

    # Get image format
    fmt = ctx.triggered_id.split('-')[-1]
    if fmt == 'svg':
        return {'type': fmt, 'action': 'download'}
    else:
        return {'type': fmt, 'action': 'store'}


######## zoom on selection callbacks ###############

# @app.callback(
#     Output('selected-nodes', 'children'),
#     Input('graph', 'selectedNodeData')
# )
# def display_selected_nodes(data):
#     if not data:
#         return "No nodes selected."

#     # Affiche les infos utiles des nœuds sélectionnés
#     return "\n".join(f"- {node.get('label')} (id: {node.get('id')})" for node in data)
