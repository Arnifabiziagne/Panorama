Panorama_graph project : tools to manipulate and visualize pangenomes graphs.

# Overview
Panorama_graph is based on a graph database modelisation (neo4j) from a gfa graph.
It allows the following functionalities : 
- Load a GFA file into a neo4j database (including optionnal annotations if present)
- Compute shared regions between a set of selected haplotypes
- Compute a phylogenetic tree from a selected region (neighbour joining with a distance matrix based on Jaccard index)
- Visualize a region and annotation of the pangenome

## Install neo4
    See the readme in neo4j_install directory

## Create environment
    TODO

## Use the tool
    You can use the tool from the graphical interface. To launch it, just launch the index.py file.