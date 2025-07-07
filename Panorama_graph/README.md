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

## Generate the database
    There are 2 ways to generate database :
    - From a dump file : this is the fastest way but the dump must be available.
    - From the GFA file (and gtf / gff if available) : use the construct_DB function of neoj4_DB_contruction.py script. 
        If the pangenome is big this can take a long time. According to the memory available, it is necessary to limit the 
        batch size. If a gtf / gff file is present, the genome_ref must be set in order to link annotations nodes with the main nodes of pangenome.

## Use the tool
    You can use the tool from the graphical interface. To launch it, just launch the index.py file.