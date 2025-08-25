Panorama_graph project : tools to manipulate and visualize pangenomes graphs.

# Overview
Panorama_graph is based on a graph database modelisation (neo4j) from a gfa graph.
It allows the following functionalities : 
- Load a GFA file into a neo4j database (including optionnal annotations if present)
- Compute shared regions between a set of selected haplotypes
- Compute a phylogenetic tree from a selected region (neighbour joining with a distance matrix based on Jaccard index)
- Visualize a region and annotation of the pangenome

## Prerequisites
- Docker available
- 16 GB RAM (32 Go+ recommanded for big data)
- sufficient hard disk space (approximately 10 times the size of the GFA, ideally SSD)
- Conda
- The version of panorama_graph used to construct the database must be compatible with the version to visualize and analyse data (the version of the database is indicated in the Stats node)

## Quickstart
- Installation (only the first time) :
  - Copy all the files / directories in the Panorama_graph project where do you want to store your data
  - If there is a dump file available (file named neo4j.dump) or csv files to import (nodes.csv, sequences.csv and relations.csv), move or copy it into ./import directory
  - Go into neo4j_install directory and run the script (replace $container_name with the desired name) : bash ./setup_neo4j.sh --container-name $container_name
  - Create the conda environnement : conda env create -f panorama_graph.yaml 
  - Make launcher executable : chmod +x launcher.sh
- Launch the tool : 
  - Easy way : just execute launcher : ./launcher.sh and go to http://localhost:8050
  - Manual way : 
      - Load conda environnement : conda activate panorama_graph
      - Run index.py (python index.py) and go to http://localhost:8050
  Prepare database
  - Go to "DB management" page, give a name to your container (e.g. DB_my_species_PGGB) and click on "Create new DB"
  - If no data (dump or csv) where in the ./data/import directory, load GFA data by selecting GFA file (if the file concerns only a single chromosome it is required to set the chromosome name)
  - Load annotations by selecting the file and the genome related. Before to load annotations it is required that indexex are fully created (after creating data or loading GFA the indexes are automatically created but if data are big it requires some time)
  Once data are loaded the tool can be ued (see quick pages description).

## Quick pages description
  The menu allow to navigate on differents pages :

  - DB management : this page is used only on the start to create DB and load data.
  - Home page : page to vizualise data (by defining the chromosome, start and stop or genome and gene_name / gene_id).
  - GWAS : page to detect the nodes shared by a selection of haplotypes. It computes the list of identified regions that can be exported in csv. Sequence associated to a region can be seen by clicking on "size" column.
  - Phylogenetic : on left it is possible to load a reference phylogenetic tree. On right, by clicking on the "Plot tree..." button it computes the tree of the region defined in the Home page.
  - Sequences : by clicking on the button it computes the sequence for each haplotype of the region selected in the home page.

## Install neo4
    See the readme in neo4j_install directory

## Create environment
- Create conda env : conda env create -f panorama_graph.yaml
- Load conda env : conda activate panorama_graph

## Generate the database
    There are 3 ways to generate database :
    - From a dump file : this is the fastest way but the dump must be available. It uses the neo4j-admin load functionnality.
    - From csv file : it is a fast way to create rthe database if the csv files are available. It uses the neo4j-admin import functionnality. The difference with the dump is that annotations and indexes won't be created.
    - From the GFA file (and gtf / gff if available) : use the construct_DB function of neoj4_DB_contruction.py script or the DB management page of IHM. 
        If the pangenome is big this can take a long time, in this case it is recommanded to use the load_gfa_data_to_csv function to generate csv files. According to the memory available, it is necessary to limit the 
        batch size. If a gtf / gff file is present, the genome_ref must be set in order to link annotations nodes with the main nodes of pangenome.

## Use the tool 
    You can use the tool from the graphical interface. To launch it, just launch the index.py file : python index.py