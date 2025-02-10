Projet panorama

# Installation 

## Install miniconda
    panorama use miniconda environment (https://docs.anaconda.com/miniconda/miniconda-install/)

## Create environment (local)
    conda env create -f panorama.yaml 

## Create environment (genobio)
    conda env create -f panorama_cluster.yaml

# Usage 
    conda activate panorama

    module load bioinfo/RAxML/8.2.12 (only on genobio cluster) 


    python pangenome_analysis.py -f pangenome_gfa_filename -p project_name -n 100000 -w True -m "random" -r True -s True -c color_filename

# Options
    -f / --graphfilename : Specify the gfa file location.
    -p / --project_name : Specify the project name : a directory with this name will be created and the out files will be generated in this directory.
    -n / --nodes_number : Number of nodes to select. These nodes will be used to compute the morphologic matrix (presence or absence of the selected nodes). This morphologic matrix will be used by raxml. Default = 10000.
    -w / --weigthed : True => Jaccard distance will be computed with the weight (= size of the nodes), False => Jaccard distance will not be weighted. Default = True.
    -m / --method	: "random" => the nodes will be chosen randomly, other value => the nodes will be chosen by length (the bigger nodes will be chosen). Default = "random".
    -r / --redundancy : True => the redundant nodes will be taken into account. False => the redundant nodes will only be taken into account once. Default = True.
    -s / --strand : True => the strand will be taken into account. False => the strand is ignored. Default = True.
    -c / --colorfilename : file containing information to color the tree. The file must contain a column “sample” = sample name, “category” = category abbreviation which will be added to the sample name, “color” = color code in #aaaaaa format.
    -h / --help : Print this Help.