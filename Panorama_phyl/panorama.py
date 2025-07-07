#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:57:26 2025

@author: fgraziani
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import matplotlib
import random as rand
import os
from tqdm import tqdm
import subprocess
import re
import time
from Bio.Seq import Seq

"""
Function to load a GFA file
The function will extract the nodes and paths, then return for each genome:
The list of traversed nodes and their size.
The function also allows for selecting nodes either randomly
or based on a size criterion and indicates for each genome whether it contains
these nodes or not.

Parameters
    @file_name: name of the GFA file to load
    @raxml_sample_percent : specifies the desired number of nodes in percent of total nodes
        This value is set to 0.1% by default
    @strand: True to use strand, False otherwise
    @chromosome_file : if the GFA concerns a single chromosome, specify the chromosome number (or X / Y)
    @masked_nodes_file_names : if some nodes must be masked, specify the nodes id in a text file (each line contains on id)
    @private_haplotypes_filter : the nodes shared less than this value won't be considered
    

Returns
    genome_dic: dictionary
        Dictionary with genome names as keys and arrays as values. Each array has a size equal to the number of nodes and contains the size of each node.
    genome_redondant_dic: dictionary
        Dictionary with genome names as keys and dictionaries as values. Each sub-dictionary has node names as keys and arrays as values, where each array represents the number of occurrences of each node.
        If strand is used, then the array has a size of 2 x number_of_nodes (one for - nodes and one for + nodes).
    pav_dic: dictionary
        Dictionary with genome names as keys and, for each genome, an array indicating presence (1) or absence (0) of selected nodes.
    stats: dictionary
        Dictionary containing genome names as keys and, for each genome, a dictionary of chromosomes with the number of forward and reverse strands.
        This object can then be exported using the export_stats function.
"""


def load_gfa(file_name, raxml_sample_percent=0.1, strand=True,
                        chromosome_file=None, masked_nodes_file_name = None, private_haplotypes_filter = 0):

    nodes_nb = 0
    nb_liens = 0
    nb_chemins = 0
    stats = {}
    minus_nb = 0
    plus_nb = 0
    set_masked_nodes = set()
    filtered_nodes = set()
    set_genomes = set()

    # Getting the masked nodes if masked_nodes_file_names is not empty
    if masked_nodes_file_name is not None and masked_nodes_file_name != "":
        file_masked_node = open(masked_nodes_file_name, "r")
        line = file_masked_node.readline()
        while line:
            set_masked_nodes.add(line.strip())
            line = file_masked_node.readline()
    print("Number of masked nodes : " + str(len(set_masked_nodes)))

    # Definition of the separators according to the line's type (P ou W)
    # sep[0] for P lines
    # sep[1] for W lines
    sep = ["[,;.*]", "(<|>)"]
    dic_count_direct_reverse_strand = {}
    file = open(file_name, "r")
    with file:
        line = file.readline()
        print("File opening " + str(file_name))
        genome_dic = {}
        genome_redondant_dic = {}
        node_lentgh_dic = {}
        node_index = {}
        pav_dic = {}
        index_noeuds_raxml_dic = {}
        nodes_length_list = []
        nodes = {}
        # First GFA browsing to find nodes to create the data frame and to get the nodes size
        tps1 = time.time()
        while line:
            ligne_dec = line.split()
            if ligne_dec[0] == 'S' and ligne_dec[1] not in node_lentgh_dic :
                    nodes[ligne_dec[1]] = nodes_nb
                    nodes_nb += 1
                    node_lentgh_dic[ligne_dec[1]] = len(ligne_dec[2])
                    nodes_length_list.append([ligne_dec[1], len(ligne_dec[2])])

            if ligne_dec[0] == 'L':
                nb_liens += 1

            if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                nb_chemins += 1

                if ligne_dec[0] == 'P':
                    ind = 2
                else:
                    ind = 6
                chromosome, genome = get_chromosome_genome(line, haplotype=False, chromosome_file=chromosome_file)
                set_genomes.add(genome)
                # the dic_count_direct_reverse_strand dictionary is used to detect haplotypes
                # for which there are more reverse nodes than direct and which should probably be reversed.
                if genome not in dic_count_direct_reverse_strand:
                    dic_count_direct_reverse_strand[genome] = {}

                if chromosome not in dic_count_direct_reverse_strand[genome]:
                    dic_count_direct_reverse_strand[genome][chromosome] = {"+": 0, "-": 0}

                if ligne_dec[0] == 'P':
                    dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count("+")
                    dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("-")
                else:
                    dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count(">")
                    dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("<")

            line = file.readline()

        dic_to_check = {}
        raxml_sample_size = min(nodes_nb,int(raxml_sample_percent * nodes_nb /100))
        print("Number of sampled nodes : " + str(raxml_sample_size))
        # Search of genomes / chromosomes potentialy inversed
        for g in dic_count_direct_reverse_strand:
            for c in dic_count_direct_reverse_strand[g]:
                if dic_count_direct_reverse_strand[g][c]["-"] > dic_count_direct_reverse_strand[g][c]["+"]:
                    if g in dic_to_check:
                        dic_to_check[g].append(c)
                    else:
                        dic_to_check[g] = [c]

        print("Inversed Genomes / chromosomes detected : " + str(dic_to_check))
        print("Number of paths : " + str(nb_chemins))
        dic_count_direct_reverse_strand = {}
        # Sorting by node length
        nodes_length_list.sort(key=lambda x: x[1])
        
        if private_haplotypes_filter > 0 :
            gi = 0
            genomes_index = {}
            for g in set_genomes:
                genomes_index[g] = gi
                gi += 1
            filtered_nodes = mask_nodes(file_name, nodes, genomes_index, private_haplotypes_filter,chromosome_file)
            set_masked_nodes = set_masked_nodes | filtered_nodes
        print("Number of filtered nodes : " + str(len(filtered_nodes)))
        print("Number of total masked nodes : " + str(len(set_masked_nodes)))
        nodes_nb = 0
        for n in nodes:
            if n not in set_masked_nodes:
                node_index[n] = nodes_nb
                nodes_nb += 1
        nodes = {}

        raxml_nodes_list = rand.sample(list(node_index.keys()), raxml_sample_size)

        nodes_length_list = []    
        raxml_nodes_nb = len(raxml_nodes_list)
        i = 0
        for node in raxml_nodes_list:
            index_noeuds_raxml_dic[node] = i
            i += 1


        print("Number of selected nodes to compute tree : " + str(raxml_nodes_nb))
        tps2 = time.time()
        print("Time parsing nodes : " + str(tps2 - tps1))
        # 2nd GFA browsing to get paths and construct the dictionnary returned by the function
        file.seek(0, 0)
        line = file.readline()
        num_chemin = 0
        nb_noeuds_chemin_max = 0
        with tqdm(total=nb_chemins) as bar:
            while line:
                ligne_dec = line.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    num_chemin += 1

                    # split path according to P lines or W lines separators
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0

                    else:
                        ind = 6
                        walk = 1
                    chromosome, genome = get_chromosome_genome(line, haplotype=False, chromosome_file=chromosome_file)
                    print("chromosome : " +str(chromosome) +  " genome : " + str(genome))
                    if genome not in genome_dic:
                        if strand:
                            genome_redondant_dic[genome] = np.zeros(2 * nodes_nb, dtype=np.int32)
                            genome_dic[genome] = np.zeros(2 * nodes_nb, dtype=np.int32)
                            
                            pav_dic[genome] = np.zeros(2 * raxml_nodes_nb, dtype=int)
                        else:
                            genome_redondant_dic[genome] = np.zeros(nodes_nb, dtype=np.int32)
                            genome_dic[genome] = np.zeros(nodes_nb, dtype=np.int32)
                            
                            pav_dic[genome] = np.zeros(raxml_nodes_nb, dtype=int)

                    if genome not in stats:
                        stats[genome] = {"total+": 0, "total-": 0}
                    if chromosome not in stats[genome]:
                        stats[genome][chromosome] = {"+": 0, "-": 0}

                    nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                    strand_list = []
                    nodes_list = []
                    i = 0
                    while i < len(nodes_list_) and len(nodes_list_[i]) > 0:
                        if walk == 1 and nodes_list_[i] in ["<", ">"]:
                            strand_list.append(nodes_list_[i])
                            nodes_list.append(nodes_list_[i + 1])
                            i += 2
                        else:
                            if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                                strand_list.append(nodes_list_[i][-1])
                                nodes_list.append(nodes_list_[i][:-1])
                                i += 1
                            else:
                                i += 1
                    for n in range(0, len(nodes_list)):
                        node = nodes_list[n]
                        s = strand_list[n]

                        if node != ""  and node in node_index:
                            if strand:
                                if s in ["<", "-"]:
                                    genome_redondant_dic[genome][node_index[node]] += 1
                                    genome_dic[genome][node_index[node]] = node_lentgh_dic[node]
                                    stats[genome][chromosome]["-"] += 1
                                    stats[genome]["total-"] += 1
                                    minus_nb += 1
                                else:
                                    genome_redondant_dic[genome][node_index[node] + nodes_nb] += 1
                                    genome_dic[genome][node_index[node] + nodes_nb] = node_lentgh_dic[node]
                                    stats[genome][chromosome]["+"] += 1
                                    stats[genome]["total+"] += 1
                                    plus_nb += 1
                            else:
                                genome_redondant_dic[genome][node_index[node]] += 1
                                genome_dic[genome][node_index[node]] = node_lentgh_dic[node]
                                if s in ["<", "-"]:
                                    stats[genome][chromosome]["-"] += 1
                                    minus_nb += 1
                                else:
                                    stats[genome][chromosome]["+"] += 1
                                    plus_nb += 1
                            if node in index_noeuds_raxml_dic:
                                if strand:
                                    if s in ["<", "-"]:
                                        pav_dic[genome][index_noeuds_raxml_dic[node]] = int(1)
                                    else:
                                        pav_dic[genome][index_noeuds_raxml_dic[node] + raxml_nodes_nb] = int(1)
                                else:
                                    pav_dic[genome][index_noeuds_raxml_dic[node]] = int(1)

                    if nb_noeuds_chemin_max < len(nodes_list):
                        nb_noeuds_chemin_max = len(nodes_list)
                    bar.update(1)
                line = file.readline()

    file.close()
    print("Number of nodes : " + str(nodes_nb) + "\nLinks : " + str(nb_liens) + "\nPaths : " + str(
        nb_chemins) + "\nMaximum number of nodes in a path : " + str(nb_noeuds_chemin_max))
    tps3 = time.time()
    print("Time parsing paths : " + str(tps3 - tps2) + "Total time : " + str(tps3 - tps1))
    print("Strand's statistics : reverse : " + str(minus_nb) + " direct : " + str(plus_nb))
    return genome_dic, genome_redondant_dic, pav_dic, stats


"""
Function that returns a list of nodes to mask
The nodes are selected according to this rule :
    if a node if present in < m haplotype or in > n-m haplotype (n is the total number of haplotypes) then it is masked
Parameters
    @file_name : name of output file
    @private_haplotypes_filter : min number of haplotypes to take a node into account
Returns 
"""
def mask_nodes(file_name, nodes_index, genomes_index, private_haplotypes_filter = 0, chromosome_file=None):
    set_masked_nodes = set()
    
    if private_haplotypes_filter == 0 :
        return set_masked_nodes
    else:
        # Definition of the separators according to the line's type (P ou W)
        # sep[0] for P lines
        # sep[1] for W lines
        sep = ["[,;.*]", "(<|>)"]
        file = open(file_name, "r")
        matrix = np.zeros((len(nodes_index), len(genomes_index)), dtype=bool)
        with file:
            line = file.readline()
            print("File opening " + str(file_name))
            dic_nodes = {}
            # First GFA browsing to find nodes to create the data frame and to get the nodes size
            tps1 = time.time()
            while line:
                ligne_dec = line.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':

                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0

                    else:
                        ind = 6
                        walk = 1
                    chromosome, genome = get_chromosome_genome(line,haplotype=False,chromosome_file=chromosome_file)
                    print("chromosome : " +str(chromosome) +  " genome : " + str(genome))
                    nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                    strand_list = []
                    nodes_list = []
                    i = 0
                    while i < len(nodes_list_) and len(nodes_list_[i]) > 0:
                        if walk == 1 and nodes_list_[i] in ["<", ">"]:
                            strand_list.append(nodes_list_[i])
                            nodes_list.append(nodes_list_[i + 1])
                            i += 2
                        else:
                            if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                                strand_list.append(nodes_list_[i][-1])
                                nodes_list.append(nodes_list_[i][:-1])
                                i += 1
                            else:
                                i += 1
                    for n in range(0, len(nodes_list)):
                        node = nodes_list[n]
                        i = nodes_index[node]-1
                        j = genomes_index[genome]
                        matrix[i,j] = True
                line = file.readline()
            if private_haplotypes_filter < len(genomes_index) :
                for n in nodes_index:
                    if sum(matrix[nodes_index[n]-1,:]) <= private_haplotypes_filter:
                        set_masked_nodes.add(n)
        matrix = None
        print("nodes masked in " + str(time.time()-tps1))
        return set_masked_nodes
                    

"""
Function to save stats in csv file
Parameters
    @file_name : name of output file
    @stats : dictionnary from fonction load_gfa
Returns 
"""
def export_stats(file_name, stats):
    file = open(file_name, "w")
    with file:
        file.write("genome,chromosome,direct strand number, reverse strand number\n")
        for g in stats:
            for chr in stats[g]:
                if chr != "total+" and chr != "total-":
                    file.write(
                        g + "," + str(chr) + "," + str(stats[g][chr]["+"]) + "," + str(stats[g][chr]["-"]) + "\n")
    file.close()


"""
Function to find the genome and chromosome of a P or W line
The P line is supposed to match one of the following pattern :
    - genome#haplotype#chromosome
    - genome.haplotype.chromosome
The W line (to be preferred) : genome = 2nd element of the line +"_" + 3rd element of the line, chromosome = 4th element of the line
If the file concerned only one chromosome with multiple haplotype, chromosome_file parameter should contain the chromosome reference (string)
Parameters :
    - WP_line : W or P line
    - haplotype : is haplotype is true the genome will be named : genome_haplotype, else only genome
    - chromosome_file : if set, the chromosome value is fixed to this value
"""


def get_chromosome_genome(WP_line, haplotype=True, chromosome_file=None):
    ligne_dec = WP_line.split()
    if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
        chromosome = "0"
        if ligne_dec[0] == 'P':
            if (len(ligne_dec[1].split("#")) > 1):
                name_dec = ligne_dec[1].split("#")
            else:
                name_dec = ligne_dec[1].split(".")
            if haplotype and len(name_dec) > 0:
                genome = name_dec[0] + "_" + name_dec[1]
            else:
                genome = name_dec[0]
            if len(name_dec) > 0:
                chromosome = re.sub("^0*", "", name_dec[-1].upper().replace("CHR", ""))
        else:
            chromosome = str(ligne_dec[3])
            if haplotype:
                genome = ligne_dec[1] + "_" + ligne_dec[2]
            else:
                genome = ligne_dec[1]
    if chromosome_file != None and chromosome_file != "":
        chromosome = chromosome_file
    return chromosome, genome


"""
This function computes the Jaccard distance from the structure returned by the GFA

Parameters
    @genome_dic: the structure returned by the GFA loading function 
        (dictionary containing genomes, the nodes present in each genome, and their size)
    @genome_redondant_dic: the structure returned by the GFA loading function 
        (dictionary containing genomes, the nodes present in each genome, and their number of occurrences)
    @project_directory: directory where output files will be stored
    @redundancy: if True, Jaccard distance is calculated on all nodes including redundant ones; 
        otherwise, only on presence/absence of nodes
    @color_filename: this file allows coloring of the leaves of the generated tree
        The file must be in CSV format, separated by "," and must contain a header row
        Headers must include at least:
            - sample: sample name, which must match genome names in the pangenome
            - category: category name to be added to the leaf labels
            - color: color code associated with the sample

Returns
    - dfP: weighted distance matrix
    - dfNP: unweighted distance matrix
"""


def calculate_distance_jaccard(genome_dic, genome_redondant_dic, project_directory, redundancy=True,
                               color_file_name=None):
    print("Computing the distances")
    genome_list = list(genome_dic.keys())
    dfP = pd.DataFrame(np.ones(shape=(len(genome_list), len(genome_list))), columns=genome_list, index=genome_list)
    dfNP = pd.DataFrame(np.ones(shape=(len(genome_list), len(genome_list))), columns=genome_list, index=genome_list)
    print("Genomes number : " + str(len(genome_list)))
    for i in range(0, len(genome_list)):
        print("Genome : " + str(i) + " name : " + genome_list[i])
        j = i
        genome = genome_list[i]
        for j in range(j, len(genome_list)):
            genome2 = genome_list[j]
            if i == j:
                dfP.loc[genome2, genome] = 0
                dfNP.loc[genome2, genome] = 0
            else:
                # calcul of weighted distance
                if redundancy:
                    CP = np.dot(np.minimum(genome_redondant_dic[genome], genome_redondant_dic[genome2]),
                                genome_dic[genome])
                    DP = np.dot(np.maximum(genome_redondant_dic[genome], genome_redondant_dic[genome2]),
                                genome_dic[genome])
                else:
                    CP = np.dot(np.nan_to_num(genome_dic[genome] / genome_dic[genome]), genome_dic[genome2])
                    DP = np.sum(np.maximum(genome_dic[genome], genome_dic[genome2]))
                if DP > 0:
                    JP = float(CP) / float(DP)
                else:
                    JP = 0
                dfP.loc[genome2, genome] = 1 - JP
                dfP.loc[genome, genome2] = dfP[genome][genome2]

                # Calcul of non weighted distance
                if redundancy:
                    CNP = np.sum(np.minimum(genome_redondant_dic[genome], genome_redondant_dic[genome2]))
                    DNP = np.sum(np.maximum(genome_redondant_dic[genome], genome_redondant_dic[genome2]))
                else:
                    CNP = np.dot(np.nan_to_num(genome_dic[genome] / genome_dic[genome]),
                                 np.nan_to_num(genome_dic[genome2] / genome_dic[genome2]))
                    DNP = np.sum(np.maximum(np.nan_to_num(genome_dic[genome] / genome_dic[genome]),
                                            np.nan_to_num(genome_dic[genome2] / genome_dic[genome2])))
                if DNP > 0:
                    JNP = float(CNP) / float(DNP)
                else:
                    JNP = 0
                dfNP.loc[genome2, genome] = 1 - JNP
                dfNP.loc[genome, genome2] = dfNP[genome][genome2]

        i += 1
    dfP.to_csv(project_directory + "/" + "weighted_distance_matrix.csv")
    calcul_arbre_nj(project_directory + "/" + "weighted_distance_matrix.csv",
                    project_directory + "/" + "weighted_nj_tree", color_file_name)
    dfNP.to_csv(project_directory + "/" + "non_weighted_distance_matrix.csv")
    calcul_arbre_nj(project_directory + "/" + "non_weighted_distance_matrix.csv",
                    project_directory + "/" + "non_weighted_nj_tree", color_file_name)
    return dfP, dfNP


"""
Main function to load a GFA graph and compute Jaccard distances
This function also builds a PHYLIP-format matrix containing, for each genome, the presence (1) or absence (0) 
of each selected node.

Parameters
    @filename: name of the GFA file containing the pangenome
    @project_directory: will create a directory with this name to store the results
    @raxml_sample_percent : specifies the desired number of nodes in percent of total nodes
        This value is set to 0.1% by default
    @redundancy: if True, nodes present multiple times in a sample are considered; otherwise, 
        only the presence or absence of the node is used
    @strand: if True, strand information is taken into account; otherwise,
        only node presence is considered, regardless of direction
    @color_filename: this file allows coloring of the generated tree leaves
        The file must be in CSV format, separated by ",", and must contain a header row
        Headers must include at least:
            - sample: sample names that must match genome names in the pangenome
            - category: category name that will be added to leaf labels
            - color: color code associated with each sample

Returns
Distance computation methods:
    - Weighted mode: the Jaccard distance takes node size into account — larger nodes have more weight
    - Unweighted mode: standard Jaccard distance based on node presence/absence
    - Redundancy handling: if enabled, the number of occurrences of a node per genome is used.
        For example, if one genome passes through a node twice and another only once, their distance will not be zero
    - Strand handling: Jaccard distance is computed based on the direction in which nodes are traversed
        In this case, if one genome passes through s1+ and another through s1-, the distance will not be zero

"""


def analyser_pangenome(file_name, project_directory, project_name, raxml_sample_percent=0.1,
                       redundancy=True, strand=True, color_file_name=None, chromosome_file=None, masked_nodes_file_name = None, private_haplotypes_filter = 0):
    print("Launch with args : \nfile_name : " + str(file_name)
          + "\nproject_directory : " + str(project_directory)
          + "\nproject_name : " + str(project_name)
          + "\nnodes sample (%) : " + str(raxml_sample_percent)
          + "\nredundancy : " + str(redundancy)
          + "\nstrand : " + str(strand)
          + "\ncolor filename : " + str(color_file_name)
          + "\nchromosome : " + str(chromosome_file)
          + "\nmasked_nodes_file_name : " + str(masked_nodes_file_name)
          + "\nprivate_haplotypes_filter : " + str(private_haplotypes_filter)
          )
    start_time = time.time()
    rep = project_directory + "/" + project_name
    if not os.path.exists(rep):
        os.makedirs(rep)

    raxml_dir = rep + "/raxml"
    if not os.path.exists(raxml_dir):
        os.mkdir(raxml_dir)

    matrice_pav = rep + "/pav_matrix.phy"

    genome_dic, genome_redondant_dic, pav_dic, stats = load_gfa(file_name, raxml_sample_percent, strand,
                                                                           chromosome_file, masked_nodes_file_name,private_haplotypes_filter)

    export_stats(rep + "/stats.csv", stats)

    pav_to_phylip(pav_dic, matrice_pav)
    # Calcul des distances de Jaccard entre les paires
    dfP, dfNP = calculate_distance_jaccard(genome_dic, genome_redondant_dic, rep, redundancy, color_file_name)
    # conversion de la matrice de distance csv vers le format phylip
    distance_matrix_to_phylip(dfP, rep + "/weighted_distance_matrix.phy")
    distance_matrix_to_phylip(dfNP, rep + "/non_weighted_distance_matrix.phy")
    # RAxML command
    raxml_command = [
        "raxmlHPC",  # Le nom de l'exécutable RAxML
        "-s", "../pav_matrix.phy",  # Le fichier d'entrée
        "-m", "BINGAMMA",  # Le modèle d'évolution (ici binaire)
        "-p", "12345",  # Graine aléatoire pour la reproductibilité
        # '-#', '100',  # Nombre d'itérations pour bootstrap
        "-n", project_name  # Nom du fichier de sortie
    ]

    # launching RaxML command
    subprocess.run(raxml_command, cwd=raxml_dir)
    newick_file_name = raxml_dir + "/RAxML_bestTree." + project_name
    plot_newick(newick_file_name, rep + "/tree_raxml.png", color_file_name)
    print("Distances computed in " + str(time.time()-start_time) + " s")


"""
Fonction to convert PAV matrix to phylip format

Parameters
    @pav_matrix_dic : PAV matrix (dictionnary)
    @matrice_distance_phylip_filename : output filename (PHYLIP)
returns
"""


def pav_to_phylip(pav_matrix_dic, matrice_distance_phylip_filename):
    # Creating phylip file
    file = open(matrice_distance_phylip_filename, "w")

    with file:
        if len(list(pav_matrix_dic.keys())) > 0:
            entete = str(len(list(pav_matrix_dic.keys()))) + " " + str(
                len(pav_matrix_dic[list(pav_matrix_dic.keys())[0]])) + "\n"
            file.write(entete)
            for item in pav_matrix_dic:
                ligne = str(item) + " "
                for p in pav_matrix_dic[item]:
                    ligne += str(p)
                ligne += "\n"
                file.write(ligne)

    file.close()


"""
Fonction to convert distance matrix dataframe to phylip format

Parameters
    @pav_matrix_dic : PAV matrix (dictionnary)
    @phylip_distance_matrix_file_name : output filename (PHYLIP)
returns
"""


def distance_matrix_to_phylip(df, phylip_distance_matrix_file_name):
    sample_list = list(df.columns)
    file_dm = open(phylip_distance_matrix_file_name, "w")
    with file_dm:
        entete = str(len(sample_list)) + "\n"
        file_dm.write(entete)
        for echantillon in sample_list:
            line = str(echantillon) + " "
            for dist in df[echantillon]:
                line += str(dist) + " "
            line += "\n"
            file_dm.write(line)
    file_dm.close()


"""
Function to plot dendrogram

Parameters
    @newick_file: tree to be plotted in Newick format
    @output_file: name of the output file
    @color_filename: this file allows coloring of the leaves in the generated tree
        The file must be in CSV format, separated by "," and must contain a header row
        Headers must include at least:
            - sample: sample names that must match the genome names in the pangenome
            - category: category name to be added to the leaf labels
            - color: color code associated with each sample
"""


def plot_newick(newick_file, output_file, color_filename=None):
    # Load newick tree
    tree = Phylo.read(newick_file, "newick")

    label_colors = {}  # mapping label -> color

    if color_filename and os.path.exists(color_filename):
        df = pd.read_csv(color_filename)

        # Create dictionnary sample → color and category
        colors = dict(zip(df['sample'].astype(str), df['color']))
        categories = dict(zip(df['sample'].astype(str), df['category']))

        # Modify leafs names to add category
        for clade in tree.get_terminals():
            orig_name = clade.name
            if orig_name in categories:
                new_label = f"{orig_name} - {categories[orig_name]}"
                clade.name = new_label
                label_colors[new_label] = colors.get(orig_name, "black")
            else:
                label_colors[clade.name] = "black"
        for clade in tree.get_nonterminals():
            clade.name = None
    else:
        for clade in tree.get_terminals():
            label_colors[clade.name] = "black"
        for clade in tree.get_nonterminals():
            clade.name = None

    # Draw tree with labels and color
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(
        tree,
        axes=ax,
        label_colors=label_colors,
        do_show=False
    )

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


"""
Function to compute a phylogenetic tree from a Jaccard distance matrix

Parameters:
    @matrice_distance_filename: name of the file containing the distance matrix
    @file_name: name of the output file
    @color_filename: file used to colorize the tree

Returns:
    - Tree in Newick format
"""


def calcul_arbre_nj(matrice_distance_filename, file_name, color_filename=None):
    df = pd.read_csv(matrice_distance_filename, index_col=0)
    # dataframe to lower distance matrix conversion
    df_val = df.values
    triangulaire_inf = []
    for i in range(0, df_val.shape[0]):
        triangulaire_inf.append([])
        for j in range(i + 1):
            triangulaire_inf[i].append(df_val[i, j])

    dm = Phylo.TreeConstruction.DistanceMatrix(list(df.columns), triangulaire_inf)
    # computing tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    newick_tree = tree.format('newick')
    Phylo.write(tree, file_name + ".nwk", 'newick')
    plot_newick(file_name + ".nwk", file_name + ".png", color_filename)
    return newick_tree



'''
Fonction to reverse inversed haplotype (chromosomes with more reversed strands than direct strands)
The Nodes that are shared with other chromosomes (not detected as reversed) are not reversed
'''


def inverse_reverse_haplotype(file_name, output_file_name):
    file = open(file_name, "r")
    output_file = open(output_file_name, "w")
    nb_chemins = 0

    # Definition of P lines or W lines separators
    # sep[0] for P lines
    # sep[1] for W lines
    sep = ["[,;.*]", "(<|>)"]
    dic_count_direct_reverse_strand = {}
    conserved_nodes_set = set()
    temps_depart = time.time()
    with file, output_file:
        total_path = sum(1 for line in file if line.startswith(('P', 'W')))
        print("Debut du parsing, nombre de chemins : " + str(total_path))
        file.seek(0, 0)
        line = file.readline()
        print("File opening " + str(file_name))

        #  First GFA browsing to find nodes to create the data frame and to get the nodes size
        tps1 = time.time()

        with tqdm(total=total_path) as bar:
            while line:
                ligne_dec = line.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    # First GFA browsing => get haplotypes with higher reversed nodes number than direct nodes
                    # These haplotypes will be inversed after
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0

                    else:
                        ind = 6
                        walk = 1

                    chromosome, genome = get_chromosome_genome(line)

                    # the dic_count_direct_reverse_strand dictionary is used to detect haplotypes
                    # for which there are more reverse nodes than direct and which should probably be reversed.
                    if genome not in dic_count_direct_reverse_strand:
                        dic_count_direct_reverse_strand[genome] = {}

                    if chromosome not in dic_count_direct_reverse_strand[genome]:
                        dic_count_direct_reverse_strand[genome][chromosome] = {"+": 0, "-": 0}

                    if ligne_dec[0] == 'P':
                        dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count("+")
                        dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("-")
                    else:
                        dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count(">")
                        dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("<")

                line = file.readline()
        dic_to_check = {}
        # Search for inversed haplotypes
        for g in dic_count_direct_reverse_strand:
            for c in dic_count_direct_reverse_strand[g]:
                if dic_count_direct_reverse_strand[g][c]["-"] > dic_count_direct_reverse_strand[g][c]["+"]:
                    if g in dic_to_check:
                        dic_to_check[g].append(c)
                    else:
                        dic_to_check[g] = [c]
        print("List of inversed haplotypes : " + str(dic_to_check))
        # 2nd GFA browsing to get nodes and strand on inversed chromosome
        file.seek(0, 0)
        line = file.readline()
        with tqdm(total=total_path) as bar:
            while line:
                ligne_dec = line.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0

                    else:
                        ind = 6
                        walk = 1

                    chromosome, genome = get_chromosome_genome(line)

                    if genome in dic_to_check and chromosome in dic_to_check[genome]:
                        nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                        i = 0
                        while i < len(nodes_list_):
                            node = ""
                            strand = ""
                            if walk == 1 and nodes_list_[i] in ["<", ">"]:
                                node = nodes_list_[i + 1]
                                strand = nodes_list_[i]
                                i += 2
                            else:
                                if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                                    node = nodes_list_[i][:-1]
                                    strand = nodes_list_[i][-1]
                                i += 1

                            if node != "" and strand != "" and node not in conserved_nodes_set:
                                conserved_nodes_set.add(node)
                line = file.readline()

        # 3rd GFA browsing to get uniques nodes
        file.seek(0, 0)
        ligne = file.readline()
        redundancy_nodes_set = set()
        with tqdm(total=total_path) as bar:
            while ligne:
                ligne_dec = ligne.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0

                    else:
                        ind = 6
                        walk = 1

                    chromosome, genome = get_chromosome_genome(ligne)

                    if genome not in dic_to_check or (
                            genome in dic_to_check and chromosome not in dic_to_check[genome]):
                        nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                        i = 0
                        while i < len(nodes_list_):
                            node = ""
                            strand = ""
                            if walk == 1 and nodes_list_[i] in ["<", ">"]:
                                node = nodes_list_[i + 1]
                                i += 2
                            else:
                                if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                                    node = nodes_list_[i][:-1]
                                i += 1
                            # Check if node is shared with a non inversed haplotype, if it is then this node will
                            # be removed form the list of nodes to inverse (nodes shared between only inversed haplotypes are inversed)
                            if node in conserved_nodes_set:
                                redundancy_nodes_set.add(node)
                ligne = file.readline()
        print("Number of nodes to inverse in the paths : " + str(len(conserved_nodes_set)))
        conserved_nodes_set = conserved_nodes_set - redundancy_nodes_set

        print("Number of unique node to inverse in nodes and links : " + str(len(conserved_nodes_set)))

        print("Number of paths : " + str(nb_chemins))

        tps2 = time.time()
        print("Time parsing nodes : " + str(tps2 - tps1))

        # Inversion du GFA
        print("GFA construction")
        total_lines = sum(1 for line in file)
        file.seek(0, 0)
        line = file.readline()
        inversed_nodes_number = 0
        nb_noeuds_conserves = 0
        with tqdm(total=total_lines) as bar:
            while line:
                ligne_to_print = line
                ligne_dec = line.split()
                if ligne_dec[0] == 'S':
                    if ligne_dec[1] in conserved_nodes_set:
                        # Inversing node's sequence
                        ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t" + str(
                            Seq(ligne_dec[2]).reverse_complement())
                        for j in range(3, len(ligne_dec)):
                            if (j < len(ligne_dec) - 1):
                                ligne_to_print += ligne_dec[j] + "\t"
                            else:
                                ligne_to_print += ligne_dec[j]
                        ligne_to_print += "\n"
                if ligne_dec[0] == 'L':
                    if ligne_dec[1] in conserved_nodes_set or ligne_dec[3] in conserved_nodes_set:
                        ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t"
                        if ligne_dec[1] in conserved_nodes_set:
                            if ligne_dec[2] == "+":
                                ligne_to_print += "-"
                            else:
                                ligne_to_print += "+"
                        else:
                            ligne_to_print += ligne_dec[2]
                        ligne_to_print += "\t" + ligne_dec[3] + "\t"
                        if ligne_dec[3] in conserved_nodes_set:
                            if ligne_dec[4] == "+":
                                ligne_to_print += "-"
                            else:
                                ligne_to_print += "+"
                        else:
                            ligne_to_print += ligne_dec[4]

                        for j in range(5, len(ligne_dec)):
                            ligne_to_print += "\t" + ligne_dec[j]
                        ligne_to_print += "\n"

                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0


                    else:
                        ind = 6
                        walk = 1

                    chromosome, genome = get_chromosome_genome(line)

                    if (genome in dic_to_check and chromosome in dic_to_check[genome]):
                        if walk == 1:
                            ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t" \
                                             + ligne_dec[2] + "\t" + ligne_dec[3] + "\t" \
                                             + ligne_dec[4] + "\t" + ligne_dec[5] + "\t"
                        else:
                            ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t"
                        nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                        nodes_list = []
                        liste_strand = []
                        i = 0
                        while i < len(nodes_list_):
                            if walk == 1 and nodes_list_[i] in ["<", ">"]:
                                nodes_list.append(nodes_list_[i + 1])
                                liste_strand.append(nodes_list_[i])
                                i += 2
                            else:
                                if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                                    nodes_list.append(nodes_list_[i][:-1])
                                    liste_strand.append(nodes_list_[i][-1])
                                i += 1

                        for i in range(len(nodes_list) - 1, -1, -1):
                            if walk == 1:
                                if nodes_list[i] not in conserved_nodes_set:
                                    inversed_nodes_number += 1
                                    if (liste_strand[i] == "<"):
                                        ligne_to_print += ">"
                                    else:
                                        ligne_to_print += "<"
                                else:
                                    nb_noeuds_conserves += 1
                                    ligne_to_print += liste_strand[i]
                                ligne_to_print += nodes_list[i]
                            else:
                                ligne_to_print += nodes_list[i]
                                if nodes_list[i] not in conserved_nodes_set:
                                    inversed_nodes_number += 1
                                    if (liste_strand[i] == "+"):
                                        ligne_to_print += "-"
                                    else:
                                        ligne_to_print += "+"
                                else:
                                    nb_noeuds_conserves += 1
                                    ligne_to_print += liste_strand[i]
                                if i > 0:
                                    ligne_to_print += ","
                        ligne_to_print += "\n"

                output_file.write(ligne_to_print)
                line = file.readline()
            bar.update(1)
        print(
            "GFA generated, total time : " + str(time.time() - temps_depart) + " Inversed nodes number : " + str(
                inversed_nodes_number) + " Nb of conserved nodes : " + str(nb_noeuds_conserves))
        file.close()
        output_file.close()
        return conserved_nodes_set, dic_count_direct_reverse_strand


'''
Fonction to reverse a GFA (+ or > becomes - or < and link / node sequences are inversed)
'''


def reverse_gfa(file_name, output_file_name):
    file = open(file_name, "r")
    output_file = open(output_file_name, "w")

    # Definition of P lines or W lines separators
    # sep[0] for P lines
    # sep[1] for W lines
    sep = ["[,;.*]", "(<|>)"]
    inversed_nodes_number = 0
    with file, output_file:
        print("Start parsing")
        temps_depart = time.time()
        line = file.readline()
        while line:
            ligne_to_print = line
            ligne_dec = line.split()
            if ligne_dec[0] == 'S':
                # Inversing node's sequence
                ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t" + str(Seq(ligne_dec[2]).reverse_complement())
                for j in range(3, len(ligne_dec)):
                    if (j < len(ligne_dec) - 1):
                        ligne_to_print += ligne_dec[j] + "\t"
                    else:
                        ligne_to_print += ligne_dec[j]
                ligne_to_print += "\n"
            if ligne_dec[0] == 'L':
                ligne_to_print = ligne_to_print.replace("+", "/plus/").replace("-", "+").replace("/plus/", "-")

            if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                if ligne_dec[0] == 'P':
                    ind = 2
                    walk = 0
                else:
                    ind = 6
                    walk = 1

                if walk == 1:
                    ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t" \
                                     + ligne_dec[2] + "\t" + ligne_dec[3] + "\t" \
                                     + ligne_dec[4] + "\t" + ligne_dec[5] + "\t"
                else:
                    ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t"
                nodes_list_ = re.split(sep[walk], ligne_dec[ind])
                nodes_list = []
                liste_strand = []
                i = 0
                while i < len(nodes_list_):
                    if walk == 1 and nodes_list_[i] in ["<", ">"]:
                        nodes_list.append(nodes_list_[i + 1])
                        liste_strand.append(nodes_list_[i])
                        i += 2
                    else:
                        if walk == 0 and nodes_list_[i][-1] in ["+", "-"]:
                            nodes_list.append(nodes_list_[i][:-1])
                            liste_strand.append(nodes_list_[i][-1])
                        i += 1

                for i in range(len(nodes_list) - 1, -1, -1):
                    if walk == 1:
                        inversed_nodes_number += 1
                        if (liste_strand[i] == "<"):
                            ligne_to_print += ">"
                        else:
                            ligne_to_print += "<"
                        ligne_to_print += nodes_list[i]
                    else:
                        ligne_to_print += nodes_list[i]
                        inversed_nodes_number += 1
                        if (liste_strand[i] == "+"):
                            ligne_to_print += "-"
                        else:
                            ligne_to_print += "+"
                        if i > 0:
                            ligne_to_print += ","
                for j in range(ind + 1, len(ligne_dec)):
                    ligne_to_print += "\t" + ligne_dec[j]
                ligne_to_print += "\n"

            output_file.write(ligne_to_print)
            line = file.readline()
        print("GFA generated, total time : " + str(time.time() - temps_depart) + "\nNumber of nodes inversed : " + str(
            inversed_nodes_number))
        file.close()
        output_file.close()
