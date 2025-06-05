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
from ete3 import Tree, TextFace
import os
from tqdm import tqdm
import subprocess
import re
import time
import sys
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
    @nb_noeuds_arbre_objectif: specifies the desired number of nodes
        This value is set to 10,000 by default, and the actual number of returned nodes will
        be slightly different
    @type_selection="random" or other: node selection mode for the PAV matrix, either selected randomly or by size (selecting the largest nodes)
    @strand: True to use strand, False otherwise

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
def charger_fichier_gfa(file_name, nb_noeuds_arbre_objectif=10000, type_selection="random", strand=True):
    file = open(file_name, "r")
    nb_noeuds = 0
    nb_liens = 0
    nb_chemins = 0
    stats = {}
    nb_moins = 0
    nb_plus = 0
    
    #Définition des séparateurs possible selon le type de chemin (P ou W) et si on garde le strand ou non
    #sep[0] dans un chemin de type Path
    #sep[1] dans un chemin de type Walk
    sep = ["[,;.*]","(<|>)"]
    dic_count_direct_reverse_strand = {}
    dic_reversed = {}
    with file:
        ligne = file.readline()
        print("File opening " + str(file_name))
        genome_dic = {}
        genome_redondant_dic= {}
        node_lentgh_dic = {}
        node_index = {}
        pav_dic = {}
        liste_noeuds_raxml = []
        index_noeuds_raxml_dic = {}
        longueurs_noeuds = []
        #Premier parcours du GFA pour trouver les noeuds
        #Ce parcours sert à compter les noeuds pour ensuite créer les structures à la taille adaptée
        #Il sert également à stocker la taille des noeuds dans le cas où l'on séelctionne les noeuds
        #sur leur taille
        tps1 = time.time()
        while ligne :
            ligne_dec = ligne.split()
            if ligne_dec[0] == 'S' and ligne_dec[1] not in node_lentgh_dic:
                node_lentgh_dic[ligne_dec[1]] = len(ligne_dec[2])
                longueurs_noeuds.append([ligne_dec[1], len(ligne_dec[2])])
                node_index[ligne_dec[1]] = nb_noeuds
                nb_noeuds += 1

            if ligne_dec[0] == 'L':
                nb_liens += 1
            
            if ligne_dec[0]=='P' or ligne_dec[0]=='W':
                nb_chemins += 1

                #Pour le premier parcours du fichier on va préparer un dictionnaire
                #pour chaque génome et chaque chromosome on va compter le nombre de 
                #strand direct et reverse de façon à inverser si le nombre de reverse
                #est supérieur au nombre de direct => détection des chromosomes potentiellement inversés
                chromosome = ""
                if ligne_dec[0]=='P' :
                    ind = 2
                    if (len(ligne_dec[1].split("#")) > 1) :
                        chromosome = str(ligne_dec[1].split("#")[1])
                    else :
                        if (len(ligne_dec[1].split(".")) > 1) :
                            chromosome = str(ligne_dec[1].split(".")[1])
                        else : 
                            chromosome = "0"
                    
                    #Les noms des génomes sont composés de façon différentes
                    #Ils peuvent contenir un ensemble de contigs et on va devoir les regrouper
                    #Le regroupement se fait à la racine du nom (on coupe sur le séparateur # ou .)
                    if (len (ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]
                    
                else:
                    ind = 6
                    chromosome = str(ligne_dec[3])
                    genome = ligne_dec[1]+"_"+ligne_dec[2]
                
                

                
                
                #préparation de la structure dic_count_direct_reverse_strand
                #on va compter le nombre de strand direct et reverse
                #si le reverse est supérieur au direct on inversera le chromosome
                #c'est du au fait que le séquençage des chromosomes n'est pas forcément dans le bon sens
                if genome not in dic_count_direct_reverse_strand :
                    dic_count_direct_reverse_strand[genome] = {}
                
                if chromosome not in dic_count_direct_reverse_strand[genome] :
                    dic_count_direct_reverse_strand[genome][chromosome] = {"+":0, "-":0}
                
                if ligne_dec[0]=='P' :
                    dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count("+")
                    dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("-")
                else:
                    dic_count_direct_reverse_strand[genome][chromosome]["+"] += ligne_dec[ind].count(">")
                    dic_count_direct_reverse_strand[genome][chromosome]["-"] += ligne_dec[ind].count("<")
                    
            
            ligne = file.readline()
        
        
        dic_to_check = {}

        #Recherche des genomes / chromosomes inversés
        for g in dic_count_direct_reverse_strand :
            for c in dic_count_direct_reverse_strand[g]:
                if dic_count_direct_reverse_strand[g][c]["-"] > dic_count_direct_reverse_strand[g][c]["+"]:
                    if g in dic_to_check :
                        dic_to_check[g].append(c)
                    else :
                        dic_to_check[g] = [c]
        
        print("Inversed Genomes / chromosomes detected : " + str(dic_to_check))
        print("Number of paths : " + str(nb_chemins))
        #On trie la liste des longueurs de noeud
        #la liste triée va servir à récupérer la taille à partir de laquelle
        #on va sélectionner les noeuds pour obtenir un nombre proche de l'objectif du 
        #nombre de noeuds
        
        longueurs_noeuds.sort(key=lambda x: x[1])
        
        
        if nb_noeuds_arbre_objectif > nb_noeuds :
            nb_noeuds_arbre_objectif = nb_noeuds
        
        if type_selection=="random" :
            liste_noeuds_raxml = rand.sample(list(node_lentgh_dic.keys()), nb_noeuds_arbre_objectif) 
        else:
            liste_noeuds_raxml = [longueurs_noeuds[i][0] for i in range(len(longueurs_noeuds)-1, len(longueurs_noeuds)-1-nb_noeuds_arbre_objectif,-1)]
        
        nb_noeuds_arbre = len(liste_noeuds_raxml)
        i = 0
        for noeud in liste_noeuds_raxml:
            index_noeuds_raxml_dic[noeud] = i
            i += 1
         
        liste_noeuds_raxml = []
        longueurs_noeuds = []
        
        print("Number of selected nodes to compute tree : " + str(nb_noeuds_arbre))
        tps2 = time.time()
        print("Time parsing nodes : " + str(tps2 - tps1))
        #Second parcours du GFA pour lire les chemins et construire la structure renvoyée
        file.seek(0,0)
        ligne = file.readline()
        num_chemin = 0
        nb_noeuds_chemin_max = 0

            
        with tqdm(total=nb_chemins) as bar:
            while ligne :
                ligne_dec = ligne.split()
                if ligne_dec[0]=='P' or ligne_dec[0]=='W':
                    bar.update(1)
                    num_chemin += 1

                    walk = 0
                    


                    chromosome = ""
                    #Le fichier GFA a plusieurs structure pour les chemins
                    #on vérifie si les chemins sont indiqués sous la balise P ou W
                    #le chemin est ensuite défini dans une partie différente de la ligne selon le cas
                    if ligne_dec[0]=='P' :
                        ind = 2
                        walk = 0

                        #Les noms des génomes sont composés de façon différentes
                        #Ils peuvent contenir un ensemble de contigs et on va devoir les regrouper
                        #Le regroupement se fait à la racine du nom (on coupe sur le séparateur # ou .)
                        if (len (ligne_dec[1].split("#")) > 1):
                            genome = ligne_dec[1].split("#")[0]
                        else:
                            genome = ligne_dec[1].split(".")[0]

                        if (len(ligne_dec[1].split("#")) > 1) :
                            chromosome = str(ligne_dec[1].split("#")[1])
                        else :
                            if (len(ligne_dec[1].split(".")) > 1) :
                                chromosome = str(ligne_dec[1].split(".")[1])
                            else : 
                                chromosome = "0"
                    else:
                        ind = 6
                        walk = 1
                        chromosome = str(ligne_dec[3])
                        genome = ligne_dec[1]+"_"+ligne_dec[2]
                    
                    
                    if genome not in genome_dic :
                        if strand :
                            genome_dic[genome] = np.zeros(2*nb_noeuds)
                            genome_redondant_dic[genome]=np.zeros(2*nb_noeuds)
                            pav_dic[genome] = np.zeros(2*nb_noeuds_arbre, dtype=int)
                        else :
                            genome_dic[genome] = np.zeros(nb_noeuds)
                            genome_redondant_dic[genome]=np.zeros(nb_noeuds)
                            pav_dic[genome] = np.zeros(nb_noeuds_arbre, dtype=int)
                            
                    if genome not in stats:
                        stats[genome] = {"total+" : 0, "total-" : 0}
                    if chromosome not in stats[genome]:
                        stats[genome][chromosome] =  {"+" : 0, "-" : 0}


                    liste_noeuds_ = re.split(sep[walk],ligne_dec[ind])
                    liste_strand = []
                    liste_noeuds = []
                    i = 0
                    while i < len(liste_noeuds_) :
                        if walk == 1 and liste_noeuds_[i] in ["<", ">"]:
                            liste_strand.append(liste_noeuds_[i])
                            liste_noeuds.append(liste_noeuds_[i+1])
                            i += 2
                        else :
                            if walk == 0 and liste_noeuds_[i][-1] in ["+", "-"]:
                                liste_strand.append(liste_noeuds_[i][-1])
                                liste_noeuds.append(liste_noeuds_[i][:-1])
                                i += 1  
                            else :
                                i += 1
  

                    for n in range(0,len(liste_noeuds)) :
                        noeud = liste_noeuds[n]
                        s = liste_strand[n]
                        
                        if noeud != "":
                            if strand :
                                if s in ["<", "-"]:
                                    genome_redondant_dic[genome][node_index[noeud]] += 1
                                    genome_dic[genome][node_index[noeud]] = node_lentgh_dic[noeud]
                                    stats[genome][chromosome]["-"] += 1
                                    stats[genome]["total-"] += 1
                                    nb_moins += 1
                                else:
                                    genome_redondant_dic[genome][node_index[noeud]+nb_noeuds] += 1
                                    genome_dic[genome][node_index[noeud]+nb_noeuds] = node_lentgh_dic[noeud]
                                    stats[genome][chromosome]["+"] += 1
                                    stats[genome]["total+"] += 1
                                    nb_plus += 1
                            else :
                                genome_redondant_dic[genome][node_index[noeud]] += 1
                                genome_dic[genome][node_index[noeud]] = node_lentgh_dic[noeud]
                                if s in ["<", "-"]:
                                    stats[genome][chromosome]["-"] += 1
                                    nb_moins += 1
                                else:
                                    stats[genome][chromosome]["+"] += 1
                                    nb_plus += 1
                            if noeud in index_noeuds_raxml_dic :
                                if strand :
                                    if s in ["<", "-"]:
                                        pav_dic[genome][index_noeuds_raxml_dic[noeud]] = int(1)
                                    else:
                                        pav_dic[genome][index_noeuds_raxml_dic[noeud]+nb_noeuds_arbre] = int(1)
                                else:
                                    pav_dic[genome][index_noeuds_raxml_dic[noeud]] = int(1)
                    
                    if nb_noeuds_chemin_max < len(liste_noeuds) :
                         nb_noeuds_chemin_max = len(liste_noeuds)
    
                ligne = file.readline()
                
    file.close() 
    print("Number of nodes : " + str(nb_noeuds) + "\nLinks : " + str(nb_liens) + "\nPaths : " + str(nb_chemins) + "\nMaximum number of nodes in a path : " + str(nb_noeuds_chemin_max))
    tps3 = time.time()
    print("Time parsing paths : " + str(tps3 - tps2) + "Total time : " + str(tps3 - tps1))
    print("Strand's statistics : reverse : " + str(nb_moins) + " direct : " + str(nb_plus))
    return genome_dic, genome_redondant_dic, pav_dic, stats

"""
Function to save stats in csv file
Parameters
    @file_name : name of output file
    @stats : dictionnary from fonction charger_fichier_gfa
Returns 
"""
def export_stats(file_name, stats):
    file = open(file_name, "w")
    with file:
        file.write("genome,chromosome,direct strand number, reverse strandnumber\n")
        for g in stats :
            for chr in stats[g]:
                if chr != "total+" and chr != "total-":
                    file.write(g+","+str(chr)+","+str(stats[g][chr]["+"])+","+str(stats[g][chr]["-"])+"\n")
    file.close()



"""
This function computes the Jaccard distance from the structure returned by the GFA

Parameters
    @genome_dic: the structure returned by the GFA loading function 
        (dictionary containing genomes, the nodes present in each genome, and their size)
    @genome_redondant_dic: the structure returned by the GFA loading function 
        (dictionary containing genomes, the nodes present in each genome, and their number of occurrences)
    @project_directory: directory where output files will be stored
    @redondance: if True, Jaccard distance is calculated on all nodes including redundant ones; 
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
def calculate_distance_jaccard(genome_dic, genome_redondant_dic, project_directory, redondance=True, color_file_name=None):
    print("Computing the distances")
    genome_list = list(genome_dic.keys())
    dfP = pd.DataFrame(np.ones(shape=(len(genome_list),len(genome_list))), columns=genome_list, index=genome_list)
    dfNP = pd.DataFrame(np.ones(shape=(len(genome_list),len(genome_list))), columns=genome_list, index=genome_list)
    print("Genomes number : " + str(len(genome_list)))
    for i in range (0,len(genome_list)):
        print("Genome : " + str(i) + " name : " + genome_list[i])
        j = i
        genome = genome_list[i]
        for j in range (j, len(genome_list)):
            genome2 = genome_list[j]
            if i == j :
                dfP.loc[genome2,genome] = 0
                dfNP.loc[genome2,genome] = 0
            else :
                #calcul de la distance pondérée
                if redondance :
                    CP = np.dot(np.minimum(genome_redondant_dic[genome],genome_redondant_dic[genome2]),genome_dic[genome])
                    DP = np.dot(np.maximum(genome_redondant_dic[genome],genome_redondant_dic[genome2]),genome_dic[genome])
                else:
                    CP = np.dot(np.nan_to_num(genome_dic[genome]/genome_dic[genome]),genome_dic[genome2])
                    DP = np.sum(np.maximum(genome_dic[genome],genome_dic[genome2]))
                if DP > 0:    
                    JP = float(CP)/float(DP)
                else:
                    JP = 0
                dfP.loc[genome2,genome] = 1 - JP
                dfP.loc[genome,genome2] = dfP[genome][genome2]

                #Calcul de la distance non pondérée
                if redondance:
                    CNP = np.sum(np.minimum(genome_redondant_dic[genome],genome_redondant_dic[genome2]))
                    DNP = np.sum(np.maximum(genome_redondant_dic[genome],genome_redondant_dic[genome2]))
                else :
                    CNP = np.dot(np.nan_to_num(genome_dic[genome]/genome_dic[genome]),np.nan_to_num(genome_dic[genome2]/genome_dic[genome2]))
                    DNP = np.sum(np.maximum(np.nan_to_num(genome_dic[genome]/genome_dic[genome]),np.nan_to_num(genome_dic[genome2]/genome_dic[genome2])))
                if DNP > 0 :
                    JNP = float(CNP)/float(DNP)
                else :
                    JNP = 0
                dfNP.loc[genome2,genome] = 1 - JNP
                dfNP.loc[genome,genome2] = dfNP[genome][genome2]

        i += 1
    dfP.to_csv(project_directory+"/"+"weighted_distance_matrix.csv")
    calcul_arbre_nj(project_directory+"/"+"weighted_distance_matrix.csv", project_directory+"/"+"weighted_nj_tree", color_file_name)
    dfNP.to_csv(project_directory+"/"+"non_weighted_distance_matrix.csv")
    calcul_arbre_nj(project_directory+"/"+"non_weighted_distance_matrix.csv", project_directory+"/"+"non_weighted_nj_tree", color_file_name)
    print("Distances computed")
    return dfP, dfNP

"""
Main function to load a GFA graph and compute Jaccard distances
This function also builds a PHYLIP-format matrix containing, for each genome, the presence (1) or absence (0) 
of each selected node.

Parameters
    @filename: name of the GFA file containing the pangenome
    @project_directory: will create a directory with this name to store the results
    @nb_noeuds_cible: approximate number of nodes desired in the presence/absence matrix
    @methode: "random" or "size", selects the nodes for analysis (either randomly or by choosing the largest nodes)
    @redondance: if True, nodes present multiple times in a sample are considered; otherwise, 
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
def analyser_pangenome(file_name, project_directory, project_name, nb_noeuds_cible = 10000, methode="random", redondance = True, strand = True, color_file_name=None):

    print("Launch with args : \nfile_name : " + str(file_name)
          + "\nproject_directory : " + str(project_directory)
          + "\nproject_name : " + str(project_name)
          + "\nnodes number : " + str(nb_noeuds_cible)
          + "\nmethod : " + str(methode)
          +"\nredundancy : " + str(redondance)
          +"\nstrand : " + str(strand)
          +"\ncolor filename : " + str(color_file_name))
    
    rep = project_directory+"/"+project_name
    if not os.path.exists(rep):
        os.makedirs(rep)
        
    raxml_dir = rep+"/raxml" 
    if not os.path.exists(raxml_dir):
        os.mkdir(raxml_dir)
    
    matrice_pav = rep+"/pav_matrix.phy"
    
    genome_dic, genome_redondant_dic, pav_dic, stats = charger_fichier_gfa(file_name, nb_noeuds_cible, methode, strand)

    export_stats(rep+"/stats.csv",stats)

    pav_to_phylip(pav_dic, matrice_pav)
    #Calcul des distances de Jaccard entre les paires
    dfP, dfNP = calculate_distance_jaccard(genome_dic, genome_redondant_dic, rep, redondance, color_file_name)
    #conversion de la matrice de distance csv vers le format phylip
    distance_matrix_to_phylip(dfP,rep+"/weighted_distance_matrix.phy")
    distance_matrix_to_phylip(dfNP,rep+"/non_weighted_distance_matrix.phy")
    # Commande RAxML
    raxml_command = [
        "raxmlHPC",  # Le nom de l'exécutable RAxML
        "-s", "../pav_matrix.phy",  # Le fichier d'entrée
        "-m", "BINGAMMA",  # Le modèle d'évolution (ici binaire)
        "-p", "12345",  # Graine aléatoire pour la reproductibilité
        #'-#', '100',  # Nombre d'itérations pour bootstrap
        "-n", project_name  # Nom du fichier de sortie
    ]

    # Lancer la commande
    subprocess.run(raxml_command, cwd=raxml_dir)
    newick_file_name = raxml_dir + "/RAxML_bestTree."+project_name
    plot_newick(newick_file_name, rep+"/tree_raxml.png", color_file_name)




"""
Fonction to convert PAV matrix to phylip format

Parameters
    @pav_matrix_dic : PAV matrix (dictionnary)
    @matrice_distance_phylip_filename : output filename (PHYLIP)
returns
"""

def pav_to_phylip(pav_matrix_dic, matrice_distance_phylip_filename):
    #Création du fichier phylip
    file = open(matrice_distance_phylip_filename, "w")

    with file : 
        if len(list(pav_matrix_dic.keys())) > 0:
            entete = str(len(list(pav_matrix_dic.keys()))) + " " + str(len(pav_matrix_dic[list(pav_matrix_dic.keys())[0]])) + "\n"
            file.write(entete)
            for item in pav_matrix_dic :
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
    liste_echantillons = list(df.columns)
    file_dm = open(phylip_distance_matrix_file_name,"w")
    with file_dm : 
        entete = str(len(liste_echantillons)) + "\n"
        file_dm.write(entete)
        for echantillon in liste_echantillons :
            ligne = str(echantillon) + " "
            for dist in df[echantillon]:
                ligne += str(dist) + " "
            ligne += "\n"
            file_dm.write(ligne)
    file_dm.close()    


"""
Function to plot a dendrogram

Parameters
    @newick_tree: tree to be plotted in Newick format
    @file_name: name of the output file
    @color_filename: this file allows coloring of the leaves in the generated tree
        The file must be in CSV format, separated by "," and must contain a header row
        Headers must include at least:
            - sample: sample names that must match the genome names in the pangenome
            - category: category name to be added to the leaf labels
            - color: color code associated with each sample
""" 
def plot_newick(newick_tree,file_name,color_filename=None):
    ete3_tree = Tree(newick_tree, format=1)
    if color_filename != None and os.path.exists(color_filename):
        cat_df = pd.read_csv(color_filename, sep=',')
        colors = {}
        category = {}
        for i in cat_df.index:
            colors[cat_df['sample'][i]] = cat_df['color'][i]
            category[cat_df['sample'][i]] = cat_df['category'][i]
        for leaf in ete3_tree:
            if leaf.is_leaf():
                name_face = TextFace(leaf.name + " - " + str(category.get(leaf.name, "none")), fgcolor=colors.get(leaf.name, "none"))
                leaf.name= None
                leaf.add_face(name_face, column=0)
                
    ete3_tree.render(file_name)

"""
Function to compute a phylogenetic tree from a Jaccard distance matrix

Parameters:
    @matrice_distance_filename: name of the file containing the distance matrix
    @file_name: name of the output file
    @color_filename: file used to colorize the tree

Returns:
    - Tree in Newick format
"""   

def calcul_arbre_nj(matrice_distance_filename, file_name,color_filename=None):
    df = pd.read_csv(matrice_distance_filename, index_col=0)
    #conversion dataframe en matrice de distance inférieure
    df_val = df.values
    triangulaire_inf = []
    for i in range(0, df_val.shape[0]):
        triangulaire_inf.append([])
        for j in range(i+1):
            triangulaire_inf[i].append(df_val[i, j])
    
    dm = Phylo.TreeConstruction.DistanceMatrix(list(df.columns), triangulaire_inf)
    #calcul de l'arbre
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    newick_tree = tree.format('newick')
    Phylo.write(tree, file_name+".nwk",'newick')
    plot_newick(newick_tree, file_name+".png", color_filename)
    return newick_tree



#Old fonction to plot dendogram depuis une matrice de distance
def plot_dendogramme(csv_df_filename, color_filename=None):
    df = pd.read_csv(csv_df_filename, index_col=0)
    #c_dist = pdist(df) # computing the distance
    #methodes possibles pour linkage : single, complete, average, weighted, centroid, median, ward
    c_link = linkage(df,  method='ward')# computing the linkage
    dendrogram(c_link,leaf_rotation=90, labels=df.columns, truncate_mode="level")
    cat_df = None
    if color_filename != None :
        cat_df = pd.read_csv(color_filename, sep=',')
        #cat = list(cat_df.iloc[:][1].unique)
        ax = plt.gca()
        xlbls = ax.get_xmajorticklabels()
        newlbls = []
        for lbl in xlbls:
            #lbl.set_color('#8dd3c7')
            #if (len(cat_df[cat_df['genome']==lbl.get_text()]['color'].values) > 0):
            if (len(cat_df[cat_df['genome']==lbl.get_text()]['pays_court'].values) > 0):
                pays = str(cat_df[cat_df['genome']==lbl.get_text()]['pays_court'].values[0])
            else :
                pays = ""
            if (len(cat_df[cat_df['genome']==lbl.get_text()]['color'].values) > 0):
                color = cat_df[cat_df['genome']==lbl.get_text()]['color'].values[0]
            else:
                color='#000000'
            lbl.set_text(lbl.get_text() + " - " + pays)
            lbl.set_color(color)
            newlbls.append(lbl)
        xlbls = ax.set_xticklabels(newlbls)
    plt.show()
    return df, cat_df



'''
Fonction to reverse inversed chromosomes (chromosomes with more reversed strands than direct strands)
The Nodes that are shared with other chromosomes (not detected as reversed) are not reversed
'''

def inverser_reverse_chromosome(file_name, output_file_name):
    file = open(file_name, "r")
    output_file = open(output_file_name, "w")
    nb_chemins = 0
    stats = {}


    # Définition des séparateurs possible selon le type de chemin (P ou W) et si on garde le strand ou non
    # sep[0] dans un chemin de type Path
    # sep[1] dans un chemin de type Walk
    sep = ["[,;.*]", "(<|>)"]
    dic_count_direct_reverse_strand = {}
    compteur_noeuds_chromosome = {}
    noeuds_a_conserver = set()
    temps_depart = time.time()
    with file, output_file:
        total_path = sum(1 for line in file if line.startswith(('P', 'W')))
        print("Debut du parsing, nombre de chemins : " + str(total_path))
        file.seek(0, 0)
        ligne = file.readline()
        print("File opening " + str(file_name))

        # Premier parcours du GFA pour trouver les noeuds
        # Ce parcours sert à compter les noeuds pour ensuite créer les structures à la taille adaptée
        # Il sert également à stocker la taille des noeuds dans le cas où l'on séelctionne les noeuds
        # sur leur taille
        tps1 = time.time()


        with tqdm(total=total_path) as bar:
            while ligne:
                ligne_dec = ligne.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    # Pour le premier parcours du fichier on va préparer un dictionnaire
                    # pour chaque génome et chaque chromosome on va compte le nombre de
                    # strand direct et reverse de façon à inverser si le nombre de reverse
                    # est supérieur au nombre de direct
                    chromosome = ""
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        if (len(ligne_dec[1].split("#")) > 1):
                            chromosome = str(ligne_dec[1].split("#")[1])
                        else:
                            chromosome = str(ligne_dec[1].split(".")[1])

                    else:
                        ind = 6
                        walk = 1
                        chromosome = str(ligne_dec[3])

                    # Les noms des génomes sont composés de façon différentes
                    # Ils peuvent contenir un ensemble de contigs et on va devoir les regrouper
                    # Le regroupement se fait à la racine du nom (on coupe sur le séparateur # ou .)
                    if (len(ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]



                    # préparation de la structure dic_count_direct_reverse_strand
                    # on va compter le nombre de strand direct et reverse
                    # si le reverse est supérieur au direct on inversera le chromosome
                    # c'est du au fait que le séquençage des chromosomes n'est pas forcément dans le bon sens
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




                ligne = file.readline()
        dic_to_check = {}
        #Recherche des genomes / chromosomes inversés
        for g in dic_count_direct_reverse_strand :
            for c in dic_count_direct_reverse_strand[g]:
                if dic_count_direct_reverse_strand[g][c]["-"] > dic_count_direct_reverse_strand[g][c]["+"]:
                    if g in dic_to_check :
                        dic_to_check[g].append(c)
                    else :
                        dic_to_check[g] = [c]
        print("Liste des chromosomes à inverser : " + str(dic_to_check))
        #Second parcours pour des questions de mémoires => on va récupérer les noeuds et le strand 
        #sur les chromosomes inversés
        file.seek(0, 0)
        ligne = file.readline()
        print("Debut du 2eme parcours")
        with tqdm(total=total_path) as bar:
            while ligne:
                ligne_dec = ligne.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    # Pour le premier parcours du fichier on va préparer un dictionnaire
                    # pour chaque génome et chaque chromosome on va compte le nombre de
                    # strand direct et reverse de façon à inverser si le nombre de reverse
                    # est supérieur au nombre de direct
                    chromosome = ""
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        if (len(ligne_dec[1].split("#")) > 1):
                            chromosome = str(ligne_dec[1].split("#")[1])
                        else:
                            chromosome = str(ligne_dec[1].split(".")[1])

                    else:
                        ind = 6
                        walk = 1
                        chromosome = str(ligne_dec[3])

                    # Les noms des génomes sont composés de façon différentes
                    # Ils peuvent contenir un ensemble de contigs et on va devoir les regrouper
                    # Le regroupement se fait à la racine du nom (on coupe sur le séparateur # ou .)
                    if (len(ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]
                    
                    set_noeuds = set()  
                    
                    if genome in dic_to_check and chromosome in dic_to_check[genome]:
                        liste_noeuds_ = re.split(sep[walk], ligne_dec[ind])
                        i = 0
                        while i < len(liste_noeuds_):
                            noeud = ""
                            strand = ""
                            if walk == 1 and liste_noeuds_[i] in ["<", ">"]:
                                noeud = liste_noeuds_[i + 1]
                                strand = liste_noeuds_[i]
                                i += 2
                            else:
                                if walk == 0 and liste_noeuds_[i][-1] in ["+", "-"]:
                                    noeud = liste_noeuds_[i][:-1]
                                    strand = liste_noeuds_[i][-1]
                                i += 1
                                
                            if noeud != "" and strand != "" and noeud not in noeuds_a_conserver :
                                noeuds_a_conserver.add(noeud)
                ligne = file.readline()

            
        #3ème parcours pour des questions de mémoires => on va vérifier les noeuds uniques
        file.seek(0, 0)
        ligne = file.readline()
        print("Debut du 3eme parcours")
        set_noeuds_redondants = set()
        with tqdm(total=total_path) as bar:
            while ligne:
                ligne_dec = ligne.split()
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    nb_chemins += 1
                    bar.update(1)
                    # Pour le premier parcours du fichier on va préparer un dictionnaire
                    # pour chaque génome et chaque chromosome on va compte le nombre de
                    # strand direct et reverse de façon à inverser si le nombre de reverse
                    # est supérieur au nombre de direct
                    chromosome = ""
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        if (len(ligne_dec[1].split("#")) > 1):
                            chromosome = str(ligne_dec[1].split("#")[1])
                        else:
                            chromosome = str(ligne_dec[1].split(".")[1])

                    else:
                        ind = 6
                        walk = 1
                        chromosome = str(ligne_dec[3])
        
                    if (len(ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]
                    if genome not in dic_to_check or (genome in dic_to_check and chromosome not in dic_to_check[genome]) :    
                        liste_noeuds_ = re.split(sep[walk], ligne_dec[ind])
                        i = 0
                        while i < len(liste_noeuds_):
                            noeud = ""
                            strand = ""
                            if walk == 1 and liste_noeuds_[i] in ["<", ">"]:
                                noeud = liste_noeuds_[i + 1]
                                strand = liste_noeuds_[i]
                                i += 2
                            else:
                                if walk == 0 and liste_noeuds_[i][-1] in ["+", "-"]:
                                    noeud = liste_noeuds_[i][:-1]
                                    strand = liste_noeuds_[i][-1]
                                i += 1    
                            #if noeud in noeuds_a_conserver and genome not in noeuds_a_conserver[noeud]["genomes"]:
                            #On regarde si le noeud est partagé avec un génomé non inversé, si oui on supprimera ce noeud
                            #de la liste des noeuds à inverser (les noeuds partagés par des chromosomes inversés sont à inverser)
                            if noeud in noeuds_a_conserver:
                                set_noeuds_redondants.add(noeud)
                ligne = file.readline()
        print("Nombre de noeuds à inverser dans le path : " + str(len(noeuds_a_conserver)))
        noeuds_a_conserver = noeuds_a_conserver - set_noeuds_redondants
        
        print("Nombre de noeuds uniques à inverser dans les noeuds et liens : " + str(len(noeuds_a_conserver)))
        
        
        
        print("Number of paths : " + str(nb_chemins))

        tps2 = time.time()
        print("Time parsing nodes : " + str(tps2 - tps1))
        
        #Inversion du GFA
        print("Génération du GFA")
        total_lignes = sum(1 for line in file )
        file.seek(0, 0)
        ligne = file.readline()
        nb_noeuds_inverses = 0
        nb_noeuds_conserves = 0
        with tqdm(total=total_lignes) as bar:
            while ligne:
                ligne_to_print = ligne
                ligne_dec = ligne.split()
                if ligne_dec[0] == 'S' :
                    if ligne_dec[1] in noeuds_a_conserver:
                        #inversion de la séquence du noeud
                        ligne_to_print = ligne_dec[0] +"\t" + ligne_dec[1] + "\t" + str(Seq(ligne_dec[2]).reverse_complement())
                        for j in range(3, len(ligne_dec)):
                            if (j < len(ligne_dec) -1):
                                ligne_to_print += ligne_dec[j] + "\t"
                            else :
                                ligne_to_print += ligne_dec[j]
                        ligne_to_print += "\n"
                if ligne_dec[0] == 'L' :
                    if ligne_dec[1] in noeuds_a_conserver  or ligne_dec[3]  in noeuds_a_conserver :
                        ligne_to_print = ligne_dec[0] + "\t" + ligne_dec[1] + "\t"
                        if ligne_dec[1] in noeuds_a_conserver :
                            if ligne_dec[2] == "+" :
                                ligne_to_print += "-"
                            else :
                                ligne_to_print += "+"
                        else :
                            ligne_to_print += ligne_dec[2]
                        ligne_to_print += "\t" + ligne_dec[3] + "\t"
                        if ligne_dec[3]  in noeuds_a_conserver :
                            if ligne_dec[4] == "+" :
                                ligne_to_print += "-"
                            else :
                                ligne_to_print += "+"
                        else :
                            ligne_to_print += ligne_dec[4]
                            
                        for j in range(5, len(ligne_dec)):
                            ligne_to_print += "\t" + ligne_dec[j]
                        ligne_to_print += "\n"
                    
                if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
                    chromosome = ""
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        if (len(ligne_dec[1].split("#")) > 1):
                            chromosome = str(ligne_dec[1].split("#")[1])
                        else:
                            chromosome = str(ligne_dec[1].split(".")[1])

                    else:
                        ind = 6
                        walk = 1
                        chromosome = str(ligne_dec[3])
        
                    if (len(ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]
                    if (genome in dic_to_check and chromosome in dic_to_check[genome]):  
                        if walk == 1 :
                            ligne_to_print = ligne_dec[0]+"\t"+ligne_dec[1] + "\t" \
                                   + ligne_dec[2] + "\t" +ligne_dec[3] + "\t" \
                                   + ligne_dec[4] + "\t" +ligne_dec[5] + "\t"
                        else :
                            ligne_to_print = ligne_dec[0]+"\t"+ligne_dec[1] + "\t"
                        liste_noeuds_ = re.split(sep[walk], ligne_dec[ind])
                        liste_noeuds = []
                        liste_strand = []
                        i = 0
                        while i < len(liste_noeuds_):
                            noeud = ""
                            strand = ""
                            if walk == 1 and liste_noeuds_[i] in ["<", ">"]:
                                liste_noeuds.append(liste_noeuds_[i + 1])
                                liste_strand.append(liste_noeuds_[i])
                                i += 2
                            else:
                                if walk == 0 and liste_noeuds_[i][-1] in ["+", "-"]:
                                    liste_noeuds.append(liste_noeuds_[i][:-1])
                                    liste_strand.append(liste_noeuds_[i][-1])
                                i += 1 
                        
                        for i in range(len(liste_noeuds)-1,-1,-1):
                            if walk == 1 :
                                if liste_noeuds[i] not in noeuds_a_conserver :
                                    nb_noeuds_inverses += 1
                                    if(liste_strand[i] ==  "<"):
                                        ligne_to_print += ">"
                                    else :
                                        ligne_to_print += "<"
                                else :
                                    nb_noeuds_conserves += 1
                                    ligne_to_print += liste_strand[i]
                                ligne_to_print += liste_noeuds[i]
                            else: 
                                ligne_to_print += liste_noeuds[i]
                                if liste_noeuds[i] not in noeuds_a_conserver :
                                    nb_noeuds_inverses += 1
                                    if(liste_strand[i] ==  "+"):
                                        ligne_to_print += "-"
                                    else :
                                        ligne_to_print += "+"
                                else :
                                    nb_noeuds_conserves += 1
                                    ligne_to_print += liste_strand[i]
                        ligne_to_print += "\n"
                                
                output_file.write(ligne_to_print)
                ligne = file.readline()
            bar.update(1)
        print("GFA généré, durée du traitement : " + str(time.time()-temps_depart) + " Nb de noeuds inversés : " + str(nb_noeuds_inverses) + " Nb de noeuds conservé : " + str(nb_noeuds_conserves))
        file.close()    
        output_file.close()
        return noeuds_a_conserver, dic_count_direct_reverse_strand



