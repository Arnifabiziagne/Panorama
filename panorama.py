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
from ete3 import Tree, TreeStyle,NodeStyle, faces, TextFace,AttrFace
import os
from tqdm import tqdm
import subprocess
import re


"""
Fonction permettant de charger un fichier GFA
La fontion va extraire les noeuds et les chemins puis retourner pour chaque genome 
La liste des noeuds traversés et leur taille
La fonction permet également de sélectiionner des noeuds soit de façon aléatoire
soit par rapport à un critère de taille et indique pour chaque génome s'il contient ou non
ces noeuds
Parameters
    @file_name : nom du fichier gfa à charger
    @nb_noeuds_arbre_objectif : permet d'indiquer le nombre de noeuds souhaité 
       cette donnée est fixée à 10000 par défaut et le nomnbre réel de noeuds renvoyé sera 
       un peu différent
    @type_selection="aleatoire" ou "taille" : mode de sélection des noeuds pour le calcul basé sur la morphologie (présence ou absence des noeuds)
    @strand : True
Returns
   genome_dic : dictionnary
       Dictionnaire avec pour clés la liste des génomes et pour valeur un tableau de la taille du nombre de noeud et pour chaque noeud la taille du noeud
   genome_redondant_dic : dictionnary
       Dictionnaire avec pour clés la liste des génomes et pour valeur un dictionnaire dont les clés = noms des noeuds et valeur = tableau de taille = nb de noeud et pour chaque noeud le nombre d'occurence du noeud
       Si le strand est utilisé alors le tableau a une taile 2 x nb_noeuds (un pour les noeuds - et un pour les noeuds +)
    pav_dic : dictionnary 
       Dictionnaire avec pour clés la liste des génomes et pour chaque génome un tableau contenant la présence (1) ou absence (0) des noeuds retenus
"""
def charger_fichier_gfa(file_name, nb_noeuds_arbre_objectif=10000, type_selection="random", strand=True):
    file = open(file_name, "r")
    nb_noeuds = 0
    nb_liens = 0
    nb_chemins = 0
    with file:
        ligne = file.readline()
        print("File opening " + str(file_name))
        genome_dic = {}
        genome_redondant_dic= {}
        node_lentgh_dic = {}
        node_index = {}
        pav_dic = {}
        arbre_noeuds_dic = {}
        longueurs_noeuds = []
        #Premier parcours du GFA pour trouver les noeuds
        #Ce parcours sert à compter les noeuds pour ensuite créer les structures à la taille adaptée
        #Il sert également à stocker la taille des noeuds dans le cas où l'on séelctionne les noeuds
        #sur leur taille
        while ligne :
            ligne_dec = ligne.split()
            if ligne_dec[0] == 'S' and ligne_dec[1] not in node_lentgh_dic:
                node_lentgh_dic[ligne_dec[1]] = len(ligne_dec[2])
                longueurs_noeuds.append(len(ligne_dec[2]))
                node_index[ligne_dec[1]] = nb_noeuds
                nb_noeuds += 1

            if ligne_dec[0] == 'L':
                nb_liens += 1
            
            if ligne_dec[0]=='P' or ligne_dec[0]=='W':
                nb_chemins += 1
            
            ligne = file.readline()
        
        print("Number of paths : " + str(nb_chemins))
        #On trie la liste des longueurs de noeud
        #la liste triée va servir à récupérer la taille à partir de laquelle
        #on va sélectionner les noeuds pour obtenir un nombre proche de l'objectif du 
        #nombre de noeuds
        longueurs_noeuds.sort()
        
        
        if nb_noeuds_arbre_objectif > nb_noeuds :
            nb_noeuds_arbre_objectif = nb_noeuds
        
        if (nb_noeuds_arbre_objectif < len(longueurs_noeuds)):
            taille_selection_noeud = longueurs_noeuds[-nb_noeuds_arbre_objectif]
        else :
            taille_selection_noeud = 0
        
        #taille_selection_noeud contient la taille minimale que le noeud doit avoir pour être retenu
        
        #second parcours du fichier pour sélectionner les noeuds pour construction d'un arbre de coalescence
        #les noeuds sont sélectionnés de façon aléatoire et le nombre attendu est en paramètre
        file.seek(0,0)
        ligne = file.readline()
        nb_noeuds_arbre = 0

        while ligne :
            ligne_dec = ligne.split()
            if ligne_dec[0] == 'S' and ligne_dec[1] not in arbre_noeuds_dic:
                #Choix des noeuds retenus pour le calcul de distance
                #si aléatoire => on tire au hasard pour savoir si le noeud est retenu ou non
                #le taux est fixé par rapport au nombre de noeuds obejctif
                if (type_selection == 'random') :
                    if rand.random() <= nb_noeuds_arbre_objectif / nb_noeuds :
                        arbre_noeuds_dic[ligne_dec[1]] = nb_noeuds_arbre
                        nb_noeuds_arbre += 1
                    
                else :
                    #Le mode n'est pas aléatoire, on va sélectionner les plus gros noeuds
                    #avec l'idée d'en retenir approximativement le nombre de noeuds objectif
                    if len(ligne_dec[2]) >= taille_selection_noeud :
                        arbre_noeuds_dic[ligne_dec[1]] = nb_noeuds_arbre
                        nb_noeuds_arbre += 1
            ligne = file.readline()
        
        print("Number of selected nodes to compute tree : " + str(nb_noeuds_arbre))
        
        #3ème parcours du GFA pour lire les chemins et construire la structure renvoyée
        file.seek(0,0)
        ligne = file.readline()
        num_chemin = 0
        nb_noeuds_chemin = 0
        nb_noeuds_chemin_max = 0
        #Définition des séparateurs possible selon le type de chemin (P ou W) et si on garde le strand ou non
        #sep[0] dans un chemin de type Path
        #sep[1] dans un chemin de type Walk

        sep = ["[,;.*]","(<|>)"]
        with tqdm(total=nb_chemins) as bar:
            while ligne :
                ligne_dec = ligne.split()
                if ligne_dec[0]=='P' or ligne_dec[0]=='W':
                    bar.update(1)
                    num_chemin += 1
                    walk = 0
                    #Le fichier GFA a plusieurs structure pour les chemins
                    #on vérifie si les chemins sont indiqués sous la balise P ou W
                    #le chemin est ensuite défini dans une partie différente de la ligne selon le cas
                    if ligne_dec[0]=='P' :
                        ind = 2
                        walk = 0
                    else:
                        ind = 6
                        walk = 1

                    #Les noms des génomes sont composés de façon différentes
                    #Ils peuvent contenir un ensemble de contigs et on va devoir les regrouper
                    #Le regroupement se fait à la racine du nom (on coupe sur le séparateur # ou .)
                    if (len (ligne_dec[1].split("#")) > 1):
                        genome = ligne_dec[1].split("#")[0]
                    else:
                        genome = ligne_dec[1].split(".")[0]
                    if genome not in genome_dic :
                        if strand :
                            genome_dic[genome] = np.zeros(2*nb_noeuds)
                            genome_redondant_dic[genome]=np.zeros(2*nb_noeuds)
                        else :
                            genome_dic[genome] = np.zeros(nb_noeuds)
                            genome_redondant_dic[genome]=np.zeros(nb_noeuds)
                        pav_dic[genome] = np.zeros(nb_noeuds_arbre, dtype=int)
                    

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
                                else:
                                    genome_redondant_dic[genome][node_index[noeud]+nb_noeuds] += 1
                                    genome_dic[genome][node_index[noeud]+nb_noeuds] = node_lentgh_dic[noeud]
                            else :
                                genome_redondant_dic[genome][node_index[noeud]] += 1
                                genome_dic[genome][node_index[noeud]] = node_lentgh_dic[noeud]
                            if noeud in arbre_noeuds_dic :
                                pav_dic[genome][arbre_noeuds_dic[noeud]] = int(1)
                    
                    if nb_noeuds_chemin_max < len(liste_noeuds) :
                         nb_noeuds_chemin_max = len(liste_noeuds)
    
                ligne = file.readline()
    file.close() 
    print("Number of nodes : " + str(nb_noeuds) + "\nLinks : " + str(nb_liens) + "\nPaths : " + str(nb_chemins) + "\nMaximum number of nodes in a path : " + str(nb_noeuds_chemin_max))

    return genome_dic, genome_redondant_dic, pav_dic

"""
Cette fonction calcule la distance de jaccard à partir de la structure renvoyée par le GFA
Parameters
    @genome_dic : la structure renvoyée par le chargement du GFA (dictionnaire contenanty les génomes et les noeuds présents dans le génome et leur taille)
    @genome_redondant_dic : la structure renvoyée par le chargement du GFA (dictionnaire contenanty les génomes et les noeuds présents dans le génome et leur nombre d'apparition)
    @project_name : répertoire ou seront stockés les fichiers de sortie
    @redondance : si True alors on calcule la distance de Jaccard sur l'ensemble des noeuds, y compris redondants, sinon uniquement sur la présence / absence de noeud
    @color_filename : ce fichier permet de coloriser les feuilles de l'arbre généré
        le fichier doit être au format csv séparé par "," et doit contenir une ligne d'entêtes
        ces entêtes doivent au moins contenir :
            - sample : nom des échantillons qui doivent correspondre aux noms des génômes du pangénome
            - category : nom de catégorie qui sera rajouté au label des feuilles
            - color : code couleur associé à l'échantillon
Returns 
    - dfP : matrice de distance pondérée
    - dfNP : matrice de distance non pondérée
"""
def calculate_distance_jaccard(genome_dic, genome_redondant_dic, project_name, redondance=True, color_file_name=None):
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
    dfP.to_csv(project_name+"/"+"weighted_distance_matrix.csv")
    calcul_arbre_nj(project_name+"/"+"weighted_distance_matrix.csv", project_name, project_name+"/"+"weighted_nj_tree.png", color_file_name)
    dfNP.to_csv(project_name+"/"+"non_weighted_distance_matrix.csv")
    calcul_arbre_nj(project_name+"/"+"non_weighted_distance_matrix.csv", project_name, project_name+"/"+"non_weighted_nj_tree.png", color_file_name)
    print("Distances computed")
    return dfP, dfNP

"""
Fonction principale permettant de charger un graphe GFA et de calculer la distance de Jaccard dessus
On construit également une matrice au format PHYLIP contenant pour chaque genome la présence (1) ou absence (0) pour chaque bnoeud retenu
Parameters
    @filename : nom du fichier GFA contenant le pangénome
    @project_name : va créer un répertoire avec ce nom pour stocker les résultats
    @nb_noeuds_cible : nombre approximatifs de noeuds souhaités dans la matrice d'absence / presence
    @pondere : si True alors calcule la distance de Jaccard pondérée par la taille du noeud, sinon calcule la distance de Jaccard non pondérée
    @methode : "aleatoire" ou "taille", permet de sélectionner les noeuds pour l'analyse (sélection soit aléatoire, soit on sélectionne les plus gros noeuds)
    @redondance : si True alors on va utiliser les noeuds présents en plusieurs exemplaires sur un échantillon, sinon on utilisera juste l'absence ou la présence du noeud
    @strand : si True alors on utilise l'indication du strand, sinon on regarde juste si le noeud est présent sans regarder le sens
    @color_filename : ce fichier permet de coloriser les feuilles de l'arbre généré
        le fichier doit être au format csv séparé par "," et doit contenir une ligne d'entêtes
        ces entêtes doivent au moins contenir :
            - sample : nom des échantillons qui doivent correspondre aux noms des génômes du pangénome
            - category : nom de catégorie qui sera rajouté au label des feuilles
            - color : code couleur associé à l'échantillon
Returns
"""
def analyser_pangenome(file_name, project_name, nb_noeuds_cible = 10000, pondere = True, methode="random", redondance = True, strand = True, color_file_name=None):
    
    print("Launch with args : \nfile_name : " + str(file_name)
          + "\nproject_name : " + str(project_name)
          + "\nnodes number : " + str(nb_noeuds_cible)
          + "\nweighted : " + str(pondere)
          + "\nmethod : " + str(methode)
          +"\nredundancy : " + str(redondance)
          +"\nstrand : " + str(strand)
          +"\ncolor filename : " + str(color_file_name))
    if not os.path.exists(project_name):
        os.mkdir(project_name)
        
    raxml_dir = project_name+"/raxml" 
    if not os.path.exists(raxml_dir):
        os.mkdir(raxml_dir)
    
    matrice_distance_csv_filename = project_name+"/distance_matrix"
    matrice_pav = project_name+"/pav_matrix.phy"
    
    genome_dic, genome_redondant_dic, pav_dic = charger_fichier_gfa(file_name, nb_noeuds_cible, methode, strand)
    
    pav_to_phylip(pav_dic, matrice_pav)
    #Calcul des distances de Jaccard entre les paires
    dfP, dfNP = calculate_distance_jaccard(genome_dic, genome_redondant_dic, project_name, redondance, color_file_name)
    #conversion de la matrice de distance csv vers le format phylip
    distance_matrix_to_phylip(dfP,project_name+"/weighted_distance_matrix.phy")
    distance_matrix_to_phylip(dfNP,project_name+"/non_weighted_distance_matrix.phy")
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
    subprocess.run(raxml_command, cwd="./"+raxml_dir)
    newick_file_name = raxml_dir + "/RAxML_bestTree."+project_name
    plot_newick(newick_file_name, project_name+"/tree_raxml.png", color_file_name)




"""
Fonction de conversion de la matrice PAV sous forme dictionnary vers format phylip

Parameters
    @pav_matrix_dic : matrice PAV au format dictionnary
    @matrice_distance_phylip_filename : nom du fichier qui contiendra la matrice PAV (format PHYLIP)
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
Fonction de conversion de la matrice de distance dataframe vers format phylip

Parameters
    @pav_matrix_dic : matrice PAV au format dictionnary
    @phylip_distance_matrix_file_name : nom du fichier qui contiendra la matrice des distances de Jaccard (format PHYLIP)
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
Fonction permettant de tracer un dendograme
Parameters 
    @newick_tree : arbre à tracer au format newick
    @file_name : nom du fichier de sortie
    @color_filename : ce fichier permet de coloriser les feuilles de l'arbre généré
        le fichier doit être au format csv séparé par "," et doit contenir une ligne d'entêtes
        ces entêtes doivent au moins contenir :
            - sample : nom des échantillons qui doivent correspondre aux noms des génômes du pangénome
            - category : nom de catégorie qui sera rajouté au label des feuilles
            - color : code couleur associé à l'échantillon'
""" 
def plot_newick(newick_tree,file_name,color_filename=None):
    ete3_tree = Tree(newick_tree, format=1)
    if color_filename != None :
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
Fonction permettant de calculer un arbre phylogénétique à partir de la matrice de distance de Jaccard
Paramètres : 
    @matrice_distance_filename : nom du fichier contenant la matrice de distance
Retour :
    - Arbre au format newick
"""   

def calcul_arbre_nj(matrice_distance_filename, project_name, file_name,color_filename=None):
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



#Ancienne fonction permettant de construire un dendogramme depuis une matrice de distance
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



    

