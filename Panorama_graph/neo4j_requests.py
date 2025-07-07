#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 16:03:15 2025

@author: fgraziani
"""

import re
from tqdm import tqdm
from math import *
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import time
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import webbrowser
import os
from typing import List, Dict
import logging
import csv
import json
import subprocess
from config import get_driver
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import pandas as pd
import numpy as np
from Bio import Phylo




#This value is used to limit the nodes number when seeking for regions :
#If the region contains more nodes than this value then it is ignored
max_bp_seeking = 500000

logging.getLogger("neo4j").setLevel(logging.ERROR)




#This function take a region (chromosome, start and stop) of a given haplotype (search_genome)
#and it returns all the nodes in this region and the other related regions : 
#for each haplotype the start and stop are given by the start and stop found on the search_genome
#then the function get all the nodes between this start and stop
def get_nodes_by_region(search_genome, chromosome, start, end ):

    #TODO : les annotations peuvent etre sur plusieurs chromosome => voir comment le traiter
    temps_depart = time.time()
    with get_driver() as driver:

        data = {}
        nodes_data = []
        shared_genomes = []
    
        with driver.session() as session:
    
            # Step 1 : find the region
            if start is not None and end is not None :
                print("Look for region : " + str(start) + " - " + str(end) + " - chromosome " + str(chromosome))
                query_1 = """
                MATCH (n:Noeud)
                WHERE n.chromosome = $chromosome and n.`"""+str(search_genome)+"""_position` >= $start and n.`"""+str(search_genome)+"""_position` <= $end
                RETURN DISTINCT n
                ORDER BY n.`"""+str(search_genome)+"""_noeud`
                """
                result_1 = session.run(query_1, chromosome=chromosome, start=start, end=end)
                find_nodes = [record["n"] for record in result_1]
                print("Number of nodes found : " + str(len(find_nodes)))
                # Setp 2 : shared_genomes construct
                for n in find_nodes:
                    for genome in n["genomes"]:
                        champ_noeud = f"`{genome}_noeud`"
                        position = n.get(f"{genome}_noeud")
        
                        if position is None:
                            continue
        
                        # Check if the genome is already in shared_genomes
                        existe = False
                        for g in shared_genomes:
                            if g["genome"] == genome:
                                g["end"] = max(g["end"], position)
                                existe = True
                                break
                        if not existe:
                            shared_genomes.append({
                                "genome": genome,
                                "start": position,
                                "end": position
                        })

                nodes_data = {}
                # Step 3 : Get the nodes and annotations for each genomes in shared_genomes
                for g in shared_genomes:
                    genome = g["genome"]
                    champ_noeud = f"`{genome}_noeud`"
                    start = g["start"]
                    end = g["end"]
                    if end - start < max_bp_seeking and end - start > 0:
                        print("Search nodes for genome : " + str(g))
                        query_genome = f"""
                        MATCH (m:Noeud)
                        WHERE  m.chromosome = $chromosome and m.{champ_noeud} >= $start AND m.{champ_noeud} <= $end 
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        RETURN m, collect(a.gene_name) as annotations
                        ORDER BY m.{champ_noeud}
                            """
                        result_genome = session.run(query_genome, chromosome=chromosome, start=start, end=end)


                        for record in result_genome:
                            if record["m"]["name"] not in nodes_data :
                                nodes_data[record["m"]["name"]] = dict(record["m"]) |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))}
                            else:
                                for a in range(len(record["annotations"])):
                                    nodes_data[record["m"]["name"]]["annotations"].add(record["annotations"][a])
                    else :
                        if end - start >= max_bp_seeking:
                            print("Region too wide for genome : " + str(g) + " start : " + str(start) + " end : " + str(end))
                for elt in nodes_data:
                    nodes_data[elt]["annotations"] = list(nodes_data[elt]["annotations"])
            print("searching annotations")
        print("Total time : " + str(time.time() - temps_depart))
        return nodes_data





def get_nodes_by_gene(genome_ref, gene_id=None, gene_name=None, chromosome = None):
    with get_driver() as driver:

        nodes_data = {}
        shared_genomes = []
    
        with driver.session() as session:
    
            # Step 1 : find nodes with gene annotation
            if gene_name is not None :
                print("Looking for gene name : " + str(gene_name))
                if chromosome is None or chromosome == "" :
                    query_1 = """
                    MATCH (a:Annotation {gene_name: $gene_name})<-[:A_POUR_ANNOTATION]-(n:Noeud)
                    RETURN DISTINCT n
                    ORDER BY n.`"""+str(genome_ref)+"""_noeud`
                    """
                    result_1 = session.run(query_1, gene_name=gene_name)
                else:
                    query_1 = """
                    MATCH (a:Annotation {chromosome:$chromosome, gene_name: $gene_name})<-[:A_POUR_ANNOTATION]-(n:Noeud)
                    RETURN DISTINCT n
                    ORDER BY n.`"""+str(genome_ref)+"""_noeud`
                    """
                    result_1 = session.run(query_1, chromosome=chromosome, gene_name=gene_name)
                
            else:
                print("Looking for gene id : " + str(gene_id))
                if chromosome is None or chromosome == "" :
                    query_1 = """
                    MATCH (a:Annotation {gene_id: $gene_id})<-[:A_POUR_ANNOTATION]-(n:Noeud)
                    RETURN DISTINCT n
                    ORDER BY n.`"""+str(genome_ref)+"""_noeud`
                    """
                    result_1 = session.run(query_1, gene_id=gene_id)
                    
                else:
                    query_1 = """
                    MATCH (a:Annotation {chromosome:$chromosome, gene_id: $gene_id})<-[:A_POUR_ANNOTATION]-(n:Noeud)
                    RETURN DISTINCT n
                    ORDER BY n.`"""+str(genome_ref)+"""_noeud`
                    """
                    result_1 = session.run(query_1, chromosome=chromosome, gene_id=gene_id)
            noeuds_annotes = [record["n"] for record in result_1]
    
            #Select first chromosome
            if chromosome is None :
                chromosome = str(noeuds_annotes[0]["chromosome"])
                print("chromosome : " + str(chromosome))
            if len(noeuds_annotes) == 0:
                print("No nodes found")
            else:
                # Step 2 : construction of shared_genomes
                for n in noeuds_annotes:
                    for genome in n["genomes"]:
                        champ_position = f"`{genome}_position`"
                        position = n.get(f"{genome}_position")
        
                        if position is None:
                            continue
        
                        # Search if this genome already exists in shared_genomes
                        existe = False
                        for g in shared_genomes:
                            if g["genome"] == genome:
                                g["end"] = max(g["end"], position)
                                existe = True
                                break
                        if not existe:
                            shared_genomes.append({
                                "genome": genome,
                                "start": position,
                                "end": position
                            })
                # Step 3 : for each genome, get nodes of the found region
                start = 0
                end = 0
                for g in shared_genomes:
                    genome = g["genome"]
                    champ_position = f"`{genome}_position`"
                    start = g["start"]
                    end = g["end"]
                    if g["genome"] == genome_ref:
                        start = start
                        end = end
                    print("Search for genome : " + str(g))
                    if end - start < max_bp_seeking and end - start > 0:
                        query_genome = f"""
                        MATCH (m:Noeud)
                        WHERE m.chromosome = $chromosome and m.{champ_position} >= $start AND m.{champ_position} <= $end
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        RETURN m, collect(a.gene_name) as annotations
                        ORDER BY m.{champ_position}
                        """
                        result_genome = session.run(query_genome, chromosome=chromosome, start=start, end=end)
                        for record in result_genome:
                            if record["m"]["name"] not in nodes_data :
                                nodes_data[record["m"]["name"]] = dict(record["m"]) |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))}
                            else:
                                for a in range(len(record["annotations"])):
                                    nodes_data[record["m"]["name"]]["annotations"].add(record["annotations"][a])
                    else :
                        if end - start >= max_bp_seeking:
                            print("Region too wide for genome : " + str(g) + " start : " + str(start) + " end : " + str(end))
            for elt in nodes_data:
                nodes_data[elt]["annotations"] = list(nodes_data[elt]["annotations"])
        return nodes_data



#This function get first annotation on nodes of a reference_genome on a given chromosome after the given position
def get_first_annotation_after_position(genome_ref, chromosome="1", position=0):
    with get_driver() as driver:

        query = """
        MATCH (n:Noeud)-[:A_POUR_ANNOTATION]->(a:Annotation)
        WHERE n.chromosome = $chromosome and n."""+str(genome_ref)+"""_position > $position
        WITH n, a
        ORDER BY n.HER_position ASC
        LIMIT 1
        RETURN a.gene_id AS gene_id, a.gene_name AS gene_name, a.start AS start, a.end AS end
        """  
        
        
        with driver.session() as session:
            result = session.run(query, chromosome=chromosome, position=position)
            record = result.single()
            return dict(record) if record else None


#This function get all annotations on nodes of a reference_genome on a given chromosome between start_position and end_position
def get_annotations_in_position_range(genome_ref, chromosome="1", start_position=0, end_position=0):
    with get_driver() as driver:

        query = """
        MATCH (n:Noeud)-[:A_POUR_ANNOTATION]->(a:Annotation)
        WHERE n.chromosome = $chromosome and n.`"""+str(genome_ref)+"""_position` >= $start AND n.`"""+str(genome_ref)+"""_position` <= $end
        RETURN DISTINCT a.gene_id AS gene_id, a.gene_name AS gene_name
        """
        #print("Query : " + query)
        with driver.session() as session:
            result = session.run(query, chromosome=chromosome, start=start_position, end=end_position)
            annotations = [dict(record) for record in result]
    
    return annotations

#This function will get all chromosomes present in the pangenome graph
def get_chromosomes():
    with get_driver() as driver:
    
        query = """
        MATCH (s:Stats) 
        RETURN s.chromosomes as all_chromosomes
        """
        all_chromosomes = []
        with driver.session() as session:
            result = session.run(query)
            for record in result:
                all_chromosomes = record["all_chromosomes"]
    
    return all_chromosomes

#This function will get all genomes present in the pangenome graph
def get_genomes():
    with get_driver() as driver:
    
        query = """
        MATCH (s:Stats) 
        RETURN s.genomes as all_genomes
        """
        all_genomes = []
        with driver.session() as session:
            result = session.run(query)
            for record in result:
                all_genomes = record["all_genomes"]
    
    return all_genomes



def analyse_to_csv(analyse, output_file):
    with open(output_file, mode='w', newline='', encoding='utf-8') as file:
        genome_ref = list(analyse.keys())[0]
        file.write(genome_ref+"\n")
        fields = analyse[genome_ref][0]
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(analyse[genome_ref])


#Function get_shared_regions : this function is used to get shared regions (nodes) between a list of genomes (GWAS)
#genomes_list : list of the genomes for which the function will look for shared regions
#genome_ref will be used to get annotations on this genome
#chromosomes : if defined the function will only look for shared region on theses chromosomes
#node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
#nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
#deletion : if True the function will look for nodes where no one of the genome set is present
#region_trim : the shared region will be expanded in order to visualise a small region before and after. Set to 0 if strict shared regions are desired.
def get_shared_regions(genomes_list, genome_ref=None, chromosomes=None, node_min_size = 10, nodes_max_gap = 100, deletion=False, region_trim = 1000):
    dic_regions, analyse = find_shared_regions(genomes_list, genome_ref, chromosomes, node_min_size, nodes_max_gap, deletion = deletion)
    shared_regions_dict = {}
    annotations_by_regions = {}
    if genome_ref in dic_regions :
        g = genome_ref
    else :
        g = list(dic_regions.keys())[0]
    for c,items in dic_regions[g].items():
        for r in items["regions"]:
            shared_regions_dict[c+"-"+str(r["start"]-region_trim) + "-"+str(r["stop"]+region_trim)] = get_nodes_by_region(g, c, r["start"]-region_trim, r["stop"]+region_trim)
            
    return shared_regions_dict, analyse
            

def find_first_ref_node_node(genome, genome_ref, genome_position, type_search = "before", chromosome="1"):
    with get_driver() as driver:
        if type_search == "before" :

            query = """
            MATCH (n:Noeud)
            WHERE n.chromosome = $chromosome and n.`"""+str(genome)+"""_position` <= $genome_position AND $genome_ref in n.genomes
            return max(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        else:

            query = """
            MATCH (n:Noeud)
            WHERE n.chromosome = $chromosome and n.`"""+str(genome)+"""_position` >= $genome_position AND $genome_ref in n.genomes
            return min(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        #print("Query : " + query)
        with driver.session() as session:
            result = session.run(query, chromosome=chromosome, genome_ref=genome_ref, genome_position=genome_position)
            for record in result :
                ref_position = record["ref_position"]
    
    return ref_position

#Function get_shared_regions : this function is used to get shared regions (positions) between a list of genomes (GWAS)
#On peut préciser sur quel chromosome chercher
#Usage exemple (cattle white spot) : dic_regions, analyse = find_shared_regions(["HER","SIM"],chromosome="6")
#genomes_list : list of the genomes for which the function will look for shared regions
#genome_ref will be used to get annotations on this genome
#chromosomes : if defined the function will only look for shared region on these chromosomes
#node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
#nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
def find_shared_regions(genomes_list, genome_ref=None, chromosomes=None, node_min_size = 10, nodes_max_gap = 100, deletion=False):
    dic_regions = {}
    temps_depart = time.time()
    if (len(genomes_list) > 1):
        print("finding shared regions for " + str(genomes_list))
        with get_driver() as driver:
            with driver.session() as session:
                query = """
                MATCH (s:Stats)
                RETURN s.genomes AS genomes
                LIMIT 1
                """
                result = session.run(query)
                for record in result:
                    genomes = record["genomes"]
                
                nb_genomes = len(genomes)
                nb_regions_total = 0
                nb_associated_genomes = len(genomes_list)
                
                flux_max = nb_associated_genomes / nb_genomes
                flux_max += 0.1 * flux_max
                print("genomes number : " + str(nb_genomes))
                
                if chromosomes != None :
                    chromosome_list = chromosomes
                else :
                    chromosome_list = get_chromosomes()

                    

                for c in chromosome_list :
                    dic_regions[c] = {}
                    query = 'MATCH (n:Noeud) USING INDEX n:Noeud(chromosome) where n.chromosome = "' + str(c) + '" AND n.flux <= ' + str(flux_max)
                    dic_number_to_position = {}  
                    dic_number_to_taille = {}
                    for g in genomes_list :
                        dic_regions[c][g] = {"nodes_number_list":[], "regions" : []}
                        dic_number_to_position[g] = {}
                        dic_number_to_taille[g] = {}
                        query += ' AND "' + str(g) + '" IN n.genomes'

                    query += " WITH n where n.taille > " + str(node_min_size)
                    query += " RETURN n AS noeuds order by n.HER_noeud ASC"
   
                    result1 = list(session.run(query))
                    if deletion :

                        query = "MATCH (n:Noeud)"
                        query += """
                            WHERE ALL(g IN $genomes_list WHERE g IN n.genomes)
                              AND size(n.genomes) > size($genomes_list)
                            
                            WITH n, [g IN n.genomes WHERE NOT g IN $genomes_list] AS autres_genomes
                            MATCH (n)-[]->(m:Noeud)
                            WHERE ALL(g IN autres_genomes WHERE g IN m.genomes)
                              AND NONE(g IN $genomes_list WHERE g IN m.genomes)
                        """
                        query += " WITH n where n.taille > " + str(node_min_size)
                        query += " RETURN n AS noeuds order by n.HER_noeud ASC"

                        result = result1 + list(session.run(query, genomes_list=genomes_list))
                    else:
                        result = result1
                    
    
                    print("Results number for chromosome " + str(c) + " : " + str(len(result)))
                    nb_regions_total += len(result)
                    for r in result:
                        for g in genomes_list:
                            dic_regions[c][g]["nodes_number_list"].append(r["noeuds"][g+"_noeud"]) 
                            dic_number_to_position[g][r["noeuds"][g+"_noeud"]] = r["noeuds"][g+"_position"]
                            dic_number_to_taille[g][r["noeuds"][g+"_noeud"]] = r["noeuds"]["taille"]

                    #Group regions if they are separated vy less than nodes_max_gap
                    for g in genomes_list :
                        dic_regions[c][g]["nodes_number_list"].sort()
                        for i in range(len(dic_regions[c][g]["nodes_number_list"])):
                            if i == 0 :
                                region_start = dic_number_to_position[g][dic_regions[c][g]["nodes_number_list"][0]]
                                region_stop = region_start + dic_number_to_taille[g][dic_regions[c][g]["nodes_number_list"][0]]
                            else :
                                if dic_regions[c][g]["nodes_number_list"][i] < dic_regions[c][g]["nodes_number_list"][i-1] + nodes_max_gap :
                                    region_stop = dic_number_to_position[g][dic_regions[c][g]["nodes_number_list"][i]]
                                else :
                                    dic_regions[c][g]["regions"].append({"start" : region_start, "stop" : region_stop, "taille" : region_stop - region_start})
                                    region_start = dic_number_to_position[g][dic_regions[c][g]["nodes_number_list"][i]]
                                    region_stop = region_start + dic_number_to_taille[g][dic_regions[c][g]["nodes_number_list"][i]]
            print("Total number of regions : " +str(nb_regions_total))
            dic_regions_2 = {}
            for c, genomes in dic_regions.items():
                for genome, valeur in genomes.items():
                    if genome not in dic_regions_2:
                        dic_regions_2[genome] = {}
                    dic_regions_2[genome][c] = valeur
            
            #print("Genomes : " + str(genomes))
            analyse = {}
            for g in tqdm(genomes):
                analyse[g] = []
                for c in tqdm(dic_regions_2[g]):
                    for r in tqdm(dic_regions_2[g][c]['regions']):
                        r["chromosome"] = c
                        r["genome"] = g
                        if genome_ref is not None and g == genome_ref:
                            r["annotations"] = get_annotations_in_position_range(genome_ref=g,chromosome=c, start_position=r["start"],end_position=r["stop"])
                            r["first_annotation_after_region"] = get_first_annotation_after_position(genome_ref=g,chromosome=c, position=r["stop"])
                        analyse[g].append(r)
                analyse[g] = sorted(analyse[g], key=lambda d: d['taille'], reverse=True)
            
            
            print("Total time : "+ str(time.time()-temps_depart))
    return dic_regions_2, analyse



def calculer_variabilite(chromosome_list=None, ref_genome=None, window_size=1000, output_html="pangenome_variability.html"):
    with get_driver() as driver:
        if ref_genome == None :
            ref_noeud = "moyenne_noeud"
            ref_position = "moyenne_position"
        else :
            ref_noeud = f"{ref_genome}_noeud"
            ref_position = f"{ref_genome}_position"
        # Définir le renderer pour ouvrir dans le navigateur
        pio.renderers.default = 'browser'
        
        html_parts = []
        if chromosome_list is None:
            chromosomes = get_chromosomes()
        else:
            chromosomes = chromosome_list
  
        for chromosome in chromosomes:
            print("Compute variability on chromosome " + str(chromosome))
            query = f"""
            MATCH (n:Noeud)
            WHERE n.chromosome = $chromosome
            WITH FLOOR(n.{ref_position} / {window_size}) AS Window, n.{ref_position} AS position, n.flux AS flux
            WITH Window, 
                 collect(position) AS positions, 
                 avg(flux) AS moyenne_flux, 
                 count(*) AS nodes_number
            RETURN 
                Window,
                reduce(min_pos = head(positions), p IN positions | CASE WHEN p < min_pos THEN p ELSE min_pos END) AS start_position,
                reduce(max_pos = head(positions), p IN positions | CASE WHEN p > max_pos THEN p ELSE max_pos END) AS end_position,
                moyenne_flux,
                nodes_number
            ORDER BY Window
            """
        
            with driver.session() as session:
                result = session.run(query, chromosome=chromosome)
                data = result.data()
        
            df = pd.DataFrame(data)
            if df.empty:
                print(f"No data found for chromosome {chromosome}.")
                continue
        
            # Add window size and hover text
            df['window_size'] = df['end_position'] - df['start_position']
            df['hover_text'] = (
                "Window: " + df['Window'].astype(str) +
                "<br>Start: " + df['start_position'].astype(str) +
                "<br>End: " + df['end_position'].astype(str) +
                "<br>Size: " + df['window_size'].astype(str) +
                "<br>Moyenne flux: " + df['moyenne_flux'].round(2).astype(str)
            )
        
            # Create figure
            fig = px.scatter(
                df,
                x='Window',
                y='moyenne_flux',
                hover_name='hover_text',
                labels={
                    'Window': 'Window',
                    'moyenne_flux': 'Moyenne du flux'
                },
                title=f"Flow variability - Chromosome {chromosome} ({ref_genome}) - window size {window_size}",
                height=500
            )
        
            fig.update_traces(marker=dict(size=10, color='dodgerblue', line=dict(width=1, color='darkblue')))
            fig.update_layout(hovermode='closest')
        
            # Generate HTML of the figure and store it
            html_parts.append(pio.to_html(fig, include_plotlyjs='cdn', full_html=False))
        
        if not html_parts:
            print("No graph computed.")
            return
        
        # Combine all graphs into a single HTML file
        with open(output_html, 'w') as f:
            f.write("<html><head><title>Variabilité flux</title></head><body>\n")
            for part in html_parts:
                f.write(part + "<hr>\n")
            f.write("</body></html>")
        
        print(f"Graphs stored into {output_html}")
        file_path = os.path.abspath(output_html)
        webbrowser.open(f'file://{file_path}')

    
#This function takes nodes data (from get_nodes_by_region ou get_nodes_by_gene for exemple)
#it computes the Jaccard distance on these nodes
#and it returns the distance matrix and the distance matrix weighted by the nodes size
def compute_phylo_tree_from_nodes(nodes_data,output_dir = "", weighted=False):
    genomes = get_genomes()
    nodes_nb = len(nodes_data)
    genome_redondant_strand_matrix = np.zeros((len(genomes),2 * nodes_nb), dtype=np.int32)
    genome_strand_matrix = np.zeros((len(genomes),2 * nodes_nb), dtype=np.int32)
    
    i = 0
    index_genomes = {}
    genomes_names = []
    for g in genomes:
        index_genomes[g] = i
        genomes_names.append(g)
        i += 1
    i = 0

    for n in nodes_data:
        for g in genomes :
            if g in nodes_data[n]["strandM"]:
                genome_redondant_strand_matrix[index_genomes[g], i] += 1
                genome_strand_matrix[index_genomes[g], i] = nodes_data[n]["taille"]
            if g in nodes_data[n]["strandP"]:
                genome_redondant_strand_matrix[index_genomes[g], i+nodes_nb] += 1
                genome_strand_matrix[index_genomes[g], i+nodes_nb] = nodes_data[n]["taille"]
        i += 1
    #computes Jaccard distance on matrix
    jaccard_matrix = np.zeros((len(genomes), len(genomes)))
    weighted_jaccard_matrix = np.zeros((len(genomes), len(genomes)))
    for i in range(len(genomes)):
        for j in range(i, len(genomes)):
            min_sum = np.minimum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :]).sum()
            max_sum = np.maximum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :]).sum()
            
            min_counts = np.minimum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :])
            max_counts = np.maximum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :])
            
            size_matrix = np.maximum(genome_strand_matrix[i, :], genome_strand_matrix[j, :])
            # Weight by nodes size
            weighted_min = (min_counts * size_matrix).sum()
            weighted_max = (max_counts * size_matrix).sum()
            
            if max_sum == 0:
                jaccard_index = 0.0
            else:
                jaccard_index = min_sum / max_sum
                
            if weighted_max == 0:
                weighted_jaccard_index = 0.0
            else:
                weighted_jaccard_index = weighted_min / weighted_max

            jaccard_matrix[i, j] = 1-jaccard_index
            jaccard_matrix[j, i] = 1-jaccard_index
            
            weighted_jaccard_matrix[i,j] = 1-weighted_jaccard_index
            weighted_jaccard_matrix[j,i] = 1-weighted_jaccard_index
            
    df_jaccard = pd.DataFrame(jaccard_matrix, index=genomes_names, columns=genomes_names)
    df_weighted_jaccard = pd.DataFrame(weighted_jaccard_matrix, index=genomes_names, columns=genomes_names)
    if output_dir != "" and os.path.isdir(output_dir) :
        df_jaccard.to_csv(output_dir+'/distance_matrix.csv')
        df_weighted_jaccard.to_csv(output_dir+'/weighted_distance_matrix.csv')
    
    if weighted:
        df_val = df_weighted_jaccard.values
    else:
        df_val = df_jaccard.values
    triangulaire_inf = []
    for i in range(0, df_val.shape[0]):
        triangulaire_inf.append([])
        for j in range(i + 1):
            triangulaire_inf[i].append(df_val[i, j])

    dm = Phylo.TreeConstruction.DistanceMatrix(list(df_jaccard.columns), triangulaire_inf)
    # computing tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    
    
    newick_tree = tree.format('newick')
    print(newick_tree)
    return newick_tree
            


#Computes the distance matrix on the whole GFA (could take a long time for big GFA)
def compute_distance_matrix(distance_matrix_filename = "distances.csv", chromosome=None, ponderation=True, strand=False):
    temps_depart = time.time()
    with get_driver() as driver:
        with driver.session() as session:
            #Get genomes list
            query = """
            MATCH (s:Stats)
            RETURN s.genomes AS genomes
            LIMIT 1
            """
            result = session.run(query)
            for record in result:
                genomes = record["genomes"]

            distance_matrix = pd.DataFrame(data=0,index=genomes, columns=genomes)
            dic_taille_genome = {}   
            if ponderation :
                #Get the genome size
                query = """
                    MATCH (n:Noeud)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, sum(n.taille) AS total_taille
                    """
                result_taille_genomes = list(session.run(query))
                if strand == False:
                    #Calcul des tailles par genome

                    #Get intersection size
                    query = """
                        MATCH (n:Noeud)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, sum(n.taille) AS taille_intersection
                        ORDER BY taille_intersection DESC
                        """

                else:
                    query = """
                        MATCH (n:Noeud)
                        WITH n.taille AS taille, n.strandM AS strandM, n.strandP AS strandP
                        WITH taille,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, taille
                        WHERE g1 < g2
                        RETURN g1, g2, sum(taille) AS taille_intersection
                        ORDER BY taille_intersection DESC
                        """
                result_intersection = list(session.run(query))
                    
            else:
                #Get the node numbers
                query = """
                    MATCH (n:Noeud)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, count(*) AS total_taille
                """
                result_taille_genomes = list(session.run(query)) 
                if strand == False:

                    #Get intersection size
                    query = """
                        MATCH (n:Noeud)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, count(*) AS taille_intersection
                        ORDER BY taille_intersection DESC
                        """
                else:
                    query = """
                        MATCH (n:Noeud)
                        WITH n.taille AS taille, n.strandM AS strandM, n.strandP AS strandP
                        WITH taille,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, taille
                        WHERE g1 < g2
                        RETURN g1, g2, count(*) AS taille_intersection
                        ORDER BY taille_intersection DESC
                        """
                result_intersection = list(session.run(query))
            
        for r in result_taille_genomes:
            dic_taille_genome[r["genome"]] = r["total_taille"]
        print(dic_taille_genome)
        for r in result_intersection:
            g1 = r["g1"]
            g2 = r["g2"]
            inter = r["taille_intersection"]
            distance_matrix.loc[g1,g2] = 1-inter/(dic_taille_genome[g1]+dic_taille_genome[g2]-inter)
            distance_matrix.loc[g2,g1] = distance_matrix.loc[g1,g2]
                
                
        distance_matrix.to_csv(distance_matrix_filename)   
        print("Total time : "+ str(time.time()-temps_depart))
    return distance_matrix
          
 