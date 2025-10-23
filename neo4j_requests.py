#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 16:03:15 2025

@author: fgraziani
"""

import re
from tqdm import tqdm
from math import *
from numpy.polynomial.polynomial import Polynomial
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import time
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import webbrowser
import os
import glob
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
from Bio.Seq import Seq
import statistics




#This value is used to limit the nodes number when seeking for regions :
#If the region is wider than this value then it is ignored
max_bp_seeking = 800000

#Maximal number of nodes to get the whole region
max_nodes_number = 100000

logging.getLogger("neo4j").setLevel(logging.ERROR)



def get_anchor(genome, chromosome, position, before = True, use_anchor=True):
    if get_driver() is None :
        return None
    genome_position = genome+"_position"
    if use_anchor :
        window_size=1000
        max_attemps = 500
        attempt = 0
        with get_driver() as driver:
            with driver.session() as session:
                while attempt < max_attemps :
                    attempt += 1
            
                    if before:
                        if attempt == 1:
                            lower_bound = position
                        upper_bound = lower_bound
                        lower_bound = lower_bound - window_size
                        order = "DESC"
    
                    else:
                        if attempt == 1:
                            upper_bound = position
                        lower_bound = upper_bound
                        upper_bound = upper_bound + window_size
                        order = "ASC"
            
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= $lower_bound
                      AND n.`{genome_position}` <= $upper_bound
                      AND n.flow >= 1.0
                    RETURN n
                    ORDER BY n.`{genome_position}` {order}
                    LIMIT 1
                    """
    
                    result = session.run(
                            query,
                            chromosome=chromosome,
                            lower_bound=lower_bound,
                            upper_bound=upper_bound
                        )
                    record = result.single()
                    if record:
                        return dict(record["n"])
                #Anchor not found => the current node will be used
                if before:
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= $position
                    RETURN n order by n.`{genome_position}` ASC limit 1
                    """
                else :
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` <= $position
                    RETURN n order by n.`{genome_position}` DESC limit 1
                    """
    
                result = session.run(
                        query,
                        chromosome=chromosome,
                        position=position
                    )
                record = result.single()
                if record:
                    return dict(record["n"])
    else:
        if before:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` >= $position
            RETURN n order by n.`{genome_position}` ASC limit 1
            """
        else :
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` <= $position
            RETURN n order by n.`{genome_position}` DESC limit 1
            """
        with get_driver() as driver:
            with driver.session() as session:
                result = session.run(
                        query,
                        chromosome=chromosome,
                        position=position
                    )
                record = result.single()
                if record:
                    return dict(record["n"])    

    return None



#This function take a region (chromosome, start and stop) of a given haplotype (search_genome)
#and it returns all the nodes in this region and the other related regions : 
#for each haplotype the start and stop are given by the first anchor node before the start position and the first after the end position
#Anchor node is a node with all haplotypes
def get_nodes_by_region(genome, chromosome, start, end, use_anchor = True ):
    if get_driver() is None :
        return []
    temps_depart = time.time()
    nodes_data = {}
    max_sequence = 1000
    if end is None:
        stop = max_bp_seeking
    else:
        stop = end
    with get_driver() as driver:

        data = {}
        shared_genomes = []
        with driver.session() as session:
    
            # Step 1 : find the region
            if start is not None and end is not None :
                genome_position = genome+"_position"
                print("Look for region : " + str(start) + " - " + str(stop) + " - chromosome " + str(chromosome) + " - genome : " + str(genome))

                anchor_start = get_anchor(genome, chromosome, start, before = True, use_anchor=use_anchor)
                anchor_stop = get_anchor(genome, chromosome, end, before = False, use_anchor=use_anchor)
                delta = abs(stop - start)
                if anchor_start is not None and  anchor_stop is not None:
                    if anchor_start[genome_position] > anchor_stop[genome_position]:
                        anchor_start_tmp = anchor_start
                        anchor_start = anchor_stop
                        anchor_stop = anchor_start_tmp
                    print("Anchor start name : " + str(anchor_start["name"]))
                    print("Anchor stop name : " + str(anchor_stop["name"]))
                    print("Anchor region : " + str(anchor_start[genome_position]) + " - " + str(anchor_stop[genome_position]))
                
                if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - anchor_start[genome_position] < max_bp_seeking and anchor_stop[genome_position] - anchor_start[genome_position] > 0 and len(anchor_start['genomes']) > 0 :

                    # Step 3 : Get the nodes and annotations for each genomes
                    query_genome = f"""
                    MATCH (m:Node)
                    WHERE  m.chromosome = "{chromosome}" 
                    AND (
                    """
                    first = True
                    for g in set(anchor_start["genomes"]+anchor_stop["genomes"]):
                        position_field = g+"_position"
                        if position_field in anchor_start and position_field in anchor_stop:
                            start = min(anchor_start[position_field],anchor_stop[position_field])
                            stop = max(anchor_start[position_field],anchor_stop[position_field])
                            if first :
                                query_genome += f"(m.`{position_field}` >= {start} AND m.`{position_field}` <= {stop})"
                                first = False
                            else :
                                query_genome += f" OR (m.`{position_field}` >= {start} AND m.`{position_field}` <= {stop})"
                    
                    query_genome += f""")
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
                        RETURN m, substring(s.sequence, 0, {max_sequence}) as sequence, collect(a.gene_name) AS annotations, collect(a.feature) AS features
                        """

                    
                    
                    #print(query_genome)
                    result = session.run(query_genome, start=start, end=end)
                    for record in result :
                        nodes_data[record["m"]["name"]] = dict(record["m"]) |{"sequence":record["sequence"]} |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))} |{"features":set(record["features"][a] for a in range(len(record["features"])))}
                else:
                    if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - anchor_start[genome_position] >= max_bp_seeking :
                        print(f"Region too wide : {anchor_stop[genome_position] - anchor_start[genome_position]}" )
                    else:
                        print("Region not found")
                    nodes_data = {}

            else:
                #if end is set to None, get all the nodes if the number of nodes is less than max_nodes_number
                
                if start == 0 and end is None :
                    total_nodes = get_nodes_number(chromosome)
                    if total_nodes <= max_nodes_number : 
                        query_genome = f"""
                        MATCH (m:Node)
                        WHERE  m.chromosome = "{chromosome}" 
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
                        RETURN m, substring(s.sequence, 0, {max_sequence}) as sequence, collect(a.gene_name) AS annotations, collect(a.feature) AS features
                            """
                        result = session.run(query_genome)
                        for record in result:
                            nodes_data[record["m"]["name"]] = dict(record["m"]) |{"sequence":record["sequence"]}  |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))} |{"features":set(record["features"][a] for a in range(len(record["features"])))}
                    else  :
                        print("Region too wide")
                        nodes_data = {}
            if len(nodes_data) > 0 :
                for elt in nodes_data:
                    nodes_data[elt]["annotations"] = list(nodes_data[elt]["annotations"])
                    nodes_data[elt]["features"] = list(nodes_data[elt]["features"])
        print("Total time : " + str(time.time() - temps_depart))
        return nodes_data




def get_nodes_by_gene(genome, chromosome, gene_id=None, gene_name=None):
    if get_driver() is None :
        return []
    with get_driver() as driver:

        nodes_data = {}
        shared_genomes = []
    
        with driver.session() as session:
            genome_position = genome+"_position"
            # Step 1 : find nodes with gene annotation
            if gene_name is not None :
                print("Looking for gene name : " + str(gene_name))
                query = f"""
                MATCH (a:Annotation {{chromosome:"{chromosome}", gene_name: $gene_name}})<-[]-(n:Node)
                RETURN DISTINCT n
                ORDER BY n.`{genome_position}` ASC
                """
                result = session.run(query, gene_name=gene_name)
                #print(query)
            else:
                query = f"""
                MATCH (a:Annotation {{chromosome:"{chromosome}", gene_id: $gene_id}})<-[]-(n:Node)
                RETURN DISTINCT n
                ORDER BY n.`{genome_position}` ASC
                """
                result = session.run(query, gene_id=gene_id)
            noeuds_annotes = [record["n"] for record in result]
            if len(noeuds_annotes) > 0 and genome_position in noeuds_annotes[0] and genome_position in noeuds_annotes[-1]:
                start = noeuds_annotes[0][genome_position]
                stop = noeuds_annotes[-1][genome_position] + noeuds_annotes[-1]["size"]
                print(f"start : {start} - stop : {stop} - nodes number : {len(noeuds_annotes)}")
                nodes_data = get_nodes_by_region(genome, chromosome, start, stop)
            
        return nodes_data



#This function get first annotation on nodes of a reference_genome on a given chromosome after the given position
def get_annotation_before_or_after_position(genome_ref, chromosome="1", position=0, before=True):
    if get_driver() is None :
        return {}
    with get_driver() as driver:
        if before :
            query = f"""
            MATCH (a:Annotation)
            WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.end < $position and a.gene_name is not null
            ORDER BY a.end DESC
            LIMIT 1
            RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
            """  
        else:
            query = f"""
            MATCH (a:Annotation)
            WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.start > $position and a.gene_name is not null
            ORDER BY a.start ASC
            LIMIT 1
            RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
            """  
        
        with driver.session() as session:
            result = session.run(query, position=position)
            record = result.single()
            return dict(record) if record else None


#This function get all annotations on nodes of a reference_genome on a given chromosome between start_position and end_position
def get_annotations_in_position_range(genome_ref, chromosome="1", start_position=0, end_position=0):
    if get_driver() is None :
        return []
    with get_driver() as driver:
        query = f"""
        MATCH (a:Annotation)
        WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.end <= {end_position} and a.start >= {start_position} and a.gene_name is not null
        RETURN DISTINCT a.gene_id AS gene_id, a.gene_name AS gene_name, a.feature as feature
        """ 
        
        # query = f"""
        # MATCH (n:Node)-[]->(a:Annotation)
        # WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome_ref)+"""_position` >= $start AND n.`"""+str(genome_ref)+"""_position` <= $end and a.gene_name is not null and a.genome_ref = $genome_ref 
        # RETURN DISTINCT a.gene_id AS gene_id, a.gene_name AS gene_name, a.feature as feature
        # """
        #print("Query : " + query)
        with driver.session() as session:
            result = session.run(query, start=start_position, end=end_position, genome_ref=genome_ref)
            annotations = [dict(record) for record in result]
    
    return annotations

#This function will get all chromosomes present in the pangenome graph
def get_chromosomes():
    if get_driver() is None :
        return []
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


#This function will get the number of nodes in the graph
def get_nodes_number(chromosome=None):
    if get_driver() is None :
        return 0
    total = 0
    with get_driver() as driver:
        if chromosome is None and chromosome != "":
            query = """
            MATCH (n:Node) 
            RETURN count(n) as total
            """
        else:
            query = f"""
            MATCH (n:Node) 
            where n.chromosome="{chromosome}"
            RETURN count(n) as total
            """
        all_chromosomes = []
        with driver.session() as session:
            result = session.run(query, chromosome=chromosome)
            for record in result:
                total = record["total"]
    
    return total

#This function will get all genomes present in the pangenome graph
def get_genomes():
    if get_driver() is None :
        return []
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

#This function get a sequence from a list of nodes names list
def get_sequence_from_names(names):
    if get_driver() is None :
        return None
    if len(names) == 0 :
        return {}
    else:
        with get_driver() as driver:
        
            query = """
            MATCH (s:Sequence) 
            WHERE s.name IN $names
            RETURN s.name as name, s.sequence as sequence
            """
            
            with driver.session() as session:
                result = session.run(query, names=names)
                return {record["name"]: record["sequence"] for record in result}
    

#This function get a sequence from a start - end / chromosome position for a given genome
def get_sequence_from_position(genome, chromosome, start, end):
    if get_driver() is None :
        return None
    sequence = ""
    if genome is None or genome == "" or chromosome is None or chromosome == "" or start is None or end is None :
        return None
    else:
        with get_driver() as driver:
            position_key = genome+"_position"
            query =  f"""MATCH (n:Node)
             WHERE n.chromosome = "{chromosome}"
               AND n.`{position_key}` >= {start}
               AND n.`{position_key}` <= {end}
               return n.ref_node AS name, 
               coalesce(n.strandM, []) AS strandMList,
               "{genome}" IN coalesce(n.strandM, []) AS strandM
               order by n.`{position_key}` ASC
            """
            #print(query)
            with driver.session() as session:
                result = session.run(query)
                sorted_names = []
                sorted_strandM = []
                for record in result:
                    sorted_names.append(record["name"])
                    sorted_strandM.append(record["strandM"])
            sequences = get_sequence_from_names(sorted_names)
            for i in range(len(sorted_names)):
                if sorted_strandM[i]:
                    sequence += Seq(sequences[sorted_names[i]]).reverse_complement()
                else:
                    sequence += sequences[sorted_names[i]]
    return sequence
    


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
def get_shared_regions(genomes_list, genome_ref=None, chromosomes=None, node_min_size = 10, nodes_max_gap = 100, deletion=False, region_trim = 1000, min_percent_selected_genomes=0, tolerance_percentage = 10):
    dic_regions, analyse = find_shared_regions(genomes_list, genome_ref, chromosomes, node_min_size, nodes_max_gap, deletion = deletion, min_percent_selected_genomes=min_percent_selected_genomes, tolerance_percentage = tolerance_percentage)
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
    if get_driver() is None :
        return None
    with get_driver() as driver:
        if type_search == "before" :

            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome)+"""_position` <= $genome_position AND $genome_ref in n.genomes
            return max(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        else:

            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome)+"""_position` >= $genome_position AND $genome_ref in n.genomes
            return min(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        #print("Query : " + query)
        with driver.session() as session:
            result = session.run(query, genome_ref=genome_ref, genome_position=genome_position)
            for record in result :
                ref_position = record["ref_position"]
    
    return ref_position

#Function get_shared_regions : this function is used to get shared regions (positions) between a list of genomes (GWAS)
#On peut préciser sur quel chromosome chercher
#Usage exemple (cattle white spot) : dic_regions, analyse = find_shared_regions(["HER","SIM"],chromosome="6")
#genomes_list : list of the genomes for which the function will look for shared regions
#genome_ref will be used to get annotations on this genome
#chromosomes : liste of chromosomes. If defined the function will only look for shared region on these chromosomes
#node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
#nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
def find_shared_regions(genomes_list, genome_ref=None, chromosomes=None, node_min_size = 10, node_max_size = 0, nodes_max_gap = 10000, deletion=False, min_percent_selected_genomes=100, tolerance_percentage = 0, min_deletion_percentage=100):
    if get_driver is None :
        return {},{}
    dic_regions = {}
    time_0 = time.time()
    if min_percent_selected_genomes > 100:
        min_percent_selected_genomes = 100
    if tolerance_percentage > 100:
        tolerance_percentage = 100
    if min_deletion_percentage > 100:
        min_deletion_percentage = 100
    print("node_min_size : " + str(node_min_size) + " node_max_size : " + str(node_max_size) + " deletion : " + str(deletion) + " min_percent_selected_genomes : "+ str(min_percent_selected_genomes) + " tolerance_percentage : " + str(tolerance_percentage) + " min deletion percentage : " + str(min_deletion_percentage))
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
                
                #max_flow = nb_associated_genomes / nb_genomes + 0.00000001
                #min_flow = (max_flow-0.00000002) * min_percent_selected_genomes / 100  
                
                min_associated_genomes = max(int(min_percent_selected_genomes*nb_associated_genomes/100), 1)
                min_flow = min_associated_genomes/nb_genomes - 0.00000001
                max_flow = nb_associated_genomes*(1+tolerance_percentage/100)/nb_genomes + 0.00000001
                

                print(f"genomes number : {nb_genomes} - min flow : {min_flow} - max flow : {max_flow} - min associated genomes : {min_associated_genomes}")
                if deletion:
                        min_unselected_genomes = max(1,int((nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100))
                        global_min_flow_deletion = min(min_associated_genomes + min_unselected_genomes,nb_genomes)/nb_genomes - 0.00000001 
                        min_flow_deletion = min(1,((nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100)/nb_genomes - 0.00000001)
                        max_flow_deletion = min(1, (nb_genomes - nb_associated_genomes)/nb_genomes) + 0.00000001
                        print(f"Look for deletions with parameters : global min flow : {global_min_flow_deletion} - min flow : {min_flow_deletion} - max flow :  {max_flow_deletion}")
                
                if chromosomes != None :
                    chromosome_list = chromosomes
                else :
                    chromosome_list = get_chromosomes()

                if genome_ref is None or genome_ref == "":
                    genome_position_ref = genomes_list[0]
                    genome_ref = genomes_list[0]
                else :
                    genome_position_ref = genome_ref
                
                print(f"ref genome : {genome_ref}")
                for c in chromosome_list :
                    print("chromosome : " + str(c))
                    dic_regions[c] = {}
                    #Looking for shared nodes
                    if node_max_size == 0:
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = '{c}'
                              AND n.flow >= {min_flow}
                              AND n.flow <= {max_flow}
                              AND n.size >= {node_min_size}
                            WITH n, [g IN n.genomes WHERE g IN $genomes_list] AS matched_genomes
                            WHERE size(matched_genomes) >= {min_associated_genomes}
                              AND size(n.genomes) - size(matched_genomes) <= size(n.genomes) * {tolerance_percentage}/100
                            RETURN n AS nodes order by n.`{genome_position_ref}_position` ASC
                        """
                    else:
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = '{c}'
                              AND n.flow >= {min_flow}
                              AND n.flow <= {max_flow}
                              AND n.size >= {node_min_size}
                              AND n.size <= {node_max_size}
                            WITH n, [g IN n.genomes WHERE g IN $genomes_list] AS matched_genomes
                            WHERE size(matched_genomes) >= {min_associated_genomes}
                              AND size(n.genomes) - size(matched_genomes) <= size(n.genomes) * {tolerance_percentage}/100
                            RETURN n AS nodes order by n.`{genome_position_ref}_position` ASC
                        """
                    #print(query)
                    dic_number_to_position = {}  
                    dic_number_to_size = {}
                    for g in genomes_list :
                        dic_regions[c][g] = {"nodes_position_list":[], "size":[], "regions" : [], "shared_size":[], "deleted_size":[]}
                        dic_number_to_position[g] = {}
                        dic_number_to_size[g] = {}
                        #query += ' AND "' + str(g) + '" IN n.genomes'

                    #print(query)
                    result1 = list(session.run(query, genomes_list=genomes_list))
                    print("Nodes selected for chromosomes " + str(c) + " : " + str(len(result1)) + "\nTime : " + str(time.time()-time_0))
                    if deletion :

                        # query = f"""
                        #     MATCH (n:Node)-[r:gfa_link]->(m:Node)
                        #     where n.chromosome = "{c}"
                        #     AND n.flow >= {min_flow_deletion}
                        #     AND n.size >= {node_min_size}
                        #     with n, [g IN n.genomes WHERE g IN $genomes_list] AS matched_genomes, count(r) as rel_count
                        #     WHERE size(matched_genomes) >= {min_associated_genomes}
                        #     AND size(matched_genomes) < size(n.genomes)
                        #     AND rel_count = 2
                        #     WITH n, [g IN n.genomes WHERE NOT g IN $genomes_list] AS other_genomes
                        #     MATCH (n)-[r1:gfa_link]->(m:Node)
                        #     WHERE ALL(g IN other_genomes WHERE g IN m.genomes)
                        #     AND NONE(g IN $genomes_list WHERE g IN m.genomes)
                        #     AND  m.size  >= {node_min_size}
                        #     Match (n)-[r2:gfa_link]->(m2:Node)
                        #     WHERE ALL(g IN n.genomes WHERE g IN m2.genomes)
                        #     AND size(m2.genomes) = size(n.genomes)
                        #     RETURN n AS nodes, m.size as deleted_node_size order by n.`{genome_position_ref}_position` ASC
                        # """
                        
                        #Query adapted for finding only first deleted node

                        
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = "{c}"
                              AND n.flow >= {min_flow_deletion} AND n.flow <= {max_flow_deletion}
                              AND n.size >= {node_min_size}
                              """
                        if node_max_size > 0:
                             query += f" AND n.size <= {node_max_size}"
                        
                        query += f"""
                                 
                              AND NONE(g IN $genomes_list WHERE g IN n.genomes)
                            
                            OPTIONAL CALL {{
                              WITH n
                              MATCH (m:Node)-[]->(n)
                              WHERE m.chromosome = "{c}"
                                AND m.flow >= {global_min_flow_deletion}
                            
                              WITH m, n,
                                   [g IN $genomes_list WHERE g IN m.genomes AND NOT g IN n.genomes] AS added_genomes
                              WHERE size(added_genomes) >= {min_associated_genomes}
                                AND ALL(g IN n.genomes WHERE g IN m.genomes)
                            
                              WITH m, n
                              WHERE COUNT {{ (m)-[]->(:Node) }} = 2
                            
                              MATCH (m)-[]->(n2:Node)
                              WHERE n2 <> n
                                AND n2.chromosome = "{c}"
                                AND ALL(g IN n2.genomes WHERE g IN m.genomes)
                                AND size(n2.genomes) = size(m.genomes)
                            
                            
                              WITH m, collect(n2) AS nodes_tmp, count(n2) AS nodes_tmp_count
                              WHERE nodes_tmp_count = 1
                              RETURN m, nodes_tmp[0] AS n2
                              
                            }}
                            WITH m, n, n2
                            WHERE m IS NOT NULL and n2 IS NOT NULL 
                            RETURN m as nodes, n.size AS deleted_node_size, n.genomes as genomes_deleted_nodes, n2 as end_deletion_nodes  
                            ORDER BY m.`{genome_position_ref}_position` ASC
                        
                        """
                        
                        result2 = list(session.run(query, genomes_list=genomes_list))
                        result = result1 + result2
                        print("Total Nodes selected for chromosome " + str(c) + " : " + str(len(result)) + "\nTime : " + str(time.time()-time_0))
                    else:
                        result = result1
                    nb_regions_total += len(result)
                    for r in result:
                        for g in genomes_list:
                            if "deleted_nodes" not in dic_regions[c][g]:
                                dic_regions[c][g]["deleted_nodes"]= []
                            if (r["nodes"][g+"_position"] != None):
                                dic_regions[c][g]["nodes_position_list"].append(r["nodes"][g+"_position"]) 

                                if "deleted_node_size" in dict(r): 
                                    #If deleted nodes are found, we try to reconstruct the size of the deleted regions. 
                                    #To do this, we check that there is no overlap.
                                    deleted_nodes_dict = {}
                                    for hap in genomes : 
                                        if hap in r["genomes_deleted_nodes"]:
                                            start_deletion = r["nodes"][hap+"_position"]+r["nodes"]["size"]
                                            end_deletion = r["end_deletion_nodes"][hap+"_position"]
                                            if start_deletion > end_deletion:
                                                end_deletion_tmp = end_deletion
                                                end_deletion = start_deletion
                                                start_deletion = end_deletion_tmp
                                            deleted_nodes_dict[hap]={"start_deletion":start_deletion, "end_deletion":end_deletion}
                                        else:
                                            deleted_nodes_dict[hap] = {"start_deletion":-1, "end_deletion":-1}
                                    dic_regions[c][g]["deleted_nodes"].append(deleted_nodes_dict)
                                    #dic_regions[c][g]["deleted_size"].append(statistics.median(gap))
                                    dic_regions[c][g]["size"].append(0)
                                                                    
                                else:
                                    dic_regions[c][g]["size"].append(r["nodes"]["size"]) 
                                    #dic_regions[c][g]["deleted_size"].append(0)
                                    dic_regions[c][g]["deleted_nodes"].append({})

                    
                    #Group regions if they are separated vy less than nodes_max_gap
                    for g in genomes_list :
                        if len(dic_regions[c][g]["nodes_position_list"]) > 0:
                            combined = list(zip(dic_regions[c][g]["nodes_position_list"],dic_regions[c][g]["size"],dic_regions[c][g]["deleted_nodes"]))
                            combined_sorted = sorted(combined, key=lambda x: x[0])
                            nodes_position_sorted, size_sorted, deleted_nodes_sorted = zip(*combined_sorted)
                            #dic_regions[c][g]["nodes_position_list"].sort()
                            shared_size = 0
                            shared_deleted_size = 0
                            current_deletion = None
                            for i in range(len(nodes_position_sorted)):
                                if i == 0 :
                                    region_start = nodes_position_sorted[0]
                                    region_stop = region_start + size_sorted[0]
                                    shared_size = size_sorted[0]
                                    if len(deleted_nodes_sorted[0]) > 0:
                                        current_deletion = {}
                                        for dg in deleted_nodes_sorted[0] :
                                            current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[0][dg]["start_deletion"], "end_deletion":deleted_nodes_sorted[0][dg]["end_deletion"]}
                                    else:
                                        current_deletion = None

                                else :
                                    if nodes_position_sorted[i] < nodes_position_sorted[i-1] + nodes_max_gap :
                                        region_stop = nodes_position_sorted[i]
                                        shared_size += size_sorted[i]

                                        if len(deleted_nodes_sorted[i]) > 0:
                                            #If deleted nodes are found, we try to reconstruct the size of the deleted regions. 
                                            #To do this, we check that there is no overlap.
                                            same_deletion = False
                                            if current_deletion is None:
                                                same_deletion = True
                                            else:
                                                for dg in deleted_nodes_sorted[i] :
                                                    if "start_deletion" in deleted_nodes_sorted[i][dg] and deleted_nodes_sorted[i][dg]["start_deletion"]>=0 \
                                                        and current_deletion[dg]["start_deletion"] >= 0 \
                                                        and ((deleted_nodes_sorted[i][dg]["start_deletion"] >= current_deletion[dg]["start_deletion"] \
                                                        and deleted_nodes_sorted[i][dg]["start_deletion"] <= current_deletion[dg]["end_deletion"]) \
                                                        or (deleted_nodes_sorted[i][dg]["end_deletion"] >= current_deletion[dg]["start_deletion"] \
                                                        and deleted_nodes_sorted[i][dg]["end_deletion"] <= current_deletion[dg]["end_deletion"])) :
                                                            same_deletion = True
                                            if same_deletion :
                                                for dg in deleted_nodes_sorted[i] :
                                                    if current_deletion is not None :
                                                        if dg in current_deletion :
                                                            current_deletion[dg] = {"start_deletion": min(deleted_nodes_sorted[i][dg]["start_deletion"], current_deletion[dg]["start_deletion"]),
                                                                               "end_deletion": max(deleted_nodes_sorted[i][dg]["end_deletion"], current_deletion[dg]["end_deletion"])}
                                                        else:
                                                            current_deletion[dg] = {"start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                                               "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                                    else:
                                                        current_deletion = {}
                                                        current_deletion[dg] = {"start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                                           "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                            else:
                                                gap = []
                                                for dg in current_deletion :
                                                    if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                                        gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                                shared_deleted_size += statistics.median(gap)
                                                for dg in deleted_nodes_sorted[i] :
                                                    current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[i][dg]["start_deletion"], "end_deletion":deleted_nodes_sorted[i][dg]["end_deletion"]}
                                        
                                    else :
                                        if current_deletion is not None:
                                            gap = []
                                            for dg in current_deletion :
                                                if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                                    gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                            shared_deleted_size += statistics.median(gap)
                                        if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                                            region_stop = region_start + shared_deleted_size
                                        
                                        if region_start == region_stop:
                                            if shared_deleted_size > 0:
                                                region_stop = region_start + shared_deleted_size
                                            else: 
                                                region_stop += 100
                                                region_start -= 100
                                        #Minimal region to allow visualization
                                        if region_stop - region_start < 200 :
                                            min_size_region = 200
                                            gap = int((min_size_region-(region_stop - region_start))/2)
                                            region_start = max(0, region_start-gap)
                                            region_stop = region_stop+gap
                                        dic_regions[c][g]["regions"].append({"start" : region_start, "stop" : region_stop, "shared_size" : shared_size, "shared_deleted_size":shared_deleted_size, "region_size" : region_stop-region_start})
                                        shared_size = size_sorted[i]
                                        shared_deleted_size = 0
                                        if len(deleted_nodes_sorted[i]) > 0:
                                            for dg in deleted_nodes_sorted[i] :
                                                if current_deletion is None :
                                                    current_deletion = {}
                                                current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[i][dg]["start_deletion"], "end_deletion":deleted_nodes_sorted[i][dg]["end_deletion"]}
                                        else:
                                            current_deletion = None
                                        region_start = nodes_position_sorted[i]
                                        region_stop = region_start + size_sorted[i]
                            if current_deletion is not None:
                                gap = []
                                for dg in current_deletion :
                                    if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                        gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                shared_deleted_size += statistics.median(gap)
                            
                            if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                                region_stop = region_start + shared_deleted_size
                            if region_start == region_stop:
                                if shared_deleted_size > 0:
                                    region_stop = region_start + shared_deleted_size
                                else: 
                                    region_stop += 100
                                    region_start -= 100   
                            #Minimal region to allow visualization
                            if region_stop - region_start < 200 :
                                min_size_region = 200
                                gap = int((min_size_region-(region_stop - region_start))/2)
                                region_start = max(0, region_start-gap)
                                region_stop = region_stop+gap
                            dic_regions[c][g]["regions"].append({"start" : region_start, "stop" : region_stop, "shared_size" : shared_size, "shared_deleted_size":shared_deleted_size, "region_size" : region_stop-region_start})
            print("Total number of regions : " +str(nb_regions_total))
            dic_regions_2 = {}
            for c, genomes in dic_regions.items():
                for genome, valeur in genomes.items():
                    if genome not in dic_regions_2:
                        dic_regions_2[genome] = {}
                    dic_regions_2[genome][c] = valeur

            analyse = {}
            print(f"genomes : {dic_regions_2.keys()}")
            for g in tqdm(dic_regions_2):
                print(f"genome : {g}")
                analyse[g] = []
                for c in dic_regions_2[g]:
                    #print(f"genome : {g} - chromosome : {c} - regions number : {len(dic_regions_2[g][c]['regions'])}")
                    for r in dic_regions_2[g][c]['regions']:
                        r["chromosome"] = c
                        r["genome"] = g
                        if (genome_ref is not None and g == genome_ref) or (genome_ref is None and g == genomes_list[0]) :
                            #print("Search annotations")
                            r["annotations"] = get_annotations_in_position_range(genome_ref=g,chromosome=c, start_position=r["start"],end_position=r["stop"])
                            #print("Search annotation before")
                            annot_before_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c, position=r["start"], before=True)
                            annot_tmp = {}
                            if annot_before_tmp is not None and "gene_name" in annot_before_tmp :
                                annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                                annot_tmp["distance"] = annot_before_tmp["end"]-r["start"]
                            r["annotation_before"] = annot_tmp
                            
                            #print("Search annotation after")
                            annot_after_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c, position=r["stop"], before=False)
                            annot_tmp = {}
                            if annot_after_tmp is not None and "gene_name" in annot_after_tmp :
                                annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                                annot_tmp["distance"] = annot_after_tmp["start"]-r["stop"]
                            r["annotation_after"] = annot_tmp
                        analyse[g].append(r)
            if genome_ref not in analyse:
                analyse[genome_ref] = []
                for a in analyse[genomes_list[0]]:
                    r = {}
                    r["chromosome"] = a["chromosome"]
                    r["shared_size"] = a["shared_size"]
                    r["shared_deleted_size"] = a["shared_deleted_size"]
                    r["genome"] = genome_ref
                    n_start = get_anchor(genomes_list[0], a["chromosome"], a["start"], before = True)
                    n_stop = get_anchor(genomes_list[0], a["chromosome"], a["stop"], before = False)
                    position_field = genome_ref+"_position"
                    r["start"] = 0
                    if position_field in n_start :
                        r["start"] = n_start[position_field]
                    if position_field in n_stop:
                        r["stop"] = n_stop[position_field]
                    else:
                        r["stop"] = r["start"]
                    r["region_size"] = r["stop"] - r["start"]
                    r["annotations"] = get_annotations_in_position_range(genome_ref=genome_ref,chromosome=a["chromosome"], start_position=r["start"],end_position=r["stop"])
                    annot_before_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref, chromosome=a["chromosome"], position=r["start"], before=True)
                    annot_tmp = {}
                    if annot_before_tmp is not None and "gene_name" in annot_before_tmp :
                        annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                        annot_tmp["distance"] = annot_before_tmp["end"]-r["start"]
                    r["annotation_before"] = annot_tmp
                    
                    annot_after_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref, chromosome=a["chromosome"], position=r["stop"], before=False)
                    annot_tmp = {}
                    if annot_after_tmp is not None and "gene_name" in annot_after_tmp :
                        annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                        annot_tmp["distance"] = annot_after_tmp["start"]-r["stop"]
                    r["annotation_after"] = annot_tmp
                    analyse[genome_ref].append(r)
                        
            for g in analyse:
                analyse[g] = sorted(analyse[g], key=lambda d: d['shared_size'], reverse=True)
            
            print("Total time : "+ str(time.time()-temps_depart))
    return dic_regions_2, analyse



def calculer_variabilite(chromosome_list=None, ref_genome=None, window_size=1000, output_html="pangenome_variability.html"):
    with get_driver() as driver:
        if ref_genome == None :
            ref_noeud = "node_mean"
            ref_position = "position_mean"
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
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
            WITH FLOOR(n.{ref_position} / {window_size}) AS Window, n.{ref_position} AS position, n.flow AS flow
            WITH Window, 
                 collect(position) AS positions, 
                 avg(flow) AS flow_mean, 
                 count(*) AS nodes_number
            RETURN 
                Window,
                reduce(min_pos = head(positions), p IN positions | CASE WHEN p < min_pos THEN p ELSE min_pos END) AS start_position,
                reduce(max_pos = head(positions), p IN positions | CASE WHEN p > max_pos THEN p ELSE max_pos END) AS end_position,
                flow_mean,
                nodes_number
            ORDER BY Window
            """
        
            with driver.session() as session:
                result = session.run(query)
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
                "<br>Flow mean: " + df['flow_mean'].round(2).astype(str)
            )
        
            # Create figure
            fig = px.scatter(
                df,
                x='Window',
                y='flow_mean',
                hover_name='hover_text',
                labels={
                    'Window': 'Window',
                    'flow_mean': 'Flow mean'
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
            f.write("<html><head><title>Flow variability</title></head><body>\n")
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
    if get_driver() is None :
        return None
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
            if "strandM" in nodes_data[n] and g in nodes_data[n]["strandM"]:
                genome_redondant_strand_matrix[index_genomes[g], i] += 1
                genome_strand_matrix[index_genomes[g], i] = nodes_data[n]["size"]
            if "strandP" in nodes_data[n] and g in nodes_data[n]["strandP"]:
                genome_redondant_strand_matrix[index_genomes[g], i+nodes_nb] += 1
                genome_strand_matrix[index_genomes[g], i+nodes_nb] = nodes_data[n]["size"]
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
    #print(newick_tree)
    return newick_tree
            

"""
Fonction to convert PAV matrix to phylip format

Parameters
    @pav_matrix_dic : PAV matrix (dictionnary)
    @matrice_distance_phylip_filename : output filename (PHYLIP)
returns
"""


def pav_to_phylip(pav_matrix_dic, distance_matrix_phylip_filename):
    # Creating phylip file
    file = open(distance_matrix_phylip_filename, "w")

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




#This function compute a raxml tree from a random selection of nodes
#RaxML output file will be created in output_dir directory
#to take account of direct / reverse traversing set the strand to True, else to False
#To limit the tree to 1 chromosome, set chromosome value (else it will be computed on the whole graph)
#project_name is used to name the RaxML files
#max_nodes is used to limit the number of sampled nodes
#min_sample_size is the minimal sampled size for pangenome of size >= min_sample_size
#min_nodes_number : tree won't be computed on pangenomes with less than min_nodes_number nodes 
def compute_global_raxml_phylo_tree_from_nodes(output_dir = "", strand=True, chromosome = None, project_name="panorama_phylo_tree", max_nodes = 1000000, min_sample_size = 10000, min_nodes_number = 1000):
    total_nodes_number = 0
    

    dir_raxml = "./export/phylo/raxml"
    distance_matrix_filename = "distance_matrix.phy"
    tree_newick_filename = os.path.join(dir_raxml,"RAxML_bestTree."+project_name)
    distance_matrix_phylip_filename = os.path.join(dir_raxml, distance_matrix_filename)
    if not os.path.exists(dir_raxml):
        os.makedirs(dir_raxml)
    if get_driver() is None :
        return None
    
    with get_driver() as driver:
        if driver is None:
            return None
        if chromosome is None :
            query = f"""
            MATCH (n:Node)
            return count(*) as total_nodes_number
            """
        else:
            print(f"Chromosome : {chromosome}")
            query = f"""
            MATCH (n:Node)
            where n.chromosome = '{chromosome}'
            return count(*) as total_nodes_number
            """
 
        with driver.session() as session:
            result = session.run(query)
            for record in result :
                total_nodes_number = record["total_nodes_number"]
                
        #The number of sampled nodes depends on the total number of nodes. It is set to be at least min_sample_size and at most max_nodes.
        #If the pangenome contains less than min_nodes_number then no tree can be computed.
        #The number of sampled nodes is determined by a polynomial interpolation in log-space of the total number of nodes, fitted through a set of control points.
        if total_nodes_number < min_nodes_number:
            return None
        if total_nodes_number < min_sample_size:
            sample_nodes_number = total_nodes_number
        else:
            #Set of control points :
            x_vals = np.array([1e4, 1e5, 1e7, 1e8, 2e9])
            y_vals = np.array([1e4, 5e4, 3e5, 5e5, 1e6])
            logx = np.log10(x_vals)
            #Compute the polynomial interpolation in log space:
            coeffs = np.polyfit(logx, y_vals, deg=4)
            def f(x):
                return np.polyval(coeffs, np.log10(x))
            sample_nodes_number = max(min_sample_size, min(int(f(total_nodes_number)),max_nodes))
        #sample_nodes_number = max(min_sample_size,min(int(node_selection_percentage*total_nodes_number/100), max_nodes))
        
        
        
        if chromosome is None :
            query = f"""
            MATCH (n:Node)
            with n, rand() AS r
            order by r
            limit {sample_nodes_number}
            return distinct(n) as nodes
            """
        else:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = '{chromosome}'
            with n, rand() AS r
            order by r
            limit {sample_nodes_number}
            return distinct(n) as nodes
            """
        nodes_list = []
        with driver.session() as session:
            result = session.run(query)
            for record in result :
                nodes_list.append(dict(record["nodes"]))
        
        
        print(f"Number of sampled nodes : {sample_nodes_number} - Total node {total_nodes_number}")
        sample_size = len(nodes_list)
        if sample_size >= min_sample_size :
            #Prepare PAV matrix
            genomes = get_genomes()
            pav_matrix = {}
            for g in genomes : 
                if g not in pav_matrix:
                    if strand :
                        pav_matrix[g] = np.zeros(2 * sample_size, dtype=int)
                    else:
                        pav_matrix[g] = np.zeros(sample_size, dtype=int)
            
            for i in range(0, len(nodes_list)):
                for g in nodes_list[i]["genomes"]:
                    if strand:
                        if "strandP" in nodes_list[i] and g in nodes_list[i]["strandP"]:
                            pav_matrix[g][i] = int(1)
                        else:
                            pav_matrix[g][i+sample_size] = int(1)
                    else:
                        pav_matrix[g][i] = int(1)         
            pav_to_phylip(pav_matrix, distance_matrix_phylip_filename)
            
            pattern = os.path.join(dir_raxml, f"RAxML_*")
            
            for f in glob.glob(pattern):
                os.remove(f)
            
            raxml_command = [
                "raxmlHPC", 
                "-s", distance_matrix_filename, 
                "-m", "BINGAMMA",  
                "-p", "12345",  
                # '-#', '100',  # Iterations number for bootstrapping
                "-n", project_name
            ]
        
            # launching RaxML command
            result = subprocess.run(raxml_command, check=True, cwd=dir_raxml)
            try:
                with open(tree_newick_filename, 'r') as f:
                    return f.read()
            except FileNotFoundError:
                return None
        else:
            return None


            

#Computes the distance matrix on the whole GFA (could take a long time for big GFA)
def compute_distance_matrix(distance_matrix_filename = "distances.csv", chromosome=None, ponderation=True, strand=False):
    if get_driver() is None :
        return None
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
            dic_size_genome = {}   
            if ponderation :
                #Get the genome size
                query = """
                    MATCH (n:Node)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, sum(n.size) AS total_size
                    """
                result_size_genomes = list(session.run(query))
                if strand == False:
                    #computes size for each genome

                    #Get intersection size
                    query = """
                        MATCH (n:Node)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, sum(n.size) AS size_intersection
                        ORDER BY size_intersection DESC
                        """

                else:
                    query = """
                        MATCH (n:Node)
                        WITH n.size AS size, n.strandM AS strandM, n.strandP AS strandP
                        WITH size,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, size
                        WHERE g1 < g2
                        RETURN g1, g2, sum(size) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                result_intersection = list(session.run(query))
                    
            else:
                #Get the node numbers
                query = """
                    MATCH (n:Node)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, count(*) AS total_size
                """
                result_size_genomes = list(session.run(query)) 
                if strand == False:

                    #Get intersection size
                    query = """
                        MATCH (n:Node)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, count(*) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                else:
                    query = """
                        MATCH (n:Node)
                        WITH n.size AS size, n.strandM AS strandM, n.strandP AS strandP
                        WITH size,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, size
                        WHERE g1 < g2
                        RETURN g1, g2, count(*) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                result_intersection = list(session.run(query))
            
        for r in result_size_genomes:
            dic_size_genome[r["genome"]] = r["total_size"]
        print(dic_size_genome)
        for r in result_intersection:
            g1 = r["g1"]
            g2 = r["g2"]
            inter = r["size_intersection"]
            distance_matrix.loc[g1,g2] = 1-inter/(dic_size_genome[g1]+dic_size_genome[g2]-inter)
            distance_matrix.loc[g2,g1] = distance_matrix.loc[g1,g2]
                
                
        distance_matrix.to_csv(distance_matrix_filename)   
        print("Total time : "+ str(time.time()-temps_depart))
    return distance_matrix



          
 