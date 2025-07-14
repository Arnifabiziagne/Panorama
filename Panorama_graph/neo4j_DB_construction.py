import re
from tqdm import tqdm
from math import *
from neo4j import GraphDatabase
import time
import hashlib
import logging
import scipy.sparse as sp
import numpy as np
import os
from config import get_driver

#Version of BDD
version="1.0.0"

#batch_size_BDD size of batch transaction in DB
batch_size_BDD = 10000


logging.getLogger("neo4j").setLevel(logging.ERROR)




def create_nodes_batch(session, nodes_dic, node_name="Node", create = False):

    nb_transactions = max(1,ceil(len(nodes_dic)/batch_size_BDD))
    current_transaction = 0
    nodes_list = list(nodes_dic.items())
    with tqdm(total=nb_transactions) as bar :
        while len(nodes_dic)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD                
                batch = nodes_list[ind_depart:ind_depart+min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD)]
                if create == False :
                    query = (
                    "UNWIND $batch AS node "
                    "MERGE (n:"+str(node_name) +" {name: node.name}) "
                    "SET n += node.attributes "
                    )
                else :
                    query = (
                    "UNWIND $batch AS node "
                    "CREATE (n:"+str(node_name) +" {name: node.name}) "
                    "SET n += node.attributes "
                    )
                
                batch_data = [{"name": name, "attributes": attributes} for name, attributes in batch]
                tx.run(query, batch=batch_data)
                tx.commit()
    
                bar.update(1)
                current_transaction += 1
    return


def creer_stats(set_genomes, set_chromosomes):
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            with session.begin_transaction() as tx:
                query ="""
                MERGE (s:Stats)
                WITH s, coalesce(s.genomes, []) + $genomes AS all_genomes, $chromosomes as liste_chromosome
                UNWIND all_genomes AS g
                WITH s, collect(DISTINCT g) AS new_genomes, liste_chromosome
                SET s.genomes = new_genomes, s.version=$version
                
                // Mise à jour de s.chromosomes
                WITH s, coalesce(s.chromosomes, []) + liste_chromosome AS all_chromosomes
                UNWIND all_chromosomes AS c
                WITH s, collect(DISTINCT c) AS new_chromosomes
                SET s.chromosomes = new_chromosomes
                """
                tx.run(query, genomes=list(set_genomes), chromosomes = list(set_chromosomes), version=version)
                
    return


def create_base_indexes():
    indexes_queries= [
        "CREATE INDEX NodeIndexName IF NOT EXISTS FOR (n:Node) ON (n.name)",
        "CREATE INDEX NodeIndexChromosome IF NOT EXISTS FOR (n:Node) ON (n.chromosome)"
        ]
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            with session.begin_transaction() as tx:
                for query in indexes_queries :
                    tx.run(query)

def creer_indexes(genome_ref=None):
    indexes_queries = [
        "CREATE INDEX NodeIndexFlow IF NOT EXISTS FOR (n:Node) ON (n.flow)",
        "CREATE INDEX NodeIndexSize IF NOT EXISTS FOR (n:Node) ON (n.size)",
        "CREATE INDEX NodeIndexRefNode IF NOT EXISTS FOR (n:Node) ON (n.ref_node)",
        "CREATE INDEX AnnotationName IF NOT EXISTS FOR (a:Annotation) ON (a.name)",
        "CREATE INDEX AnnotationIndexChromosome IF NOT EXISTS FOR (a:Annotation) ON (a.chromosome)",
        "CREATE INDEX AnnotationIndexStart IF NOT EXISTS FOR (a:Annotation) ON (a.start)",
        "CREATE INDEX AnnotationIndexEnd IF NOT EXISTS FOR (a:Annotation) ON (a.end)",
        "CREATE INDEX AnnotationIndexGeneId IF NOT EXISTS FOR (a:Annotation) ON (a.gene_id)",
        "CREATE INDEX AnnotationIndexGeneName IF NOT EXISTS FOR (a:Annotation) ON (a.gene_name)",
        "CREATE INDEX SequenceIndexSequence IF NOT EXISTS FOR (s:Sequence) ON (s.sequence)"
        ]
    if genome_ref is not None:
        indexes_queries.append("CREATE INDEX NodeIndex"+str(genome_ref)+"_position IF NOT EXISTS FOR (n:Node) ON (n."+str(genome_ref)+"_position)")
        indexes_queries.append("CREATE INDEX NodeIndex"+str(genome_ref)+"_node IF NOT EXISTS FOR (n:Node) ON (n."+str(genome_ref)+"_node)")

                
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            with session.begin_transaction() as tx:
                for query in indexes_queries :
                    tx.run(query)



def creer_index_chromosome(chromosomes_list = []):
    query_genomes = """
    MATCH (s:Stats) 
    RETURN s.genomes as all_genomes
    """
    all_chromosomes = chromosomes_list
    all_genomes = []
    indexes_queries = []
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            result = session.run(query_genomes)
            for record in result:
                all_genomes = record["all_genomes"]

        for c in all_chromosomes:
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexName IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.name)")
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexFlow  IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.flow)")
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexSize IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.size)")
            for g in all_genomes :
                indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"Index"+str(g)+"_position IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n."+str(g)+"_position)")
                indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"Index"+str(g)+"_node IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n."+str(g)+"_node)")
    
        with get_driver() as driver:
            print("Connection established.")
            with driver.session() as session:
                with session.begin_transaction() as tx:
                    for query in indexes_queries :
                        tx.run(query)
                        
                        
def creer_index_chromosome_genomes():
    query_genomes = """
    MATCH (s:Stats) 
    RETURN s.genomes as all_genomes
    """
    indexes_queries = []
    all_genomes = []
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            result = session.run(query_genomes)
            for record in result:
                all_genomes = record["all_genomes"]

            print("Connection established.")
            nb_genomes = len(all_genomes)
            current_genome = 0
            for g in all_genomes :
                print("creating indexes for genome " + g + " ("+str(current_genome) + "/"+str(nb_genomes) +")")
                current_genome += 1
                indexes_queries = []
                indexes_queries.append("CREATE INDEX NodeIndex"+str(g).replace("-", "_")+"_position IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.`"+str(g)+"_position`)")
                indexes_queries.append("CREATE INDEX NodeIndex"+str(g).replace("-", "_")+"_node IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.`"+str(g)+"_node`)")
                with session.begin_transaction() as tx:
                    for query in indexes_queries :
                        tx.run(query)   
                if current_genome % 5 == 0:
                    time.sleep(60*10)
    print("Indexes created")



def create_labels_chromosomes(liste_chromosomes=[]):
    all_chromosomes = []
    labels_queries = []
    with get_driver() as driver :
        if len(liste_chromosomes) == 0:
            query = """
            MATCH (s:Stats) 
            RETURN s.chromosomes as all_chromosomes
            """
            with driver.session() as session:
                result = session.run(query)
                for record in result:
                    all_chromosomes = record["all_chromosomes"]
    
        else:
            all_chromosomes = liste_chromosomes
        for c in all_chromosomes :
            label = "Node_chr_"+str(c)
            labels_queries.append(
                f"""
                CALL apoc.periodic.iterate(
              "MATCH (n:Node) WHERE n.chromosome = '{c}' RETURN n",
              "CALL apoc.create.addLabels(n, ['{label}']) YIELD node RETURN node",
              {{batchSize:1000, parallel:true}}
            )
            """
            )
        with driver.session() as session:
            with session.begin_transaction() as tx:
                for query in labels_queries:
                    print(query)
                    tx.run(query)



def creer_relations_batch(session, liste_relations):
    
    nb_transactions = ceil(len(liste_relations)/batch_size_BDD)
    current_transaction = 0
    with tqdm(total=nb_transactions) as bar :
        while len(liste_relations)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD

                batch = liste_relations[ind_depart:ind_depart+min(batch_size_BDD, len(liste_relations)-current_transaction*batch_size_BDD)]
                query = f"""
                        UNWIND $batch AS pair
                        MATCH (a:Node {{name: pair.depart}})
                        MATCH (b:Node {{name: pair.arrivee}})
                        MERGE (a)-[r:{"lien_gfa"}]->(b)
                        """
                tx.run(query, batch=batch)        
                bar.update(1)
                current_transaction += 1
    return



'''
Cette fonction va créer dans la BDD neo4j les noeuds et leur séquences (uniquement nom et séquence)
Puis elle va créer les indexes : kmer et relations avec les noeuds
'''
def creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, nodes_dic):
    
    nb_transactions = ceil(len(nodes_dic)/batch_size_BDD)
    current_transaction = 0
    nodes_list = list(nodes_dic.keys())
    print("Début de création des séquences en BDD")
    with tqdm(total=len(nodes_dic)) as bar :
        while len(nodes_dic)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD
                nodes_list_transaction = nodes_list[ind_depart:ind_depart+min(batch_size_BDD, len(nodes_dic)-current_transaction*10000)]
    
                print("Transaction " + str(current_transaction))
                for t in range(min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD)):
                    nc = nodes_list_transaction[t]
                    q = """
                        CREATE (n:Sequence {name:$nom, sequence:$sequence})
                        """ 
                    result = tx.run(
                        q,
                        nom=nc,
                        sequence=nodes_dic[nc]
                    )
                bar.update(min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD))
                current_transaction += 1
    print("Fin de création des séquences en BDD")
    liste_kmer = list(dic_kmer_relation.keys())
    print("Début de création des indexes en BDD")
    current_transaction = 0
    with tqdm(total=len(dic_kmer_relation)) as bar :
        while len(dic_kmer_relation)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD
                liste_kmer_transaction = liste_kmer[ind_depart:ind_depart+min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*10000)]
    
                print("Transaction " + str(current_transaction))
                for t in range(min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*batch_size_BDD)):
                    kmer = liste_kmer_transaction[t]
                    q = """
                        CREATE (k:kmer {kmer:$kmer, size:$size})
                        WITH k 
                        MATCH (n:Sequence)
                        WHERE n.name IN $target_nodes
                        CREATE (k)-[:index]->(n)
                        """
                    result = tx.run(
                        q,
                        target_nodes=dic_kmer_relation[kmer],
                        kmer=kmer,
                        size=kmer_size
                    )
                bar.update(min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*batch_size_BDD))
                current_transaction += 1
    print("Fin de création des indexes en BDD")
    return



'''
Cette fonction permet de créer des noeuds avec la séquence du noeud et son nom
Elle créé également des indexes (kmers et noeuds associés)
'''
def charger_sequences_et_indexes(gfa_file_name, kmer_size=31):
    nodes_dic = {}
    file = open(gfa_file_name, "r", encoding='utf-8')
    kmer_set = set()
    dic_kmer_relation = {}
    with file:
        total_nodes = sum(1 for line in file if line.startswith(('S')))
        file.seek(0,0)
        with tqdm(total=total_nodes) as bar :
            ligne = file.readline()
            while ligne:
                if ligne.startswith(('S')):
                    ligne_dec = ligne.split()
                    if len(ligne_dec) > 0:
                        bar.update(1)
                        seq = ligne_dec[2]
                        node = ligne_dec[1]
                        nodes_dic[node] = seq
                        if len(seq) > kmer_size:
                            liste_kmer = [seq[i:i+kmer_size] for i in range(len(seq)-kmer_size+1)]
                            for kmer in liste_kmer:
                                if kmer in kmer_set :
                                    dic_kmer_relation[kmer].append(node)
                                else:
                                    kmer_set.add(kmer)
                                    dic_kmer_relation[kmer] = [node]
                ligne = file.readline()
    file.close()
    return dic_kmer_relation, nodes_dic

'''
Cette fonction permet de créer des noeuds avec la séquence du noeud et son nom
'''
#TODO : découper en lot le chargement des noeuds
def charger_sequences(gfa_file_name, chromosome_file = None, create=False, batch_size=20000000):
    nodes_dic = {}
    start_time = time.time()
    file = open(gfa_file_name, "r", encoding='utf-8')
    with get_driver() as driver:
        print("Connection established.")
        with file:
            total_nodes = sum(1 for line in file if line.startswith(('S')))
            file.seek(0,0)
            with tqdm(total=total_nodes) as bar :
                ligne = file.readline()
                while ligne:
                    if ligne.startswith(('S')):
                        ligne_dec = ligne.split()
                        if len(ligne_dec) > 0:
                            bar.update(1)
                            seq = ligne_dec[2]
                            if chromosome_file != None :
                                node = chromosome_file + "_"+ ligne_dec[1]
                            else : 
                                node = ligne_dec[1]
                            nodes_dic[node] = {"sequence":seq}
                    if len(nodes_dic) >= batch_size:
                        with driver.session() as session:
                            create_nodes_batch(session, nodes_dic, node_name="Sequence", create = create)
                        nodes_dic = {}
                    ligne = file.readline()
        print("Recuperation des noeud terminé en " + str(time.time()-start_time))
        print("Creation des séquences en BDD : " + str(len(nodes_dic)) + " noeuds à créer")
        if len(nodes_dic) > 0 :
            with driver.session() as session:
                create_nodes_batch(session, nodes_dic, node_name="Sequence", create = create)
    print("Creation des séquences terminée en : " + str(time.time()-start_time) + " s")
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
def get_chromosome_genome(WP_line, haplotype=True, chromosome_file = None):
    ligne_dec = WP_line.split()
    if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
        chromosome = "0"
        if ligne_dec[0] == 'P':
            if (len (ligne_dec[1].split("#")) > 1):
                name_dec = ligne_dec[1].split("#")
            else:
                name_dec = ligne_dec[1].split(".")
            if haplotype  and len(name_dec) > 0 :
                genome = name_dec[0]+"_"+name_dec[1]
            else:
                genome = name_dec[0]
            if len(name_dec) > 0 :
                chromosome = re.sub("^0*", "", name_dec[-1].upper().replace("CHR", ""))
        else:
            chromosome = str(ligne_dec[3])
            if haplotype :
                genome = ligne_dec[1]+"_"+ligne_dec[2]
            else : 
                genome = ligne_dec[1]
    if chromosome_file != None and chromosome_file != "":
        chromosome = chromosome_file
    genome = genome.replace("-", "_")
    return chromosome,genome

'''
This function creates nodes in Neo4j from a GFA file.
It processes the data in batches (maximum defined by batch_size) to avoid memory issues.
For batch processing: we iterate through the nodes and select batch_size nodes.
Then we go through the paths to check if we find the corresponding nodes; if yes, we update the node.
At the end of the paths, we create the batch of nodes in the database.
The nodes contain the following properties:
   - Node size
   - Associated chromosome
   - Node index for each sequence passing through the node
   - Position: the cumulative size of the previous nodes for each sequence on the chromosome
   - Master node: in case of sub-node creation to prevent a node from being used multiple times by the same sequence
                  the sub-nodes are named _2, _3, etc. For example, if the master node is s1, the sub-nodes 
                  will be s1_2, s1_3, etc.
   - List of genomes passing through this node
   - List of genomes passing directly through this node
   - List of genomes passing through this node in reverse
   - Flow: the percentage of genomes passing through this node
The tool works with P lines or W lines; however, since the P line format is not well-defined,
especially for chromosome and sequence naming, W lines are preferred.
For W lines: either the file relates to a single chromosome, or the chromosome is defined in the fourth column.
Chromosome information is required to properly use the tool.
Input:
    - gfa_file_name: GFA file to load
    - chromosome_file: if the file relates to a single chromosome, specify the chromosome name (string)
    - chromosome_prefix: if True, the node name is prefixed by the chromosome (this is the case if chromosome_file is provided); if False and no chromosome_file, then the node is not prefixed
    - batch_size: the splitting is done by nodes, the tool will process batches of size "batch_size"
    - start_node: in case there was an issue during loading, we can resume from a specific node
    - create: if it's the first run, this can be set to True; in this case, the tool will create nodes in the DB without checking for existence
              in other cases, it is better to set it to False
    - haplotype: indicates whether the sample name should be concatenated with the haplotype
'''
def load_gfa_node_neo4j(gfa_file_name, chromosome_file = None, chromosome_prefix = False, batch_size = 5000000, start_chromosome = None, create = False, haplotype = True, create_only_relations = False):
    sep = ["[,;.*]", "(<|>)"]
    nb_lots = 0
    set_genome = set()
    set_chromosome = set()
    total_nodes = 0
    total_path = 0
    nodes_size_dic = {}
    temps_depart = time.time()
    set_all_genomes = set()
    set_all_chromosomes = set()
    liste_relations = []
    chromosomes_list = []
    set_relations = set()
    #création des indexes de base (nom et chromosome)
    create_base_indexes() 
    index_first_chromosme = 0
    nodes_set_next_chromosome = set()
    first_chromosome = None
    file = open(gfa_file_name, "r", encoding='utf-8')
    with file:
        
        #Premier parcours pour récupérer les longueurs des noeuds et les genomes
        ligne = file.readline()
        while ligne :
            if ligne.startswith(('S')):
                ligne_dec = ligne.split()
                if (len(ligne_dec) > 0): 
                    nodes_size_dic[ligne_dec[1]]=int(len(ligne_dec[2]))
                    total_nodes += 1
            if ligne.startswith(('P',"W")):
                print(ligne[0:40])
                total_path += 1
                ligne_dec = ligne.split()
                chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                if chromosome not in set_all_chromosomes :
                    chromosomes_list.append(chromosome)
                if len(set_all_chromosomes) == 0:
                    first_chromosome = chromosome
                    if start_chromosome is not None and start_chromosome != "":
                        first_chromosome = start_chromosome
                set_all_genomes.add(genome) 
                set_all_chromosomes.add(chromosome)
                if chromosome == first_chromosome :
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                    else:
                        ind = 6
                        walk = 1
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                
            ligne = file.readline() 
        if first_chromosome is not None and first_chromosome != "":
            for k in range(len(chromosomes_list)):
                if chromosomes_list[k] == first_chromosome:
                    index_first_chromosme = k
        print("Nombre de génomes : " + str(len(set_all_genomes)) + " - Liste des chromosomes : " + str(set_all_chromosomes))
        #If create_only_relations is set to True then the nodes are not processed
        if not create_only_relations :
            print("Debut du parsing, nombre de noeuds : " + str(total_nodes) + "\nChromosome de départ : " + str(start_chromosome))
            for k in range(index_first_chromosme,len(chromosomes_list)) :
                c = chromosomes_list[k]
                nodes_set_chromosome = set(nodes_set_next_chromosome)
                nodes_set_next_chromosome = set()
                print("chromosome " + str(c) + " - number of nodes : " + str(len(nodes_set_chromosome)))
                nb_lots = ceil(len(nodes_set_chromosome)/batch_size)
                current_lot = 0
                while current_lot < nb_lots :
                    temps_0_lot = time.time()
                    nodes_batch_set = set(list(nodes_set_chromosome)[current_lot*batch_size:min(len(nodes_set_chromosome),(current_lot+1)*batch_size)])
                    current_lot += 1
                    print("chromosome " + c + " lot " + str(current_lot) + "/"+str(nb_lots) + " nombre noeuds : " + str(len(nodes_batch_set)))
                    #Parcours du fichier pour récupérer la liste des noeuds à traiter pour ce lot
                    file.seek(0,0)
                    ligne = file.readline()
    
                   #Parcours des chemins pour ce lot 
    
                    nodes_list = []
                    liste_strand = []
                    position_count = {}
                    nodes_count = {}
                    nodes_dic = {}
                    ref_nodes_dic = {}
                    set_genomes_lot = set()
                    #nodes_size_dic = {}
    
                    nodes_set = set()
                    with tqdm(total=total_path) as bar2 :
                        while ligne:
                            ligne_dec = ligne.split()
                            if len(ligne_dec) > 0:
                                #if ligne.startswith(('S')) and ligne_dec[1] in nodes_batch_set:
                                #        nodes_size_dic[ligne_dec[1]]=int(len(ligne_dec[2]))
                                if ligne[0] == 'P' or ligne[0] == 'W':
                                    chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                                    ligne = None
                                    if current_lot == nb_lots and k < len(chromosomes_list) - 1 and chromosome == chromosomes_list[k+1]:
                                        #Dernier lot du chromosome, on prépare les noeuds à traiter pour le chromosome suivant
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                                        else:
                                            ind = 6
                                            walk = 1
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                                    
                                    if chromosome == c:
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [chaine[-1] for chaine in nodes_list]
                                            nodes_list = [chaine[:-1] for chaine in nodes_list]
                                            # if chromosome_file == None :
                                            #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                            #     nodes_list = [chaine[:-1] for chaine in nodes_list]
                                            # else :
                                            #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                            #     nodes_list = [chromosome_file+"_"+chaine[:-1] for chaine in nodes_list]                              
                                        else:
                                            ind = 6
                                            walk = 1
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 != 0]
                                            nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                            # if chromosome_file == None :
                                            #     nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                            # else:
                                            #     nodes_list = [chromosome_file+"_"+nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
          
                                        if chromosome is not None and chromosome not in set_chromosome:
                                            set_chromosome.add(chromosome)
                                        if genome != "_MINIGRAPH_":
                                            if genome not in set_genome:
                                                set_genome.add(genome)
                                            if genome not in set_genomes_lot :
                                                set_genomes_lot.add(genome)
                                                nodes_count[genome] = {}
                                                position_count[genome] = {}
                                            if chromosome not in nodes_count[genome] :
                                                nodes_count[genome][chromosome] = 0
                                                position_count[genome][chromosome] = {}
                                                position_count[genome][chromosome]["current_position"] = 0
                                                position_count[genome][chromosome]["previous_position"] = 0
                                                position_count[genome][chromosome]["current_contig"] = ""
                                            #Dans le cas du Walk, la position de départ est indiquée => on l'utilise
                                            if ind == 6:
                                                if position_count[genome][chromosome]["current_contig"] != ligne_dec[3] :
                                                    #On change de contig, on va ajouter le départ du contig suivant
                                                    position_count[genome][chromosome]["current_position"] += int(ligne_dec[4])
                                                else :
                                                    #On est dans le même contig, on ajoute les éventuels gaps
                                                    if position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4]) > 0 :
                                                        position_count[genome][chromosome]["current_position"] += position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4])
                                                position_count[genome][chromosome]["current_contig"] = ligne_dec[3]
                                                position_count[genome][chromosome]["previous_position"] = int(ligne_dec[5])
        
                                            node = ""
                                            ref_node = ""
                                            strand = ""
                                            ligne_dec=None
                                            #Linéarisation du graphe
                                            # Parcours de la liste des noeuds : si un noeud existe déjà dans un chemin pour la même séquence
                                            # Alors on va créer un nouveau noeud (par exemple si c'est la 6ème itération pour le noeud S1 on va créer S1_6)
                                            for i in range(len(nodes_list)):
                                                if i > 0:
                                                    previous_node = node
                                                node = nodes_list[i]
                                                ref_node = node
                                                size = nodes_size_dic[ref_node]
                                                if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                                    node = chromosome + "_" + node
                                                #On ne traite le noeud que s'il fait parti du lot
                                                
                                                
                                                if ref_node in nodes_batch_set:
                                                    if chromosome_file is not None and chromosome_file != "":
                                                        ref_node = chromosome + "_" + ref_node
                                                    strand = ""
                                                    if (i < len(liste_strand) and liste_strand[i] in ["-", "<"]):
                                                        strand = "M"
                                                    else:
                                                        strand = "P"
                                                    
                                                    if node not in nodes_set:
                                                        if strand == "P" :
                                                            nodes_dic[node] = {"genomes":[genome], "max":1, "strandP":[genome], "strandM":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                        else :
                                                            nodes_dic[node] = {"genomes":[genome], "max":1, "strandM":[genome], "strandP":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                        nodes_set.add(node)
                                                    else:
                                                        if genome not in nodes_dic[node]["genomes"] and chromosome == nodes_dic[node]["chromosome"]:
                                                            nodes_dic[node]["genomes"].append(genome)
                                                            nodes_dic[node]["strand"+strand].append(genome)
                                                            nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]
                                                            nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]
                                                            if position_count[genome][chromosome]["current_position"] < nodes_dic[node]["position_min"]:
                                                                nodes_dic[node]["position_min"] = position_count[genome][chromosome]["current_position"]
                                                            if position_count[genome][chromosome]["current_position"] > nodes_dic[node]["position_max"]:
                                                                nodes_dic[node]["position_max"] = position_count[genome][chromosome]["current_position"]
                                                        else :
                                                            #Le noeud est redondant, on va chercher si un noeud de même séquence
                                                            #est disponible, sinon on créé un nouveau noeud
                                                            if ref_node not in ref_nodes_dic:
                                                                ref_nodes_dic[ref_node] = {}
                                                            if genome+"-"+chromosome not in ref_nodes_dic[ref_node] :
                                                                ref_nodes_dic[ref_node][genome+"-"+chromosome] = 2
                                                            else:
                                                                ref_nodes_dic[ref_node][genome+"-"+chromosome] += 1
                                                            node = node + "_" + str(ref_nodes_dic[ref_node][genome+"-"+chromosome])
                                                            if node in nodes_set:
                                                                nodes_dic[node]["genomes"].append(genome)
                                                                nodes_dic[node]["strand"+strand].append(genome)
                                                                nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]
                                                                nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]
                                                                if position_count[genome][chromosome]["current_position"] < nodes_dic[node]["position_min"]:
                                                                    nodes_dic[node]["position_min"] = position_count[genome][chromosome]["current_position"]
                                                                if position_count[genome][chromosome]["current_position"] > nodes_dic[node]["position_max"]:
                                                                    nodes_dic[node]["position_max"] = position_count[genome][chromosome]["current_position"]
                                                            else:
                                                                if strand == "P" :
                                                                    nodes_dic[node] = {"genomes":[genome], "max":1, "ref_node" : ref_node, "strandP":[genome], "strandM":[], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                                else :
                                                                    nodes_dic[node] = {"genomes":[genome], "max":1, "ref_node" : ref_node, "strandM":[genome], "strandP":[], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                                nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]   
                                                                nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]   
                                                                nodes_dic[node]["max"]+=1
                                                                nodes_set.add(node)
                                                                
                                                nodes_count[genome][chromosome] += 1
                                                position_count[genome][chromosome]["current_position"] += size
                                    bar2.update(1)
                            ligne = file.readline() 
                    nodes_list = None
                    #Flow computing
                    print("\nSize of elements to create into DB : " + str(len(list(nodes_dic.items()))))
                    if len(nodes_dic) > 0:
                        for node in nodes_dic:
                            nodes_dic[node]["flow"] = len(nodes_dic[node]["genomes"])/len(set_all_genomes)
                            node_mean = 0
                            position_mean = 0
                            nb_genomes = 0
                            for g in nodes_dic[node]["genomes"]:
                                nb_genomes += 1
                                node_mean += nodes_dic[node][g+"_node"]
                                position_mean += nodes_dic[node][g+"_position"]
                            nodes_dic[node]["position_mean"] = int(position_mean/nb_genomes)
                            nodes_dic[node]["node_mean"] = int(node_mean/nb_genomes)
                    
                        print("Time of batch analysis : "+ str(time.time()-temps_0_lot) + "\nTotal time : " + str(time.time()-temps_depart))
                        #Create batch nodes in DB
                        print("Lot " + str(current_lot) + " : Creating nodes in DB")
                        with get_driver() as driver:
                            print("Connection established.")
                            with driver.session() as session:
                                create_nodes_batch(session, nodes_dic, create=create)
                        nodes_dic = None
                        ref_nodes_dic = None
                        
                        
            creer_stats(set_genome, set_chromosome)
            nodes_size_dic = None
            print("Creation of nodes terminated\nTotal time : " + str(time.time()-temps_depart) + "\nGenomes analysed : " + str(set_genome) + "\nNodes number : "+str(total_nodes) )
        
        #Les relations sont traitées après car le découpage en noeud ne permet pas de garantir 
        #qu'on va créer les relations entre les bons noeuds en cas de noeuds doublonnés
        #Pour les relations on va donc traiter par chromosome
        print("\nDébut du traitement des relations")
        
        set_relations = set()
        liste_relations = []
        nodes_list = []
        liste_strand = []
        current_lot = 0
        for k in range(index_first_chromosme,len(chromosomes_list)) :
            c = chromosomes_list[k]
            repeat_nodes = {}
            file.seek(0,0)
            ligne = file.readline()
            current_lot += 1
            print("chromosome " + str(current_lot) + " / " + str(len(chromosomes_list)-index_first_chromosme))
            with tqdm(total=total_path) as bar3 :
                while ligne:
                    ligne_dec = ligne.split()
                    if len(ligne_dec) > 0:
                        if ligne[0] == 'P' or ligne[0] == 'W':
                            chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                            if chromosome == c:
                                if ligne_dec[0] == 'P':
                                    ind = 2
                                    walk = 0
                                    nodes_list = re.split(sep[walk], ligne_dec[ind])
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        liste_strand = [chaine[-1] for chaine in nodes_list]
                                        nodes_list = [chromosome_file+"_"+chaine[:-1] for chaine in nodes_list]   
                                    else :
                                        liste_strand = [chaine[-1] for chaine in nodes_list]
                                        nodes_list = [chaine[:-1] for chaine in nodes_list]                           
                                else:
                                    ind = 6
                                    walk = 1
                                    nodes_list = re.split(sep[walk], ligne_dec[ind])
                                    liste_strand = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 != 0]
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        nodes_list = [chromosome_file+"_"+nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                    else:
                                        nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
        
        
                                if genome != "_MINIGRAPH_":
                                    node = ""
                                    ref_node = ""
                                    strand = ""
                                    # Parcours de la liste des noeuds : si un noeud existe déjà dans un chemin pour la même séquence
                                    # Alors on va créer un nouveau noeud (par exemple si c'est la 6ème itération pour le noeud S1 on va créer S1_6)
                                    for i in range(len(nodes_list)):
                                        if i > 0:
                                            previous_node = node
                                        node = nodes_list[i]
                                        if node in repeat_nodes :
                                            if genome+"-"+chromosome in repeat_nodes[node]:
                                                repeat_nodes[node][genome+"-"+chromosome] += 1
                                                node = node+"_"+ str(repeat_nodes[node][genome+"-"+chromosome])
                                            else:
                                                repeat_nodes[node][genome+"-"+chromosome] = 1
                                        else:
                                            repeat_nodes[node] = {genome+"-"+chromosome:1}
                                        if i > 0 and previous_node+"->"+node not in set_relations :
                                            set_relations.add(previous_node+"->"+node)
                                            #liste_relations.append({"depart":noeud_precedent,"arrivee":noeud})
                            bar3.update(1)
                    ligne = file.readline() 

            for rel in set_relations:
                liste_relations.append({"depart":rel.split("->")[0],"arrivee":rel.split("->")[1]})
            set_relations = set()
        
            print("Lot : " + str(current_lot) + " création des relations en BDD, nombre à créer : " + str(len(liste_relations)))
            with get_driver() as driver:
                print("Connection established.")
                with driver.session() as session:
                    creer_relations_batch(session, liste_relations)
            liste_relations = []
        print("Fin de création des relations\nDurée du traitement : " + str(time.time()-temps_depart) + "\nRelations créées : " + str(len(set_relations)))
        set_relations = set()
        
        #Fin des lots, on créé toutes les relations                                               
    file.close()

    print("Chargement gfa terminée\nDurée totale du traitement : "+ str(time.time()-temps_depart))
    return set_genome


def get_chromosome_annotation(annotation):
    chromosome = None
    tab = re.split(r'[_,.#]', annotation)
    for i in range(len(tab)):
        t = tab[i]
        if t.lower().startswith(("chr")):
            if len(t) > 3:
                chromosome = t[3:].lstrip('0')
            else:
                if i < len(tab)-1:
                    chromosome = tab[i+1].lstrip('0')
    return chromosome
    


#This function will create the Annotation nodes in the neo4j database from a gtf file, but without creating the relationships.
def charger_annotations_neo4j(annotations_file_name, genome_ref, node_name="Annotation", single_chromosome = None):
    temps_depart = time.time()
    file = open(annotations_file_name, "r", encoding='utf-8')
    nodes_dic = {}
    file_format = "gtf"
    with file:
        n = 0
        for line in file :
            n +=1
        file.seek(0,0)
        with tqdm(total=n) as bar :
            #Creating annotations nodes
            ligne = file.readline()
            while ligne :
                if ligne[0] != '#':
                    ligne_dec = ligne.split()
                    chromosome = get_chromosome_annotation(ligne_dec[0])
                    if single_chromosome == None or single_chromosome == chromosome : 
                        name = hashlib.sha256(ligne.encode("utf-8")).hexdigest()
                        nodes_dic[name] = {}
                        
                        feature = ligne_dec[2].lower()
                        nodes_dic[name]["name"] = name
                        nodes_dic[name]["chromosome"] = chromosome
                        nodes_dic[name]["genome_ref"] = genome_ref
                        nodes_dic[name]["source"] = ligne_dec[1]
                        nodes_dic[name]["feature"] = feature
                        
                        start = int(ligne_dec[3])
                        end = int(ligne_dec[4])
                        nodes_dic[name]["start"] = start
                        nodes_dic[name]["end"] = end
                        strand = ligne_dec[6]
                        frame = ligne_dec[7]
                        if len(ligne_dec) == 9:
                            attributes = re.split(r"[;=]", ligne_dec[8])
                            file_format = "gff"
                        else:
                            attributes = ligne_dec[8:]
                        if file_format == "gtf":
                            for i in range(len(attributes)):
                                match attributes[i].lower() :
                                    case "gene_id":
                                        nodes_dic[name]["gene_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_version":
                                        nodes_dic[name]["gene_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_id":
                                        nodes_dic[name]["transcript_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_version":
                                        nodes_dic[name]["transcript_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "exon_number":
                                        nodes_dic[name]["exon_number"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_name":
                                        nodes_dic[name]["gene_name"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_source":
                                        nodes_dic[name]["gene_source"] = attributes[i+1][:-1].replace('"',"")    
                                    case "gene_biotype":
                                        nodes_dic[name]["gene_biotype"] = attributes[i+1][:-1].replace('"',"")     
                                    case "transcript_name":
                                        nodes_dic[name]["transcript_name"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_source":
                                        nodes_dic[name]["transcript_source"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_biotype":
                                        nodes_dic[name]["transcript_biotype"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_id":
                                        nodes_dic[name]["protein_id"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_version":
                                        nodes_dic[name]["protein_version"] = attributes[i+1][:-1].replace('"',"") 
                                    case "tag":
                                        nodes_dic[name]["tag"] = attributes[i+1][:-1].replace('"',"") 
                        else:
                            if feature == "gene" :
                                current_gene_id = ""
                                current_gene_name = ""
                                current_transcript_id = ""
                                current_transcript_name = ""
                            if feature == "exon":
                                exon_id = ""
                                exon_number = ""
                                
                            for i in range(len(attributes)):
                                current_attribute = attributes[i].lower()
                                if feature == "gene" and current_attribute =="id":
                                    current_gene_id = attributes[i+1]
                                elif feature == "gene" and current_attribute =="gene_id":
                                    current_gene_id = attributes[i+1]
                                elif feature == "gene" and current_attribute =="name":
                                    current_gene_name = attributes[i+1]            
                                elif feature == "mrna" and current_attribute == "id":
                                    current_transcript_id = attributes[i+1]
                                elif feature == "mrna" and current_attribute == "transcript_id":
                                    current_transcript_id = attributes[i+1]
                                elif feature == "mrna" and current_attribute == "name":
                                    current_transcript_name = attributes[i+1]
                                elif feature == "exon" and current_attribute == "id":
                                    exon_id = attributes[i+1]
                                    exon_number = exon_id.split(".")[-1]
                                    nodes_dic[name]["exon_id"] = exon_id
                                    nodes_dic[name]["exon_number"] = exon_number 
                                else :
                                    if current_attribute == "id":
                                        nodes_dic[name]["id"] = attributes[i+1]
                                    elif current_attribute == "name":
                                        nodes_dic[name]["name"] = attributes[i+1]
                                        
                            if current_gene_id != "" :
                                nodes_dic[name]["gene_id"] = current_gene_id
                                nodes_dic[name]["gene_name"] = current_gene_name
                                if current_transcript_id != "":
                                    nodes_dic[name]["transcript_id"] = current_transcript_id
                                    nodes_dic[name]["transcript_name"] = current_transcript_name
                            
                                    
                            
                
                bar.update(1)
                ligne = file.readline()

                    
           
            print("\nSize of elements to create in DB : " + str(len(list(nodes_dic.items()))))

            print("Time for batch analyzing : "+ str(time.time()-temps_depart))


            with get_driver() as driver:
                print("Connection established.")
                with driver.session() as session:
                    create_nodes_batch(session, nodes_dic, node_name=node_name)

        print("Nodes created\nTime : " + str(time.time()-temps_depart))
        #Fin des lots, on créé toutes les relations                                               
    file.close()
    

    print("Annotation load terminated.\nTotal time : "+ str(time.time()-temps_depart))


def process_annotation_simple_batch(tx, annotations, genome_ref):
    #for a in annotations:
        #print(a)
    tx.run("""
        UNWIND $annotations AS annot
        MATCH (a:Annotation {name: annot.name})
        MATCH (n:Node)
        WHERE n.chromosome = annot.chromosome 
        AND n.`"""+str(genome_ref)+"""_position` >= annot.start 
        AND n.`"""+str(genome_ref)+"""_position` <= annot.end
        MERGE (n)-[:A_POUR_ANNOTATION]->(a) """, annotations=annotations)

def process_annotation_complexe_batch(tx, last_id, genome_ref, batch_size = 100000):
    
    result = tx.run(
        """
        MATCH (n:Node)
        WHERE id(n) > $last_id
        WITH n ORDER BY id(n) ASC LIMIT $limit
        OPTIONAL MATCH (a:Annotation)
        WHERE a.chromosome = n.chromosome
          AND a.start > n.`"""+str(genome_ref)+"""_position`
          AND a.start < n.`"""+str(genome_ref)+"""_position` + n.size
        WITH n, collect(a) AS annotations
        UNWIND annotations AS a
        MERGE (n)-[:A_POUR_ANNOTATION]->(a)
        RETURN count(*) AS relations_created
        """,
        last_id=last_id, limit=batch_size).single()

    return result["relations_created"]

#This function will create relationships between nodes and annotations in the database. 
#processing is divided into two types of relationships:   
#- Simple relationships, where the position of a node (for a reference genome) lies between the start and end of an annotation.
#- Complex relationships, where the start/end of an annotation lies between the start/end of a node
#Parameters : 
#   - genome_ref : only nodes on genome_ref will be analyzed
def creer_relations_annotations_neo4j(genome_ref, chromosome=None):
    temps_depart = time.time()
    with get_driver() as driver:
        last_id = -1
        batch_size = 10000
        
        
        with driver.session() as session:
            #Traitement des annotations simples : celles pour lesquelles le début est entre le début et la fin d'un noeud
            #c'est le plus gros volume d'annotation
            print("Processing simple annotations")
            i = 0
            last_name = None
            last_id = -1
            max_id = batch_size
            min_id = -1
            current_id = -1
            

            query = """
            MATCH (a:Annotation) return min(ID(a)) as min_id, max(ID(a)) as max_id
            """

            result = session.run(query)
            for record in result:
                min_id = record["min_id"]
                max_id = record["max_id"]
            print("min id : " + str(min_id), "max_id : " + str(max_id))
            current_id = min_id
            while current_id < max_id:
                i+=1
                print("traitement du lot " + str(i) + " Current id : " + str(current_id) + " max id : " + str(max_id) )
                if chromosome is None :
                    annotations = session.run(
                        """
                        MATCH (a:Annotation)
                        WHERE ID(a) >= $min_id AND ID(a) < $max_id
                        RETURN a.name AS name, a.chromosome AS chromosome, a.start AS start, a.end AS end

                        """,
                        min_id=current_id,
                        max_id=current_id+batch_size
                        ).data()                    
                else:
                    annotations = session.run(
                        """
                        MATCH (a:Annotation)
                        WHERE ID(a) >= $min_id AND ID(a) < $max_id and a.chromosome = chromosome
                        RETURN a.name AS name, a.chromosome AS chromosome, a.start AS start, a.end AS end
                        """,
                        min_id=current_id,
                        max_id=current_id+batch_size,
                        chromosome=chromosome
                        ).data()
                current_id += batch_size
                current_id = min(max_id, current_id)
                print("Nombre annotations : " + str(len(annotations)))
                if annotations :
                    with session.begin_transaction() as tx:
                        print("creating annotation lot " + str(i))
                        process_annotation_simple_batch(tx,annotations, genome_ref)
                        #process_annotation_complexe_batch(tx,annotations)
                        tx.commit()

            print("Processing complex annotations")
            #Traitement des annotations complexes : celles pour lesquelles le début et la fin d'un noeud sont avant et après l'annotation
            #le volume est beaucoup plus faible (moins de 1%)
            total_created = 0
            last_id = -1
            batch_size = 100000
            result = session.run("MATCH (n:Node) RETURN max(id(n)) AS max_id")
            max_id = result.single()["max_id"]
            print("max id : " + str(max_id))
            while last_id is not None and last_id < max_id :

                created = session.execute_write(process_annotation_complexe_batch, last_id, genome_ref, batch_size)
                
                result = session.run(
                    """
                    MATCH (n:Node)
                    WHERE id(n) > $last_id
                    WITH n ORDER BY id(n) ASC LIMIT $limit
                    return max(id(n)) as max_id
                    """,last_id=last_id, limit=batch_size)
                last_id = result.single()["max_id"]
                print(f"Traitement jusqu'à id {last_id} sur {max_id} → {created} relations créées.")
                total_created += created


    print("End of relationships creation. Total time : " + str(time.time()-temps_depart))
    


#Main function to construct the whole db from the gfa file
#If the gfa relate to a single chromosome, chromosome_file must contains the reference of this chromosome (1, 2, X, Y, etc.)
#batch_size value is important to limit memory usage, according to the memory available it can be necessary to reduce this value for big pangenomes graphs.
#genome_ref is required if an annotation_file_name is present : this name is used to link the annotations nodes with the main nodes of the graph.
def construct_DB(gfa_file_name, annotation_file_name = None, genome_ref = None, chromosome_file = None, chromosome_prefix = False, batch_size = 5000000, start_chromosome = None, create = False, haplotype = True, create_only_relations = False):
    start_time = time.time()
    charger_sequences(gfa_file_name, chromosome_file, create=create)
    sequence_time = time.time()
    print("Sequences loaded in " + str(sequence_time-start_time) + " s")
    load_gfa_node_neo4j(gfa_file_name, chromosome_file = chromosome_file,  chromosome_prefix = chromosome_prefix, batch_size = batch_size, create = create, start_chromosome = start_chromosome, haplotype=haplotype, create_only_relations=create_only_relations)
    graph_time = time.time()
    print("Graph loaded in " + str(graph_time-sequence_time) + " s")
    creer_indexes(genome_ref)
    #Sleep time due to time construction of indexes : if too much indexes are created in the same time
    #the creation will fail
    time.sleep(600)
    creer_index_chromosome_genomes()
    index_time = time.time()
    print("Indexes created in " + str(index_time-graph_time) + " s")
    if annotation_file_name != None and genome_ref != None :
        charger_annotations_neo4j(annotation_file_name, genome_ref = genome_ref, single_chromosome = chromosome_file)
        annotation_time = time.time()
        print("Annotations loaded in " + str(annotation_time-index_time) + " s")
        creer_relations_annotations_neo4j(genome_ref)
        annotation_relation_time = time.time()
        print("Annotations relations loaded in " + str(annotation_relation_time-annotation_time) + " s")
    print("Process terminated. BDD construct in " + str(time.time()-start_time) + " s")


#This function load multiple gfa files : each filme must relate to a single chromosome
#The files must be named so that last character before .gfa extension contains the reference of the chromosome
#Exmples : chr1.gfa, chr01.gfa, chromosome_1.gfa, exemple_chr_X.gfa, etc.
def construct_db_by_chromosome(gfa_chromosomes_dir, annotation_file_name = None, genome_ref = None, chromosome_file = None, start_node = 0, batch_size = 5000000, create=False):
    start_time = time.time()
    for gfa_file_name in  os.listdir(gfa_chromosomes_dir):
        if gfa_file_name.endswith(".gfa"):
            if "_" in gfa_file_name:
                chromosome = gfa_file_name[:-4].split("_")[-1]
            else:
                chromosome = gfa_file_name[:-4]

            chromosome = chromosome.lower().removeprefix("chr")
            chromosome = chromosome.lstrip("0")
            print("Loading chromosome : " + str(chromosome))
            if chromosome != "" :
                charger_sequences(gfa_file_name, chromosome_file=chromosome, create=create)
                load_gfa_node_neo4j(gfa_file_name, chromosome_file = chromosome_file, batch_size = batch_size, start_node = start_node, create = create)
    db_time = time.time()
    print("Sequences loaded in " + str(db_time-start_time) + " s")
    creer_indexes(genome_ref)
    #Sleep time due to time construction of indexes : if too much indexes are created in the same time
    #the creation will fail
    time.sleep(600)
    creer_index_chromosome_genomes()
    index_time = time.time()
    print("Indexes created in " + str(index_time-graph_time) + " s")
    if annotation_file_name != None and genome_ref != None :
        charger_annotations_neo4j(annotation_file_name, genome_ref = genome_ref, single_chromosome = chromosome_file)
        annotation_time = time.time()
        print("Annotations loaded in " + str(annotation_time-index_time) + " s")
        creer_relations_annotations_neo4j(genome_ref)
        annotation_relation_time = time.time()
        print("Annotations relations loaded in " + str(annotation_relation_time-annotation_time) + " s")
    print("Process terminated. BDD construct in " + str(time.time()-start_time) + " s")

                
    
    
            
def contruire_sequences_et_indexes_bdd(gfa_file_name, kmer_size=31):
    dic_kmer_relation, nodes_dic = charger_sequences_et_indexes(gfa_file_name, kmer_size=31)
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, nodes_dic)
        
def launcher():
    for c in ["18", "19", "20", "21", "22", "X", "Y"]:
        print("chargement du chromosome : " + str(c))
        set_genomes = load_gfa_node_neo4j("/home/fgraziani/work/project/Pangenomique/GFA/humain/hprc-mc/chr"+str(c)+".gfa", chromosome_file = str(c), batch_size=1000000, create=False)
    