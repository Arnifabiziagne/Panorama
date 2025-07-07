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


#taille_lot_BDD fixe la taille des données envoyées en BDD
taille_lot_BDD = 10000


logging.getLogger("neo4j").setLevel(logging.ERROR)




def creer_noeuds_batch(session, dic_noeuds, node_name="Noeud", create = False):

    nb_transactions = max(1,ceil(len(dic_noeuds)/taille_lot_BDD))
    current_transaction = 0
    liste_noeuds = list(dic_noeuds.items())
    with tqdm(total=nb_transactions) as bar :
        while len(dic_noeuds)-current_transaction*taille_lot_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*taille_lot_BDD                
                batch = liste_noeuds[ind_depart:ind_depart+min(taille_lot_BDD, len(dic_noeuds)-current_transaction*taille_lot_BDD)]
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
                SET s.genomes = new_genomes
                
                // Mise à jour de s.chromosomes
                WITH s, coalesce(s.chromosomes, []) + liste_chromosome AS all_chromosomes
                UNWIND all_chromosomes AS c
                WITH s, collect(DISTINCT c) AS new_chromosomes
                SET s.chromosomes = new_chromosomes
                """
                tx.run(query, genomes=list(set_genomes), chromosomes = list(set_chromosomes))
                
    return


def create_base_indexes():
    indexes_queries= [
        "CREATE INDEX NoeudIndexName IF NOT EXISTS FOR (n:Noeud) ON (n.name)",
        "CREATE INDEX NoeudIndexChromosome IF NOT EXISTS FOR (n:Noeud) ON (n.chromosome)"
        ]
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            with session.begin_transaction() as tx:
                for query in indexes_queries :
                    tx.run(query)

def creer_indexes(genome_ref=None):
    indexes_queries = [
        "CREATE INDEX NoeudIndexFlux IF NOT EXISTS FOR (n:Noeud) ON (n.flux)",
        "CREATE INDEX NoeudIndexTaille IF NOT EXISTS FOR (n:Noeud) ON (n.taille)",
        "CREATE INDEX NoeudIndexTaille IF NOT EXISTS FOR (n:Noeud) ON (n.ref_node)",
        "CREATE INDEX AnnotationName IF NOT EXISTS FOR (a:Annotation) ON (a.name)",
        "CREATE INDEX AnnotationIndexChromosome IF NOT EXISTS FOR (a:Annotation) ON (a.chromosome)",
        "CREATE INDEX AnnotationIndexStart IF NOT EXISTS FOR (a:Annotation) ON (a.start)",
        "CREATE INDEX AnnotationIndexEnd IF NOT EXISTS FOR (a:Annotation) ON (a.end)",
        "CREATE INDEX AnnotationIndexGeneId IF NOT EXISTS FOR (a:Annotation) ON (a.gene_id)",
        "CREATE INDEX AnnotationIndexGeneName IF NOT EXISTS FOR (a:Annotation) ON (a.gene_name)"
        ]
    if genome_ref is not None:
        indexes_queries.append("CREATE INDEX NoeudIndex"+str(genome_ref)+"_position IF NOT EXISTS FOR (n:Noeud) ON (n."+str(genome_ref)+"_position)")
        indexes_queries.append("CREATE INDEX NoeudIndex"+str(genome_ref)+"_noeud IF NOT EXISTS FOR (n:Noeud) ON (n."+str(genome_ref)+"_noeud)")

                
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
            indexes_queries.append("CREATE INDEX Noeud_chr_" + str(c) +"IndexName IF NOT EXISTS FOR (n:Noeud_chr_" + str(c) +") ON (n.name)")
            indexes_queries.append("CREATE INDEX Noeud_chr_" + str(c) +"IndexFlux  IF NOT EXISTS FOR (n:Noeud_chr_" + str(c) +") ON (n.flux)")
            indexes_queries.append("CREATE INDEX Noeud_chr_" + str(c) +"IndexTaille IF NOT EXISTS FOR (n:Noeud_chr_" + str(c) +") ON (n.taille)")
            for g in all_genomes :
                indexes_queries.append("CREATE INDEX Noeud_chr_" + str(c) +"Index"+str(g)+"_position IF NOT EXISTS FOR (n:Noeud_chr_" + str(c) +") ON (n."+str(g)+"_position)")
                indexes_queries.append("CREATE INDEX Noeud_chr_" + str(c) +"Index"+str(g)+"_noeud IF NOT EXISTS FOR (n:Noeud_chr_" + str(c) +") ON (n."+str(g)+"_noeud)")
    
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
                indexes_queries.append("CREATE INDEX NoeudIndex"+str(g).replace("-", "_")+"_position IF NOT EXISTS FOR (n:Noeud) ON (n.chromosome, n.`"+str(g)+"_position`)")
                indexes_queries.append("CREATE INDEX NoeudIndex"+str(g).replace("-", "_")+"_noeud IF NOT EXISTS FOR (n:Noeud) ON (n.chromosome, n.`"+str(g)+"_noeud`)")
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
            label = "Noeud_chr_"+str(c)
            labels_queries.append(
                f"""
                CALL apoc.periodic.iterate(
              "MATCH (n:Noeud) WHERE n.chromosome = '{c}' RETURN n",
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
    
    nb_transactions = ceil(len(liste_relations)/taille_lot_BDD)
    current_transaction = 0
    with tqdm(total=nb_transactions) as bar :
        while len(liste_relations)-current_transaction*taille_lot_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*taille_lot_BDD

                batch = liste_relations[ind_depart:ind_depart+min(taille_lot_BDD, len(liste_relations)-current_transaction*taille_lot_BDD)]
                query = f"""
                        UNWIND $batch AS pair
                        MATCH (a:Noeud {{name: pair.depart}})
                        MATCH (b:Noeud {{name: pair.arrivee}})
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
def creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, dic_noeuds):
    
    nb_transactions = ceil(len(dic_noeuds)/taille_lot_BDD)
    current_transaction = 0
    liste_noeuds = list(dic_noeuds.keys())
    print("Début de création des séquences en BDD")
    with tqdm(total=len(dic_noeuds)) as bar :
        while len(dic_noeuds)-current_transaction*taille_lot_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*taille_lot_BDD
                liste_noeuds_transaction = liste_noeuds[ind_depart:ind_depart+min(taille_lot_BDD, len(dic_noeuds)-current_transaction*10000)]
    
                print("Transaction " + str(current_transaction))
                for t in range(min(taille_lot_BDD, len(dic_noeuds)-current_transaction*taille_lot_BDD)):
                    nc = liste_noeuds_transaction[t]
                    q = """
                        CREATE (n:Sequence {name:$nom, sequence:$sequence})
                        """ 
                    result = tx.run(
                        q,
                        nom=nc,
                        sequence=dic_noeuds[nc]
                    )
                bar.update(min(taille_lot_BDD, len(dic_noeuds)-current_transaction*taille_lot_BDD))
                current_transaction += 1
    print("Fin de création des séquences en BDD")
    liste_kmer = list(dic_kmer_relation.keys())
    print("Début de création des indexes en BDD")
    current_transaction = 0
    with tqdm(total=len(dic_kmer_relation)) as bar :
        while len(dic_kmer_relation)-current_transaction*taille_lot_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*taille_lot_BDD
                liste_kmer_transaction = liste_kmer[ind_depart:ind_depart+min(taille_lot_BDD, len(dic_kmer_relation)-current_transaction*10000)]
    
                print("Transaction " + str(current_transaction))
                for t in range(min(taille_lot_BDD, len(dic_kmer_relation)-current_transaction*taille_lot_BDD)):
                    kmer = liste_kmer_transaction[t]
                    q = """
                        CREATE (k:kmer {kmer:$kmer, taille:$taille})
                        WITH k 
                        MATCH (n:Sequence)
                        WHERE n.name IN $noeuds_cibles
                        CREATE (k)-[:index]->(n)
                        """
                    result = tx.run(
                        q,
                        noeuds_cibles=dic_kmer_relation[kmer],
                        kmer=kmer,
                        taille=kmer_size
                    )
                bar.update(min(taille_lot_BDD, len(dic_kmer_relation)-current_transaction*taille_lot_BDD))
                current_transaction += 1
    print("Fin de création des indexes en BDD")
    return



'''
Cette fonction permet de créer des noeuds avec la séquence du noeud et son nom
Elle créé également des indexes (kmers et noeuds associés)
'''
def charger_sequences_et_indexes(gfa_file_name, kmer_size=31):
    dic_noeuds = {}
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
                        noeud = ligne_dec[1]
                        dic_noeuds[noeud] = seq
                        if len(seq) > kmer_size:
                            liste_kmer = [seq[i:i+kmer_size] for i in range(len(seq)-kmer_size+1)]
                            for kmer in liste_kmer:
                                if kmer in kmer_set :
                                    dic_kmer_relation[kmer].append(noeud)
                                else:
                                    kmer_set.add(kmer)
                                    dic_kmer_relation[kmer] = [noeud]
                ligne = file.readline()
    file.close()
    return dic_kmer_relation, dic_noeuds

'''
Cette fonction permet de créer des noeuds avec la séquence du noeud et son nom
'''
#TODO : découper en lot le chargement des noeuds
def charger_sequences(gfa_file_name, chromosome_file = None, create=False, batch_size=20000000):
    dic_noeuds = {}
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
                                noeud = chromosome_file + "_"+ ligne_dec[1]
                            else : 
                                noeud = ligne_dec[1]
                            dic_noeuds[noeud] = {"sequence":seq}
                    if len(dic_noeuds) >= batch_size:
                        with driver.session() as session:
                            creer_noeuds_batch(session, dic_noeuds, node_name="Sequence", create = create)
                        dic_noeuds = {}
                    ligne = file.readline()
        print("Recuperation des noeud terminé en " + str(time.time()-start_time))
        print("Creation des séquences en BDD : " + str(len(dic_noeuds)) + " noeuds à créer")
        if len(dic_noeuds) > 0 :
            with driver.session() as session:
                creer_noeuds_batch(session, dic_noeuds, node_name="Sequence", create = create)
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
Cette fonction créé les noeuds dans neo4j depuis un fichier GFA
Elle traite les données par lot (maximum fixé par taille_lot) pour éviter les problèmes mémoire
Pour le traitement par lot : on parcourt les noeuds : on sélectionne taille_lot noeuds
Puis on parcourt les chemins pour voir si on trouve les noeuds correspondants, si oui on met à jour le noeud
A la fin des chemins, on créé le lot de noeud en BDD
Les noeuds contiennent les propriétés suivantes : 
   - Taille du noeud
   - Chromosome associé
   - index du noeud pour chaque séquence empruntant le noeud
   - position : la taille cumulée des noeuds précédents pour chaque séquence surle chromosome
   - noeud maitre : en cas de création de sous-noeud pour éviter qu'un noeud soit emprunté plusieurs fois par une meme séquence
                  les sous-noeuds sont notés _2, _3, etc. Par exemple si le noeud maitre est s1, les noeuds secondaires 
                  seront s1_2, s1_3, etc
   - la liste des génomes passant par ce noeud
   - la liste des génomes passant par ce noeud en direct
   - la liste des génomes passant par ce noeud en reverse
   - flux : le pourcentage de génomes passant par ce noeud
L'outil marche avec des P lines ou des W lines, néanmoins, le format des P lines n'étant pas bien défini
En particulier pour le nommage des chromosome et séquence, il faut privilégier les W lines
Pour les W lines : soit le fichier est relatif à un seul chromosome, soit le chromosome est défini dans la 4ème colonne
L'information du chromosome est nécessaire pour pouvoir utiliser correctement l'outil
Input :
    - gfa_file_name : fichier gfa à charger
    - chromosome_file : si le fichier est relatif à un seul chromosome, indiqué le nom du chromosome (chaine de caractère)
    - chromosome_prefix : si True alors le nom du noeud est préfixé par le chromosome (c'est le cas si chromosome_file ets renseigné), si False et pas de chromosome_file alors le noeud n'est pas préfixé
    - taille_lot : le découpage se fait par les noeuds, l'outil va traiter des lots de "taille_lot"
    - noeud_depart : au cas où il y a eu un problème lors d'un chargement, on peut reprendre à un noeud donné
    - create : si c'est la première exécution on peut mettre true, dans ce cas l'outil va créer en BDD sans regarder l'existance
                dans les autres cas il vaut mieux mettre à False
    - haplotype : indique si le nom de l'échantillon doit être concaténer avec l'haplotype

'''
def charger_noeuds_gfa_neo4j(gfa_file_name, chromosome_file = None, chromosome_prefix = False, taille_lot = 5000000, start_chromosome = None, create = False, haplotype = True, create_only_relations = False):
    sep = ["[,;.*]", "(<|>)"]
    nb_lots = 0
    set_genome = set()
    set_chromosome = set()
    total_noeuds = 0
    total_path = 0
    dic_longueur_noeuds = {}
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
                    dic_longueur_noeuds[ligne_dec[1]]=int(len(ligne_dec[2]))
                    total_noeuds += 1
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
                        liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([chaine[:-1] for chaine in liste_noeuds])
                    else:
                        ind = 6
                        walk = 1
                        liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0])
                
            ligne = file.readline() 
        if first_chromosome is not None and first_chromosome != "":
            for k in range(len(chromosomes_list)):
                if chromosomes_list[k] == first_chromosome:
                    index_first_chromosme = k
        print("Nombre de génomes : " + str(len(set_all_genomes)) + " - Liste des chromosomes : " + str(set_all_chromosomes))
        #If create_only_relations is set to True then the nodes are not processed
        if not create_only_relations :
            print("Debut du parsing, nombre de noeuds : " + str(total_noeuds) + "\nChromosome de départ : " + str(start_chromosome))
            for k in range(index_first_chromosme,len(chromosomes_list)) :
                c = chromosomes_list[k]
                nodes_set_chromosome = set(nodes_set_next_chromosome)
                nodes_set_next_chromosome = set()
                print("chromosome " + str(c) + " - number of nodes : " + str(len(nodes_set_chromosome)))
                nb_lots = ceil(len(nodes_set_chromosome)/taille_lot)
                current_lot = 0
                while current_lot < nb_lots :
                    temps_0_lot = time.time()
                    set_noeuds_lot = set(list(nodes_set_chromosome)[current_lot*taille_lot:min(len(nodes_set_chromosome),(current_lot+1)*taille_lot)])
                    current_lot += 1
                    print("chromosome " + c + " lot " + str(current_lot) + "/"+str(nb_lots) + " nombre noeuds : " + str(len(set_noeuds_lot)))
                    #Parcours du fichier pour récupérer la liste des noeuds à traiter pour ce lot
                    file.seek(0,0)
                    ligne = file.readline()
    
                   #Parcours des chemins pour ce lot 
    
                    liste_noeuds = []
                    liste_strand = []
                    compteur_position = {}
                    compteur_noeud = {}
                    dic_noeuds = {}
                    dic_noeuds_maitre = {}
                    set_genomes_lot = set()
                    #dic_longueur_noeuds = {}
    
                    set_noeuds = set()
                    with tqdm(total=total_path) as bar2 :
                        while ligne:
                            ligne_dec = ligne.split()
                            if len(ligne_dec) > 0:
                                #if ligne.startswith(('S')) and ligne_dec[1] in set_noeuds_lot:
                                #        dic_longueur_noeuds[ligne_dec[1]]=int(len(ligne_dec[2]))
                                if ligne[0] == 'P' or ligne[0] == 'W':
                                    chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                                    ligne = None
                                    if current_lot == nb_lots and k < len(chromosomes_list) - 1 and chromosome == chromosomes_list[k+1]:
                                        #Dernier lot du chromosome, on prépare les noeuds à traiter pour le chromosome suivant
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([chaine[:-1] for chaine in liste_noeuds])
                                        else:
                                            ind = 6
                                            walk = 1
                                            liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0])
                                    
                                    if chromosome == c:
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [chaine[-1] for chaine in liste_noeuds]
                                            liste_noeuds = [chaine[:-1] for chaine in liste_noeuds]
                                            # if chromosome_file == None :
                                            #     liste_strand = [chaine[-1] for chaine in liste_noeuds]
                                            #     liste_noeuds = [chaine[:-1] for chaine in liste_noeuds]
                                            # else :
                                            #     liste_strand = [chaine[-1] for chaine in liste_noeuds]
                                            #     liste_noeuds = [chromosome_file+"_"+chaine[:-1] for chaine in liste_noeuds]                              
                                        else:
                                            ind = 6
                                            walk = 1
                                            liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 != 0]
                                            liste_noeuds = [liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0]
                                            # if chromosome_file == None :
                                            #     liste_noeuds = [liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0]
                                            # else:
                                            #     liste_noeuds = [chromosome_file+"_"+liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0]
          
                                        if chromosome is not None and chromosome not in set_chromosome:
                                            set_chromosome.add(chromosome)
                                        if genome != "_MINIGRAPH_":
                                            if genome not in set_genome:
                                                set_genome.add(genome)
                                            if genome not in set_genomes_lot :
                                                set_genomes_lot.add(genome)
                                                compteur_noeud[genome] = {}
                                                compteur_position[genome] = {}
                                            if chromosome not in compteur_noeud[genome] :
                                                compteur_noeud[genome][chromosome] = 0
                                                compteur_position[genome][chromosome] = {}
                                                compteur_position[genome][chromosome]["current_position"] = 0
                                                compteur_position[genome][chromosome]["previous_position"] = 0
                                                compteur_position[genome][chromosome]["current_contig"] = ""
                                            #Dans le cas du Walk, la position de départ est indiquée => on l'utilise
                                            if ind == 6:
                                                if compteur_position[genome][chromosome]["current_contig"] != ligne_dec[3] :
                                                    #On change de contig, on va ajouter le départ du contig suivant
                                                    compteur_position[genome][chromosome]["current_position"] += int(ligne_dec[4])
                                                else :
                                                    #On est dans le même contig, on ajoute les éventuels gaps
                                                    if compteur_position[genome][chromosome]["previous_position"] - int(ligne_dec[4]) > 0 :
                                                        compteur_position[genome][chromosome]["current_position"] += compteur_position[genome][chromosome]["previous_position"] - int(ligne_dec[4])
                                                compteur_position[genome][chromosome]["current_contig"] = ligne_dec[3]
                                                compteur_position[genome][chromosome]["previous_position"] = int(ligne_dec[5])
        
                                            noeud = ""
                                            noeud_maitre = ""
                                            strand = ""
                                            ligne_dec=None
                                            #Linéarisation du graphe
                                            # Parcours de la liste des noeuds : si un noeud existe déjà dans un chemin pour la même séquence
                                            # Alors on va créer un nouveau noeud (par exemple si c'est la 6ème itération pour le noeud S1 on va créer S1_6)
                                            for i in range(len(liste_noeuds)):
                                                if i > 0:
                                                    noeud_precedent = noeud
                                                noeud = liste_noeuds[i]
                                                noeud_maitre = noeud
                                                taille = dic_longueur_noeuds[noeud_maitre]
                                                if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                                    noeud = chromosome + "_" + noeud
                                                #On ne traite le noeud que s'il fait parti du lot
                                                
                                                
                                                if noeud_maitre in set_noeuds_lot:
                                                    if chromosome_file is not None and chromosome_file != "":
                                                        noeud_maitre = chromosome + "_" + noeud_maitre
                                                    strand = ""
                                                    if (i < len(liste_strand) and liste_strand[i] in ["-", "<"]):
                                                        strand = "M"
                                                    else:
                                                        strand = "P"
                                                    
                                                    if noeud not in set_noeuds:
                                                        if strand == "P" :
                                                            dic_noeuds[noeud] = {"genomes":[genome], "max":1, "strandP":[genome], "strandM":[], "ref_node" : noeud_maitre, genome+"_noeud":compteur_noeud[genome][chromosome],genome+"_position":compteur_position[genome][chromosome]["current_position"], "taille" : taille, "chromosome"  : chromosome, "position_min":compteur_position[genome][chromosome]["current_position"], "position_max":compteur_position[genome][chromosome]["current_position"]}
                                                        else :
                                                            dic_noeuds[noeud] = {"genomes":[genome], "max":1, "strandM":[genome], "strandP":[], "ref_node" : noeud_maitre, genome+"_noeud":compteur_noeud[genome][chromosome],genome+"_position":compteur_position[genome][chromosome]["current_position"], "taille" : taille, "chromosome"  : chromosome, "position_min":compteur_position[genome][chromosome]["current_position"], "position_max":compteur_position[genome][chromosome]["current_position"]}
                                                        set_noeuds.add(noeud)
                                                    else:
                                                        if genome not in dic_noeuds[noeud]["genomes"] and chromosome == dic_noeuds[noeud]["chromosome"]:
                                                            dic_noeuds[noeud]["genomes"].append(genome)
                                                            dic_noeuds[noeud]["strand"+strand].append(genome)
                                                            dic_noeuds[noeud][genome+"_noeud"] = compteur_noeud[genome][chromosome]
                                                            dic_noeuds[noeud][genome+"_position"] = compteur_position[genome][chromosome]["current_position"]
                                                            if compteur_position[genome][chromosome]["current_position"] < dic_noeuds[noeud]["position_min"]:
                                                                dic_noeuds[noeud]["position_min"] = compteur_position[genome][chromosome]["current_position"]
                                                            if compteur_position[genome][chromosome]["current_position"] > dic_noeuds[noeud]["position_max"]:
                                                                dic_noeuds[noeud]["position_max"] = compteur_position[genome][chromosome]["current_position"]
                                                        else :
                                                            #Le noeud est redondant, on va chercher si un noeud de même séquence
                                                            #est disponible, sinon on créé un nouveau noeud
                                                            if noeud_maitre not in dic_noeuds_maitre:
                                                                dic_noeuds_maitre[noeud_maitre] = {}
                                                            if genome+"-"+chromosome not in dic_noeuds_maitre[noeud_maitre] :
                                                                dic_noeuds_maitre[noeud_maitre][genome+"-"+chromosome] = 2
                                                            else:
                                                                dic_noeuds_maitre[noeud_maitre][genome+"-"+chromosome] += 1
                                                            noeud = noeud + "_" + str(dic_noeuds_maitre[noeud_maitre][genome+"-"+chromosome])
                                                            if noeud in set_noeuds:
                                                                dic_noeuds[noeud]["genomes"].append(genome)
                                                                dic_noeuds[noeud]["strand"+strand].append(genome)
                                                                dic_noeuds[noeud][genome+"_noeud"] = compteur_noeud[genome][chromosome]
                                                                dic_noeuds[noeud][genome+"_position"] = compteur_position[genome][chromosome]["current_position"]
                                                                if compteur_position[genome][chromosome]["current_position"] < dic_noeuds[noeud]["position_min"]:
                                                                    dic_noeuds[noeud]["position_min"] = compteur_position[genome][chromosome]["current_position"]
                                                                if compteur_position[genome][chromosome]["current_position"] > dic_noeuds[noeud]["position_max"]:
                                                                    dic_noeuds[noeud]["position_max"] = compteur_position[genome][chromosome]["current_position"]
                                                            else:
                                                                if strand == "P" :
                                                                    dic_noeuds[noeud] = {"genomes":[genome], "max":1, "ref_node" : noeud_maitre, "strandP":[genome], "strandM":[], "taille" : taille, "chromosome"  : chromosome, "position_min":compteur_position[genome][chromosome]["current_position"], "position_max":compteur_position[genome][chromosome]["current_position"]}
                                                                else :
                                                                    dic_noeuds[noeud] = {"genomes":[genome], "max":1, "ref_node" : noeud_maitre, "strandM":[genome], "strandP":[], "taille" : taille, "chromosome"  : chromosome, "position_min":compteur_position[genome][chromosome]["current_position"], "position_max":compteur_position[genome][chromosome]["current_position"]}
                                                                dic_noeuds[noeud][genome+"_noeud"] = compteur_noeud[genome][chromosome]   
                                                                dic_noeuds[noeud][genome+"_position"] = compteur_position[genome][chromosome]["current_position"]   
                                                                dic_noeuds[noeud]["max"]+=1
                                                                set_noeuds.add(noeud)
                                                                
                                                compteur_noeud[genome][chromosome] += 1
                                                compteur_position[genome][chromosome]["current_position"] += taille
                                    bar2.update(1)
                            ligne = file.readline() 
                    liste_noeuds = None
                    #Calcul des flux
                    print("\nTaille des éléments à créer en BDD : " + str(len(list(dic_noeuds.items()))))
                    if len(dic_noeuds) > 0:
                        for noeud in dic_noeuds:
                            dic_noeuds[noeud]["flux"] = len(dic_noeuds[noeud]["genomes"])/len(set_all_genomes)
                            moyenne_noeud = 0
                            moyenne_position = 0
                            nb_genomes = 0
                            for g in dic_noeuds[noeud]["genomes"]:
                                nb_genomes += 1
                                moyenne_noeud += dic_noeuds[noeud][g+"_noeud"]
                                moyenne_position += dic_noeuds[noeud][g+"_position"]
                            dic_noeuds[noeud]["moyenne_position"] = int(moyenne_position/nb_genomes)
                            dic_noeuds[noeud]["moyenne_noeud"] = int(moyenne_noeud/nb_genomes)
                    
                        print("Durée analyse des noeuds du lots : "+ str(time.time()-temps_0_lot) + "\nDurée totale écoulée : " + str(time.time()-temps_depart))
                        #Création du lot de noeud en BDD
                        print("Lot " + str(current_lot) + " : création des noeuds en BDD")
                        with get_driver() as driver:
                            print("Connection established.")
                            with driver.session() as session:
                                creer_noeuds_batch(session, dic_noeuds, create=create)
                        dic_noeuds = None
                        dic_noeuds_maitre = None
                        
                        
            creer_stats(set_genome, set_chromosome)
            dic_longueur_noeuds = None
            print("Fin de création des noeuds\nDurée du traitement : " + str(time.time()-temps_depart) + "\nGénomes traités : " + str(set_genome) + "\nNb de noeuds traités : "+str(total_noeuds) )
        
        #Les relations sont traitées après car le découpage en noeud ne permet pas de garantir 
        #qu'on va créer les relations entre les bons noeuds en cas de noeuds doublonnés
        #Pour les relations on va donc traiter par chromosome
        print("\nDébut du traitement des relations")
        
        set_relations = set()
        liste_relations = []
        liste_noeuds = []
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
                                    liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        liste_strand = [chaine[-1] for chaine in liste_noeuds]
                                        liste_noeuds = [chromosome_file+"_"+chaine[:-1] for chaine in liste_noeuds]   
                                    else :
                                        liste_strand = [chaine[-1] for chaine in liste_noeuds]
                                        liste_noeuds = [chaine[:-1] for chaine in liste_noeuds]                           
                                else:
                                    ind = 6
                                    walk = 1
                                    liste_noeuds = re.split(sep[walk], ligne_dec[ind])
                                    liste_strand = [liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 != 0]
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        liste_noeuds = [chromosome_file+"_"+liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0]
                                    else:
                                        liste_noeuds = [liste_noeuds[j] for j in range(1, len(liste_noeuds)) if j % 2 == 0]
        
        
                                if genome != "_MINIGRAPH_":
                                    noeud = ""
                                    noeud_maitre = ""
                                    strand = ""
                                    # Parcours de la liste des noeuds : si un noeud existe déjà dans un chemin pour la même séquence
                                    # Alors on va créer un nouveau noeud (par exemple si c'est la 6ème itération pour le noeud S1 on va créer S1_6)
                                    for i in range(len(liste_noeuds)):
                                        if i > 0:
                                            noeud_precedent = noeud
                                        noeud = liste_noeuds[i]
                                        if noeud in repeat_nodes :
                                            if genome+"-"+chromosome in repeat_nodes[noeud]:
                                                repeat_nodes[noeud][genome+"-"+chromosome] += 1
                                                noeud = noeud+"_"+ str(repeat_nodes[noeud][genome+"-"+chromosome])
                                            else:
                                                repeat_nodes[noeud][genome+"-"+chromosome] = 1
                                        else:
                                            repeat_nodes[noeud] = {genome+"-"+chromosome:1}
                                        if i > 0 and noeud_precedent+"->"+noeud not in set_relations :
                                            set_relations.add(noeud_precedent+"->"+noeud)
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
    dic_noeuds = {}
    file_format = "gtf"
    with file:
        n = 0
        for line in file :
            n +=1
        file.seek(0,0)
        with tqdm(total=n) as bar :
            #Création des noeuds annotations
            ligne = file.readline()
            while ligne :
                if ligne[0] != '#':
                    ligne_dec = ligne.split()
                    chromosome = get_chromosome_annotation(ligne_dec[0])
                    if single_chromosome == None or single_chromosome == chromosome : 
                        name = hashlib.sha256(ligne.encode("utf-8")).hexdigest()
                        dic_noeuds[name] = {}
                        
                        feature = ligne_dec[2].lower()
                        dic_noeuds[name]["name"] = name
                        dic_noeuds[name]["chromosome"] = chromosome
                        dic_noeuds[name]["genome_ref"] = genome_ref
                        dic_noeuds[name]["source"] = ligne_dec[1]
                        dic_noeuds[name]["feature"] = feature
                        
                        start = int(ligne_dec[3])
                        end = int(ligne_dec[4])
                        dic_noeuds[name]["start"] = start
                        dic_noeuds[name]["end"] = end
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
                                        dic_noeuds[name]["gene_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_version":
                                        dic_noeuds[name]["gene_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_id":
                                        dic_noeuds[name]["transcript_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_version":
                                        dic_noeuds[name]["transcript_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "exon_number":
                                        dic_noeuds[name]["exon_number"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_name":
                                        dic_noeuds[name]["gene_name"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_source":
                                        dic_noeuds[name]["gene_source"] = attributes[i+1][:-1].replace('"',"")    
                                    case "gene_biotype":
                                        dic_noeuds[name]["gene_biotype"] = attributes[i+1][:-1].replace('"',"")     
                                    case "transcript_name":
                                        dic_noeuds[name]["transcript_name"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_source":
                                        dic_noeuds[name]["transcript_source"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_biotype":
                                        dic_noeuds[name]["transcript_biotype"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_id":
                                        dic_noeuds[name]["protein_id"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_version":
                                        dic_noeuds[name]["protein_version"] = attributes[i+1][:-1].replace('"',"") 
                                    case "tag":
                                        dic_noeuds[name]["tag"] = attributes[i+1][:-1].replace('"',"") 
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
                                    dic_noeuds[name]["exon_id"] = exon_id
                                    dic_noeuds[name]["exon_number"] = exon_number 
                                else :
                                    if current_attribute == "id":
                                        dic_noeuds[name]["id"] = attributes[i+1]
                                    elif current_attribute == "name":
                                        dic_noeuds[name]["name"] = attributes[i+1]
                                        
                            if current_gene_id != "" :
                                dic_noeuds[name]["gene_id"] = current_gene_id
                                dic_noeuds[name]["gene_name"] = current_gene_name
                                if current_transcript_id != "":
                                    dic_noeuds[name]["transcript_id"] = current_transcript_id
                                    dic_noeuds[name]["transcript_name"] = current_transcript_name
                            
                                    
                            
                
                bar.update(1)
                ligne = file.readline()

                    
           
            print("\nTaille des éléments à créer en BDD : " + str(len(list(dic_noeuds.items()))))

            print("Durée analyse des noeuds du lots : "+ str(time.time()-temps_depart))


            with get_driver() as driver:
                print("Connection established.")
                with driver.session() as session:
                    creer_noeuds_batch(session, dic_noeuds, node_name=node_name)

        print("Fin de création des noeuds\nDurée du traitement : " + str(time.time()-temps_depart))
        #Fin des lots, on créé toutes les relations                                               
    file.close()
    

    print("Chargement annotations terminée\nDurée totale du traitement : "+ str(time.time()-temps_depart))


def process_annotation_simple_batch(tx, annotations, genome_ref):
    #for a in annotations:
        #print(a)
    tx.run("""
        UNWIND $annotations AS annot
        MATCH (a:Annotation {name: annot.name})
        MATCH (n:Noeud)
        WHERE n.chromosome = annot.chromosome 
        AND n.`"""+str(genome_ref)+"""_position` >= annot.start 
        AND n.`"""+str(genome_ref)+"""_position` <= annot.end
        MERGE (n)-[:A_POUR_ANNOTATION]->(a) """, annotations=annotations)

def process_annotation_complexe_batch(tx, last_id, genome_ref, batch_size = 100000):
    
    result = tx.run(
        """
        MATCH (n:Noeud)
        WHERE id(n) > $last_id
        WITH n ORDER BY id(n) ASC LIMIT $limit
        OPTIONAL MATCH (a:Annotation)
        WHERE a.chromosome = n.chromosome
          AND a.start > n.`"""+str(genome_ref)+"""_position`
          AND a.start < n.`"""+str(genome_ref)+"""_position` + n.taille
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
            result = session.run("MATCH (n:Noeud) RETURN max(id(n)) AS max_id")
            max_id = result.single()["max_id"]
            print("max id : " + str(max_id))
            while last_id is not None and last_id < max_id :

                created = session.execute_write(process_annotation_complexe_batch, last_id, genome_ref, batch_size)
                
                result = session.run(
                    """
                    MATCH (n:Noeud)
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
def construct_DB(gfa_file_name, annotation_file_name = None, genome_ref = None, chromosome_file = None, noeud_depart = 0, taille_lot = 5000000, create=False):
    start_time = time.time()
    charger_sequences(gfa_file_name, chromosome_file, create=create)
    sequence_time = time.time()
    print("Sequences loaded in " + str(sequence_time-start_time) + " s")
    charger_noeuds_gfa_neo4j(gfa_file_name, chromosome_file = chromosome_file, taille_lot = taille_lot, noeud_depart = noeud_depart, create = create)
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
def construct_db_by_chromosome(gfa_chromosomes_dir, annotation_file_name = None, genome_ref = None, chromosome_file = None, noeud_depart = 0, taille_lot = 5000000, create=False):
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
                charger_noeuds_gfa_neo4j(gfa_file_name, chromosome_file = chromosome_file, taille_lot = taille_lot, noeud_depart = noeud_depart, create = create)
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
    dic_kmer_relation, dic_noeuds = charger_sequences_et_indexes(gfa_file_name, kmer_size=31)
    with get_driver() as driver:
        print("Connection established.")
        with driver.session() as session:
            creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, dic_noeuds)
        
def launcher():
    for c in ["18", "19", "20", "21", "22", "X", "Y"]:
        print("chargement du chromosome : " + str(c))
        set_genomes = charger_noeuds_gfa_neo4j("/home/fgraziani/work/project/Pangenomique/GFA/humain/hprc-mc/chr"+str(c)+".gfa", chromosome_file = str(c), taille_lot=1000000, create=False)
    