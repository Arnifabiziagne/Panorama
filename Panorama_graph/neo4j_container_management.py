#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 10:21:02 2025

@author: fgraziani
"""

import os
import shutil
import subprocess
import json
import time

from config import *
from neo4j_DB_construction import create_stats_from_nodes, create_indexes

# --- CONSTANTES ---
DOCKER_IMAGE = "neo4j:2025.05-community-bullseye"
HTTP_PORT = 7474
BOLT_PORT = 7687
NEO4J_BASE_DIR = os.path.abspath("./data")
CONF_SOURCE_FILE = os.path.abspath("./neo4j_install/conf/neo4j.conf")
CONF_FILE = os.path.abspath("./db_conf.json")
IMPORT_DIR = os.path.abspath("./data/import")
DUMP_FILE = os.path.join(IMPORT_DIR, "neo4j.dump")
#Read buffer size is important in case of csv import : 
# it limits the line's size (by default 4 * 1024 * 1024) : for long nodes sequences it could be necessary to improve this value
READ_BUFFER_SIZE = 33554432

MAX_MEM = "24g"
MAX_SWAP = "25g"
MAX_CPU = "8"

NEO4J_AUTH = "neo4j/Administrateur"
NEO4J_LOGIN = "neo4j"
NEO4J_PASSWORD = "Administrateur"


def prepare_data_directories_in_container():
    """
    Ensure /data/databases/neo4j exists inside the container (run as root to avoid permission issues).
    """
    print("🛠️ Preparing data directories inside container...")
    subprocess.run([
        "docker", "run", "--rm",
        "--user=root",  # 👈 important
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        DOCKER_IMAGE,
        "bash", "-c", "mkdir -p /data/databases/neo4j"
    ], check=True)

def remove_directories():
    for folder in ["data", "logs", "conf", "plugins"]:
        path = os.path.join(NEO4J_BASE_DIR, folder)
        if os.path.exists(path):
            shutil.rmtree(path)

def import_dump():
    print("📂 Importing dump...")
    prepare_data_directories_in_container()
    subprocess.run([
        "docker", "run", "--rm",
        f"--cpus={MAX_CPU}",
        f"--memory={MAX_MEM}",
        f"--memory-swap={MAX_SWAP}",
        "-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
        "-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
        "-v", f"{IMPORT_DIR}:/import",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        DOCKER_IMAGE,
        "neo4j-admin", "database", "load", "neo4j",
        "--from-path=/import", "--overwrite-destination=true"
    ], check=True)

def import_csv():
    print(f"📂 Importing CSV - data dir : {NEO4J_BASE_DIR}/data ...")
    prepare_data_directories_in_container()
    subprocess.run([
        "docker", "run", "--rm",
        f"--cpus={MAX_CPU}",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        "-v", f"{IMPORT_DIR}:/import",
        "-e", f"NEO4J_AUTH={NEO4J_AUTH}",
        DOCKER_IMAGE,
        "neo4j-admin", "database", "import", "full",
        "--verbose",
        f"--read-buffer-size={READ_BUFFER_SIZE}",
        f"--max-off-heap-memory={MAX_MEM}",
        "--nodes=Sequence=/import/sequences.csv",
        "--nodes=Node=/import/nodes.csv",
        "--relationships=/import/relations.csv"
    ], check=True)

def start_container(container_name):
    print("🚀 Starting Neo4j container...")
    subprocess.run([
        "docker", "run", "-d",
        "--name", container_name,
        "-e", f"NEO4J_AUTH={NEO4J_AUTH}",
        "-e", "NEO4J_ACCEPT_LICENSE_AGREEMENT=yes",
        "-e", "NEO4J_apoc_export_file_enabled=true",
        "-e", "NEO4J_apoc_import_file_enabled=true",
        "-e", "NEO4J_apoc_import_file_use__neo4j__config=true",
        "-e", "NEO4J_PLUGINS=[\"apoc\"]",
        "-p", f"{HTTP_PORT}:7474",
        "-p", f"{BOLT_PORT}:7687",
        "-v", f"{NEO4J_BASE_DIR}/data:/data",
        "-v", f"{NEO4J_BASE_DIR}/logs:/logs",
        "-v", f"{NEO4J_BASE_DIR}/conf:/conf",
        "-v", f"{NEO4J_BASE_DIR}/import:/import",
        "-v", f"{NEO4J_BASE_DIR}/plugins:/plugins",
        DOCKER_IMAGE
    ], check=True)

    time.sleep(10)
    print(f"✅ Neo4j {container_name} is ready!")
    print(f"🌍 HTTP: http://localhost:{HTTP_PORT}")
    print(f"🔗 BOLT: bolt://localhost:{BOLT_PORT}")

def write_config(container_name):
    with open(CONF_FILE, "w") as f:
        json.dump({
            "container_name": container_name,
            "http_port": HTTP_PORT,
            "bolt_port": BOLT_PORT,
            "login": NEO4J_LOGIN,
            "password": NEO4J_PASSWORD
        }, f, indent=2)


def stop_container(container_name: str):
    """
    Stops and removes a Docker container if it exists.
    Does nothing if the container does not exist.
    """
    try:
        # List all container names (including stopped ones)
        result = subprocess.run(
            ["docker", "ps", "-a", "--format", "{{.Names}}"],
            capture_output=True, text=True, check=True
        )
        existing_containers = result.stdout.strip().splitlines()

        if container_name in existing_containers:
            print(f"🛑 Stopping existing container: {container_name}")
            subprocess.run(["docker", "stop", container_name], check=True)
        else:
            print(f"ℹ️ Container '{container_name}' does not exist. Nothing to do.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while stopping/removing container: {e}")


def remove_container(container_name: str):
    """
    Stops and removes a Docker container if it exists.
    Does nothing if the container does not exist.
    """
    try:
        # List all container names (including stopped ones)
        result = subprocess.run(
            ["docker", "ps", "-a", "--format", "{{.Names}}"],
            capture_output=True, text=True, check=True
        )
        existing_containers = result.stdout.strip().splitlines()

        if container_name in existing_containers:
            print(f"🛑 Stopping and removing existing container: {container_name}")
            subprocess.run(["docker", "rm", "-f", container_name], check=True)
        else:
            print(f"ℹ️ Container '{container_name}' does not exist. Nothing to do.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while stopping/removing container: {e}")

#This function create a docker neo4j database
#If a dump file exists in ./data/import directory it will load the data into database
#If no dump but nodes.csv, sequences.csv and relations.csv exists in ./data/import directory it will load these data into database
#In other cases it will create an empty database
def create_db(container_name, docker_image=DOCKER_IMAGE):
    csv_import_mode = False
    
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    print(f"🔧 Creating database for container '{container_name}' with {DOCKER_IMAGE} neo4j image")
    
    
    
    # Stop container
    remove_container(container_name)
    
    data_db_dir = os.path.join(NEO4J_BASE_DIR, "data", "databases", "neo4j")
    csv_nodes = os.path.join(IMPORT_DIR, "nodes.csv")
    csv_relations = os.path.join(IMPORT_DIR, "relations.csv")
    csv_sequences = os.path.join(IMPORT_DIR, "sequences.csv")
    
    if os.path.exists(data_db_dir) and os.listdir(data_db_dir):
        remove_directories()
    for d in ["data", "logs", "conf", "import", "plugins", "gfa", "annotations", "data/databases", "data/databases/neo4j"]:
        os.makedirs(os.path.join(NEO4J_BASE_DIR, d), exist_ok=True)
   
    # copy conf file
    if os.path.isfile(CONF_SOURCE_FILE):
        shutil.copy(CONF_SOURCE_FILE, os.path.join(NEO4J_BASE_DIR, "conf"))
        print("🔧 Config file copied")
    else:
        print(f"⚠️ Config file {CONF_SOURCE_FILE} not found")
    
    # --- Import via dump ---
    if os.path.isfile(DUMP_FILE):
        print("📥 Detected dump file for import")
        import_dump()
    
    # --- Import via CSV ---
    elif os.path.isfile(csv_nodes) and os.path.isfile(csv_relations) and os.path.isfile(csv_sequences):
        print("📥 Detected CSV files for import")
        import_csv()
        csv_import_mode = True
    # --- New databse creation --- #
    else :
        remove_directories()
        # Directories creation
        for d in ["data", "logs", "conf", "import", "plugins", "gfa", "annotations"]:
            os.makedirs(os.path.join(NEO4J_BASE_DIR, d), exist_ok=True)
    
    
    # Save container conf in db_conj.json
    print(f"writing conf for container name : {container_name}")
    write_config(container_name)
    
    # Launch the container
    start_container(container_name)

    if csv_import_mode:
        print("creating stats")
        create_stats_from_nodes()
        print("creating indexes")
        create_indexes(base=True, extend=True, genomes_index=True)

    
    
def dump_db(container_name, docker_image=DOCKER_IMAGE):
    print("📦 Creating Neo4j dump...")
    if docker_image is not None :
        DOCKER_IMAGE = docker_image
    # Stop container
    stop_container(container_name)

    try:
        subprocess.run([
            "docker", "run", "--rm",
            f"--cpus={MAX_CPU}",
            f"--memory={MAX_MEM}",
            f"--memory-swap={MAX_SWAP}",
            "-e", f"JAVA_OPTS=-Xmx{MAX_MEM} -Xms1g",
            "-e", f"NEO4J_dbms.memory.heap.max_size={MAX_MEM}",
            "-v", f"{NEO4J_BASE_DIR}/data:/data",
            "-v", f"{IMPORT_DIR}:/import",
            DOCKER_IMAGE,
            "neo4j-admin", "database", "dump", "neo4j",
            f"--to-path=/import"
        ], check=True)

        print(f"✅ Dump successfully created at: {os.path.join(IMPORT_DIR, 'neo4j.dump')}")
        start_container(container_name)

    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to create dump: {e}")