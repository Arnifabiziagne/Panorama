#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:45:39 2025

@author: fgraziani
"""



import os
import json
import subprocess
import time
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable


# global variables with default value
DB_URL="bolt://localhost:7687"
AUTH=("neo4j", "Administrateur")

#Delay in s to wait between tests of the database connexion
neo4j_start_delay = 10
#max_iter_test_connexion is the maximum number of pooling database to test if it is up
max_iter_test_connexion = 60


def test_connection(DB_URL, AUTH):
    try:
        driver = GraphDatabase.driver(DB_URL, auth=AUTH)
        with driver.session() as session:
            session.run("RETURN 1")
        driver.close()
        return True
    except ServiceUnavailable:
        return False


def load_config_from_json():
    config_path = os.path.join(os.path.dirname(__file__), "db_conf.json")
    if not os.path.exists(config_path):
        return None

    with open(config_path, "r") as f:
        return json.load(f)


def start_neo4j_container(container_name, DB_URL, AUTH):
    print(f"Launching Neo4j docker container : {container_name}")
    # Lancer le container s'il n'existe pas déjà
    try:
        subprocess.run(["docker", "start", container_name], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Fail to launch container : {e}")
        return False

    # Waiting for neo4j up
    time.sleep(neo4j_start_delay)
    i = 0
    while not test_connection(DB_URL, AUTH) and i < max_iter_test_connexion:
        i += 1
        print(f"⏳ Trying to connect to database – attempt {i}")
        time.sleep(neo4j_start_delay)

    if i < max_iter_test_connexion:
        print("✅ Database is up.")
        return True
    else:
        print("❌ Database did not start in time.")
        return False

def stop_container():
    conf = load_config_from_json()
    if not conf:
        print("⚠️ db_conf.json not present, could not connect to neo4j DB.")
        return None
    container_name = conf.get("container_name")
    try:
        subprocess.run(["docker", "stop", container_name], check=True)
        time.sleep(10)
        print("container stopped")
    except subprocess.CalledProcessError as e:
        print(f"Fail to stop container : {e}")
        return False

        
        
def get_driver():
    global DB_URL, AUTH


    if test_connection(DB_URL, AUTH):
        driver = GraphDatabase.driver(DB_URL, auth=AUTH)
        return driver

    conf = load_config_from_json()
    if not conf:
        print("⚠️ db_conf.json not present, could not connect to neo4j DB.")
        return None
    
    bolt_port = conf.get("bolt_port")
    DB_URL = f"bolt://localhost:{bolt_port}"
    AUTH = (conf.get("login"),conf.get("password"))
    container_name = conf.get("container_name")
    
    
    
    
    if test_connection(DB_URL, AUTH):
        driver = GraphDatabase.driver(DB_URL, auth=AUTH)
        return driver

    if start_neo4j_container(container_name, DB_URL, AUTH):
        if test_connection(DB_URL, AUTH):
            driver = GraphDatabase.driver(DB_URL, auth=AUTH)
            return driver
        else:
            print("❌ Connexion to neo4j db failed.")
            return None
    else:
        print("❌ Fail to launch container.")
        return None