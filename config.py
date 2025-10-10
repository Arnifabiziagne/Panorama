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
from neo4j_container_management import *


CONF_FILE = os.path.abspath("./db_conf.json")

#Delay in s to wait between tests of the database connexion
neo4j_start_delay = 10
#max_iter_test_connexion is the maximum number of pooling database to test if it is up
max_iter_test_connexion = 60


def get_conf():
    with open(CONF_FILE) as f:
        conf = json.load(f)
    container_name = str(conf["container_name"])   
    HTTP_PORT = int(conf["http_port"])
    BOLT_PORT = int(conf["bolt_port"])
    AUTH = (conf["login"],conf["password"])
    DB_URL="bolt://localhost:"+str(BOLT_PORT)
    return {"container_name":container_name, "HTTP_PORT":HTTP_PORT, "BOLT_PORT":BOLT_PORT, "AUTH":AUTH, "DB_URL":DB_URL}

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

        
        
def get_driver(max_retries=5, retry_delay=5):
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        DB_URL = conf["DB_URL"]
        AUTH= conf["AUTH"]
        for attempt in range(1, max_retries + 1):
            if test_connection(DB_URL, AUTH):
                driver = GraphDatabase.driver(DB_URL, auth=AUTH)
                return driver
            else:
                print(f"Retry {attempt}/{max_retries} to connect to driver")
                time.sleep(retry_delay)
                
    
        if test_connection(DB_URL, AUTH):
                driver = GraphDatabase.driver(DB_URL, auth=AUTH)
                return driver
        print(f"❌ Fail to access Database of container {conf['container_name']}.")
        return None
    else:
        return None