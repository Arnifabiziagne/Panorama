#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:45:39 2025

@author: fgraziani
"""



import os
import json
import shutil
import subprocess
import time
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import logging


logger = logging.getLogger("panorama_logger")


CONF_FILE = os.path.abspath("./conf.json")
OLD_CONF_FILE = os.path.abspath("./db_conf.json")


#Delay in s to wait between tests of the database connexion
neo4j_start_delay = 10
#max_iter_test_connexion is the maximum number of pooling database to test if it is up
max_iter_test_connexion = 60


#Default value to add to the old config file
DEFAULT_ADDITIONS = {
    "server_mode": False,
    "admin_mode": False,
    "admin_users": {"admin": "1234"},
    "server_log_mode": "both",
    "log_retention_days": 7,
    "log_level": "DEBUG"
}

#Function to retrocompatibility with old conf file
def check_conf_file():
    """Check if the main config file exists; otherwise copy the old one,
    update it with default fields, and delete the old config file."""

    # 1️⃣ If the main config file already exists → nothing to do
    if os.path.exists(CONF_FILE):
        return
    
    # 2️⃣ If the old config file exists → copy and upgrade it
    if os.path.exists(OLD_CONF_FILE):
        logger.info(f"⚙️ Copying {OLD_CONF_FILE} to {CONF_FILE} ...")
        shutil.copyfile(OLD_CONF_FILE, CONF_FILE)
    
        # Try to load the JSON content
        try:
            with open(CONF_FILE, "r", encoding="utf-8") as f:
                conf_data = json.load(f)
        except json.JSONDecodeError:
            logger.error("❌ Error: old config file is invalid JSON. Creating a new one.")
            conf_data = {}
    
        # Add or update the required fields
        conf_data.update(DEFAULT_ADDITIONS)
    
        # Save the updated config
        with open(CONF_FILE, "w", encoding="utf-8") as f:
            json.dump(conf_data, f, indent=4, ensure_ascii=False)
        logger.info(f"✅ New config file created: {CONF_FILE}")
    
        # Delete the old config file
        try:
            os.remove(OLD_CONF_FILE)
            logger.info(f"🗑️  Old config file removed: {OLD_CONF_FILE}")
        except Exception as e:
            logger.error(f"⚠️ Could not remove old config file: {e}")
        return

check_conf_file()    

def get_conf(log_levels=["INFO", "DEBUG", "WARNING","ERROR", "CRITICAL", "NOTSET"]):
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    with open(CONF_FILE) as f:
        conf = json.load(f)
    container_name = str(conf["container_name"])   
    HTTP_PORT = int(conf["http_port"])
    BOLT_PORT = int(conf["bolt_port"])
    AUTH = (conf["login"],conf["password"])
    DB_URL="bolt://localhost:"+str(BOLT_PORT)
    SERVER_MODE = bool(conf.get("server_mode",False))
    ADMIN_MODE = bool(conf.get("admin_mode",False))
    users_list=conf.get("admin_users", {})
    LOG_RETENTION_DAYS = int(conf.get("log_retention_days",7))
    SERVER_LOG_MODE = str(conf.get("server_log_mode","both"))
    LOG_LEVEL_PARAM = str(conf.get("log_level","INFO")).upper()
    if LOG_LEVEL_PARAM in log_levels:
        LOG_LEVEL= "logging."+LOG_LEVEL_PARAM
    else:
        LOG_LEVEL= "logging.INFO"
    USERS = {}
    if len(users_list) > 0:
        for k, v in users_list.items():
            if not isinstance(k, str) or not isinstance(v, str):
                raise ValueError("Password must be string")
            if k.strip() == "" or v.strip() == "":
                raise ValueError("Usernames and password can't be empty")
            USERS[k.strip()] = v.strip()
    return {"container_name":container_name, "HTTP_PORT":HTTP_PORT, "BOLT_PORT":BOLT_PORT, 
            "AUTH":AUTH, "DB_URL":DB_URL, "SERVER_MODE":SERVER_MODE, "USERS":USERS, "ADMIN_MODE":ADMIN_MODE,
            "SERVER_LOG_MODE":SERVER_LOG_MODE,"LOG_RETENTION_DAYS":LOG_RETENTION_DAYS, "LOG_LEVEL":LOG_LEVEL}

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
    if not os.path.exists(CONF_FILE):
        return None

    with open(CONF_FILE, "r") as f:
        return json.load(f)

def is_server_mode():
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    server_mode = False
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        server_mode = conf["SERVER_MODE"]
    return server_mode


def is_admin_mode():
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    admin_mode = False
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        admin_mode = conf["ADMIN_MODE"]
    return admin_mode



#Authorization is set to True for local installation of panorama
#or if the mode admin is set to True
def check_authorization():
    authorization = True
    if is_server_mode() and not is_admin_mode():
        authorization = False
    return authorization
    
def get_users(): 
     if not os.path.exists(CONF_FILE):
        check_conf_file()
     users = {}
     if os.path.exists(CONF_FILE):
         conf=get_conf()
         users = conf["USERS"]
     return users

    
def get_driver(max_retries=5, retry_delay=5):
    if not os.path.exists(CONF_FILE):
        check_conf_file()
    if os.path.exists(CONF_FILE):
        conf=get_conf()
        DB_URL = conf["DB_URL"]
        AUTH= conf["AUTH"]
        for attempt in range(1, max_retries + 1):
            if test_connection(DB_URL, AUTH):
                driver = GraphDatabase.driver(DB_URL, auth=AUTH)
                return driver
            else:
                logger.warn(f"Retry {attempt}/{max_retries} to connect to driver")
                time.sleep(retry_delay)
                
    
        if test_connection(DB_URL, AUTH):
                driver = GraphDatabase.driver(DB_URL, auth=AUTH)
                return driver
        logger.error(f"❌ Fail to access Database of container {conf['container_name']}.")
        return None
    else:
        return None