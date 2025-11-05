#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:16:02 2025

@author: fgraziani
"""

import os
import requests
import zipfile
import io
import shutil
from dash import html, Input, Output, callback, State, dcc
from Bio.Seq import Seq
from app import *
from neo4j_requests import *
from auth_utils import require_authorization
import logging


logger = logging.getLogger("panorama_logger")


# Dépôt cible (format : "owner/repo")
repo = "Arnifabiziagne/Panorama"
update_dir = "./update"
panorama_dir = "."

# URL de l'API pour la dernière release
url = f"https://api.github.com/repos/{repo}/releases/latest"

@require_authorization
def copy_update_to_root():
    for root, dirs, files in os.walk(update_dir):
        rel_path = os.path.relpath(root, update_dir)
        target_dir = os.path.join(panorama_dir, rel_path)

        os.makedirs(target_dir, exist_ok=True)

        for file in files:
            src_file = os.path.join(root, file)
            dest_file = os.path.join(target_dir, file)
            shutil.copy2(src_file, dest_file)

    if os.path.exists(update_dir):
        shutil.rmtree(update_dir)
    logger.info("✅ The panorama update has been successfully completed. It is required to restart the server.")


@app.callback(
    Output('update-panorama-output', 'children'),
    Input('update-panorama-btn', 'n_clicks'),
    prevent_initial_call=True
)
@require_authorization
def update_panorama(n_clicks):

    # Requête GET
    response = requests.get(url)
    
    # Vérification du code de statut
    if response.status_code == 200:
        release = response.json()
        version = release["tag_name"]
        major_version = version.split(".")[0]
        if major_version == DB_VERSION.split(".")[0]:

            logger.info("Latest release found :", release["name"])
            logger.info("Tag :", release["tag_name"])
            logger.info("Publication date :", release["published_at"])
            zip_url = release.get("zipball_url")
            logger.info(f"Download zip file : {zip_url}")
            
            r = requests.get(zip_url)
            if r.status_code != 200:
                logger.error(f"Error {r.status_code} while downloading")
                exit(1)
            
            # Créer dossier update s'il n'existe pas
            os.makedirs(update_dir, exist_ok=True)
            
            # Dézipper le contenu de l'archive dans /update en écrasant
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                # Attention : l'archive GitHub zipball a un dossier racine unique
                root_folder = z.namelist()[0].split('/')[0]
            
                for member in z.namelist():
                    # retirer le dossier racine unique
                    filename = member[len(root_folder)+1:]
                    if not filename:
                        continue  # passer le dossier racine
            
                    dest_path = os.path.join(update_dir, filename)
                    if member.endswith('/'):
                        os.makedirs(dest_path, exist_ok=True)
                    else:
                        with open(dest_path, 'wb') as f:
                            f.write(z.read(member))
            copy_update_to_root()
            return html.Div(f"✅ The panorama update has been successfully completed. It is required to restart the server.")
        else:
            logger.info(f"Latest version {version} is not compatible with the current data. To use latest release it is required to regenerate data.")
            return html.Div(f"❌ Latest version {version} is not compatible with the current data. To use latest release it is required to regenerate data.")
            
    else:
        logger.error("Error :", response.status_code)
        return html.Div(f"❌ Error : {response.status_code}")