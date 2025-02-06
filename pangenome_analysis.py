#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:27:11 2025

@author: fgraziani
"""

import panorama
import os
import sys
import getopt


#Main program to compute the differents trees from a pangenome
if __name__ == "__main__":
    argv = sys.argv[1:]
    print("Command : " + str(sys.argv[0:]))
    # Récupérer les arguments passés depuis la ligne de commande
    if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print('Usage: python pantoolbox.py -f|--graphfilename <file_name:mandatory> -p|--project_name <project_name : mandatory> -n|--nodes_number  <nodes_number : optionnal, default=1000> -w|--weigthed <boolean : optionnal, default = True> -m|--method <methode : optionnal, default ="random"> -r|--redundancy <boolean : optionnal, default = True> -s|--strand <boolean : optionnal, default = True> -c|--colorfilename <colorfilename : optionnal, default = "">')
        sys.exit(1)
    else :
        file_name = ""
        project_name = ""
        nb_noeuds_cible = 10000
        pondere = True
        strand = True
        redudancy = True
        method = "random"
        color_file_name = ""
        
        try:
            opts, args = getopt.getopt(argv, "f:p:n:w:s:m:r:c:", ["graphfilename=", "projectname=", "nodes_number=","weighted=","strand=", "method=","redundancy=", "colorfilename="])
        except getopt.GetoptError:
            print("Bad argument")
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-f", "--graphfilename"):
                file_name = arg
            if opt in ("-p", "--projectname"):
                project_name = arg
            if opt in ("-n", "--nodes_number"):
                nb_noeuds_cible = int(arg)
            if opt in ("-w", "--weighted"):
                pondere = arg.lower() == "true"
            if opt in ("-s", "--strand"):
                strand = arg.lower() == "true"
            if opt in ("-m", "--method"):
                method = arg
            if opt in ("-c", "--colorfilename"):
                color_file_name = str(arg)
            if opt in ("-r", "--redundancy"):
                redudancy = arg.lower() == "true"

        if not isinstance(nb_noeuds_cible, int):
            print("nodes_number should be integer")
            exit(1)
        if file_name == "" or not os.path.exists(file_name):
            print("Unknown file")
            exit(1)
        if project_name == "":
            print("project name is mandatory")
            exit(1)
        if not isinstance(pondere, bool):
            print("Weighted argument must be a boolean True or False")
            exit(1)
        if not isinstance(strand, bool):
            print("Strand argument must be a boolean True or False")
            exit(1)
        if not isinstance(redudancy, bool):
            print("Redundancy argument must be a boolean True or False")
            exit(1)
        if color_file_name != "" and not os.path.exists(color_file_name):
            print("Unknown color_file_name")
            exit(1)
        panorama.analyser_pangenome(file_name, project_name, nb_noeuds_cible, pondere, method, redudancy, strand, color_file_name)