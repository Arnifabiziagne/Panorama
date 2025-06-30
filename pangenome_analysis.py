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
        print("""Usage: python pantoolbox.py -f|--graphfilename <file_name:mandatory> 
              -d|--output_directory <output_directory : mandatory> 
              -p|--project_name <project_name : mandatory> 
              -n|--nodes_sampled_percent  <nodes_sampled_percent : optionnal, default=1> 
              -r|--redundancy <boolean : optionnal, default = True> 
              -s|--strand <boolean : optionnal, default = True> 
              -c|--colorfilename <colorfilename : optionnal, default = ""> 
              -k|--chromosome <chromosome : optionnal, default = ""> 
              -x|--maskednodesfilename <maskednodesfilename : optionnal, default = ""> 
              -a|--privatehaplotypesfilter <privatehaplotypesfilter : optionnal, default = "">'
              """)
        sys.exit(1)
    else :
        file_name = ""
        project_name = ""
        output_directory=""
        RaxML_sample_percent = 0.1
        pondere = True
        strand = True
        redudancy = True
        color_file_name = ""
        private_haplotypes_filter = 0
        masked_nodes_file_name = ""
        chromosome = None
        
        try:
            opts, args = getopt.getopt(argv, "f:d:p:n:s:r:c:k:x:a:", ["graphfilename=", "output_directory=", "projectname=", "nodes_sampled_percent=","strand=","redundancy=", "colorfilename=", "chromosome=", "maskednodesfilename=", "privatehaplotypesfilter="])
        except getopt.GetoptError as err:
            print("Bad argument " + str(err))
            sys.exit(2)
        for opt, arg in opts:
            if opt in ("-f", "--graphfilename"):
                file_name = arg
            if opt in ("-d", "--output_directory"):
                output_directory = arg
            if opt in ("-p", "--projectname"):
                project_name = arg
            if opt in ("-n", "--nodes_sampled_percent"):
                RaxML_sample_percent = float(arg)
            if opt in ("-s", "--strand"):
                strand = arg.lower() == "true"
            if opt in ("-c", "--colorfilename"):
                color_file_name = str(arg)
            if opt in ("-r", "--redundancy"):
                redudancy = arg.lower() == "true"
            if opt in ("-r", "--chromosome"):
                chromosome = arg
            if opt in ("-x", "--maskednodesfilename"):
                masked_nodes_file_name = arg
            if opt in ("-a", "--privatehaploypesfilter"):
                private_haplotypes_filter = arg

        if not isinstance(RaxML_sample_percent, float):
            print("nodes_sampled_percent should be integer")
            exit(1)
        if file_name == "" or not os.path.exists(file_name):
            print("Unknown file " + str(file_name))
            exit(1)
        if output_directory == "":
            print("output directory is mandatory")
            exit(1)
        if project_name == "":
            print("project name is mandatory")
            exit(1)
        if not isinstance(strand, bool):
            print("Strand argument must be a boolean True or False")
            exit(1)
        if not isinstance(redudancy, bool):
            print("Redundancy argument must be a boolean True or False")
            exit(1)
        if color_file_name != "" and not os.path.exists(color_file_name):
            print("Unknown color_file_name " + str(color_file_name))
            exit(1)
        if masked_nodes_file_name != "" and not os.path.exists(masked_nodes_file_name):
            print("Unknown masked_nodes_file_name " + str(masked_nodes_file_name))
            exit(1)

        panorama.analyser_pangenome(file_name, output_directory, project_name, RaxML_sample_percent, redudancy, strand, color_file_name, chromosome, masked_nodes_file_name, private_haplotypes_filter)