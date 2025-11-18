# PANORAMA
**Tools to manipulate and visualize pangenomes variation graphs.**


## Overview
Panorama is based on a graph database modelisation (neo4j) from a gfa graph.
It allows the following functionalities : 
- Load a GFA file into a neo4j database 
- Load annotation files (GFF or GTF) : this will link annotations to the pangenome.
- Compute shared regions between a set of selected haplotypes
- Computes the sequences of a selected region
- Compute a global phylogenetic tree or a local phylogenetic tree from a selected region (neighbour joining with a distance matrix based on Jaccard index)
- Visualize a region and annotation of the pangenome

## Installation 
### Requirements
* Docker available : see docker documentation if not installed. Docker must be able to be launched by the $USER user; otherwise, see the procedure for launching Docker in non-root mode (adding docker group : "sudo groupadd docker" then "sudo usermod -aG docker $USER" and finally "newgrp docker")
* Miniconda 3. Please choose the installer corresponding to your OS: [Miniconda dowloads](https://docs.conda.io/en/latest/miniconda.html) 
* Mamba : this package will be automatically installed if not present.
* 20 GB RAM (32 Go+ recommanded for big data)
* Sufficient disk space : approximately 10 times the size of the GFA, ideally SSD (HDD are about 10 times slower and are not recommended for this use)
* The version of panorama used to construct the database must be compatible with the version to visualize and analyse data (the version of the database is indicated in the Stats node)

### Quickstart
* Installation (only the first time) :
  - Download the last release in the Panorama project and unzip the archive where do you want to store your data
  - Make launcher executable (linux) : chmod +x launcher.sh
* Launch the tool : 
  - Execute launcher : ./launch.sh on linux (or ./launch_gunicorn.sh to launch with gunicorn server in production) or ./launch.bat on windows and go to http://localhost:8050
  - To launch on another port (default is 8050), just specify the porty after. For example, to launch on port 8051 : ./launcher.sh 8051 or ./launcher.bat 8051.
* Prepare database (IHM)
  - Go to "DB management" page, give a name to your container (e.g. DB_my_species_PGGB) and click on "Create new DB"
  - If no data (dump or csv) present in the ./data/import directory, load GFA data by selecting GFA file. If the file concerns only a single chromosome it is required to set the chromosome name.
  - Load annotations by selecting the file (gtf or gff) and the genome associated to this file. Before to load annotations it is required that indexes are fully created : after creating data or loading GFA the indexes are automatically created but if data are big it requires some time.
  Once data are loaded the tool can be ued (see quick pages description).

## Important notes
  - To launch multiple neo4j instances, it is required to change neo4j ports. These ports are defined in the db_conf.json and can be updated here.
  - On Windows system, raxml-ng must be installed manually, see raxml documentation. If not installed, then the global phylogenetic tree could not be computed with this method (but the neighbor joining method will work).
  - The default memory used by the neo4j database is defined into the data/conf/neo4j.conf file, it requires at least 20 Go, if the system (and docker configuration) doesn't have this memory available it will be necessary to tune these values.
  - The GFA file must be properly structured for the application to correctly identify the individual name and chromosome. We strongly recommend to use W lines but according to the GFA format:
    - **For GFA files with `W` lines:**  
      The data is typically organized as follows:  
      - **Column 2:** Individual name  
      - **Column 3:** Haplotype number (the individual will then be named `individual_haplotype`)  
      - **Column 4:** Chromosome identifier  

    - **For `P` lines:**  
      The path name (`pathName`, in column 2) must follow one of the two formats below:  
      - `genome#haplotype#chromosome`  
      - `genome.haplotype.chromosome`

In all cases, if a chromosome is specified when loading the GFA file, that value will take precedence.

## Logs

Logs are displayed by default in the console and in log files located in the `./logs` directory.  

To configure logging behavior, modify the following parameters in the `./conf.json` file:

- **`log_retention_days`** — Defines the number of days to keep log files.  
- **`log_level`** — Set to `DEBUG`, `INFO`, `WARNING`, or `ERROR` to display only logs at or above the selected level.  
- **`log_server`** — Determines where logs are written:  
  - `"console"` — Log only to the console.  
  - `"file"` — Log only to a file.  
  - `"both"` — Log to both the console and file *(default value)*.  

## Quick pages description
  The menu allow to navigate on different pages :

  - DB management (available only in admin mode for server mode installation) : this page is used only on the start to create DB and load data.
  - Home page : page to visualize data (by defining the chromosome, start and stop or genome and gene_name / gene_id).
  - Shared region discovery : page to detect the nodes shared by a selection of haplotypes. It computes the list of identified regions that can be exported in csv. Sequence associated to a region can be seen by clicking on "size" column.
  - Phylogenetic : on left it is possible to load a reference phylogenetic tree. On right, by clicking on the "Plot tree..." button it computes the tree of the region defined in the Home page.
  - Sequences : by clicking on the button it computes the sequence for each haplotype of the region selected in the home page.


## Generate the database
* There are 3 ways to generate database :
    - From a dump file : this is the fastest way but the dump must be available. It uses the neo4j-admin load functionnality. If the neo4j.dump file is available, move or copy it into ./import directory, this file will be used to generate database.
    - From csv file : it is a fast way to create the database if the csv files are available. It uses the neo4j-admin import functionnality.
    - From the GFA file (and gtf / gff if available) : the database is constructed directly from the GFA file. This procedure can be usefull to add data but requires more time than other methods. It is not recommended for big gfa files.

## Server Mode Configuration

In server mode, if the server is publicly accessible, it may be necessary to restrict access and disable administration and file upload features.  

To do so, modify the configuration file (`./conf.json`):

- Set the `"server_mode"` parameter to `true`.  
- If temporary administration functions are needed (for example, for the initial data upload):  
  - Set the `"admin_mode"` parameter to `true`.  
  - In this case, the application will prompt for a username and password to access it.  
  - This login information is defined in the `"admin_users"` field — you should update it with the desired credentials.  
  - Once the data has been loaded, set `"admin_mode"` back to `false` to allow open access for all users.
- To launch application, use the ./launch_gunicorn.sh (server gunicorn is available only for linux) script.
- To stop application, use the ./stop_server.sh script.
- Gunicorn log into /logs/gunicorn but it is recommended to set a rotation file for this log. To do that :
  - Create file /etc/logrotate.d/gunicorn
- Write this into this file by replacing $ABSOLUTE_PATH_TO_PANORAMA by the aboslute path to your panorama directory :
  ``` bash
  $ABSOLUTE_PATH_TO_PANORAMA/logs/gunicorn/*.log {
      daily
      rotate 14
      missingok
      notifempty
      compress
      delaycompress
      copytruncate
  }
  ```


## Parameters file

The parameter file is named `./conf.json`. It contains the following parameters:

- `"container_name"`: set by the application, it is the name of the Docker container containing the Neo4j DB.  
- `"http_port"`: the HTTP port for the Neo4j DB (default: `7474`).  
- `"bolt_port"`: the Bolt port for the Neo4j DB (default: `7687`).  
- `"login"`: login to access the Neo4j DB (default: `"neo4j"`). If this value needs to be changed, it must first be modified in the Neo4j configuration file.  
- `"password"`: password to access the Neo4j DB (default: `"Administrateur"`). If this value needs to be changed, it must first be modified in Neo4j.  
- `"server_mode"`: set to `false` for local installation and `true` for server installation. This will block administration and file upload features.  
- `"admin_mode"`: only used in server mode. Set to `false` to block all administration features. If set to `true`, a login/password will be required to access the application.  
- `"admin_users"`: defines the login/password to access the application when `admin_mode` is set to `true`.  
- `"server_log_mode"`: `"console"`, `"file"`, or `"both"` — determines where logs are written.  
- `"log_retention_days"`: number of days to keep logs (default: `7`).  
- `"log_level"`: `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"` — logs are displayed only if their level is equal to or higher than this setting.  


## Contacts
*F Graziani, M Zytnicki*

Genotoul Bioinfo team, MIAT, INRAE Auzeville, France.

## Licence

Panorama is available under the GNU license.