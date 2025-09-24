Setup neo4j
---------
To setup neo4j, first download the neo4j_install directory. In the neo4j_install directory launch the setup_neo4j.sh script with the following parameters :
* --base-dir $neo4j_data_base_directory : by default the directory is set to ../data directory, if the data must be store in another directory, set it with this parameter.
* --dump $dump_file : optionnal. The defaut dump file is ../data/data/import/neo4j.dump. If the dump file is in another location, set it with this parameter.
* --image neo4j:$version : optionnal. Default is the 2025.05-community-bullseye version. The version of neo4j image must be compatible with the dump file.
* --http-port $HTTP_PORT : optionnal. Default is 7474. It is necessary to launch multiple database.
* --bolt-port $BOLT_PORT : optionnal. Default is 7687. It is necessary to launch multiple database.
* --container-name $CONTAINER_NAME : required : name of the docker container, must be unique.
* --max-mem : optionnal. Used when dump or import files are used to create DB. Defaut is 24G.
* --max-swap : optionnal. Used when dump or import files are used to create DB. Defaut is 25G, must be reater than max-mem.

Remarks :
---------
* The database can be construct from the GFA file instead of using a dump file, but il the dump file is available it is faster to load it than to reconstruct the whole database.
* The dump file (if required) is ideally located in ../import
* If data already exists, the script will check it and ask to delete existing data before dumping new data
* neo4j image must be compatible with the dump file
* APOC plugins for 2025-05-0-community-bullseye version is in the setup directory but if another version of neo4j is used, the apoc plugin will be downloaded
* Once database is up, you can access the data in a browser : http://localhost:$HTTP_PORT (replace $HTTP_PORT by the defined port, by default 7474). The default login is "neo4j" and password "Administrateur".

Basic usage exemple (default configuration)
---------
* the neo4j data are stored in "../data/data". 
* The dump file "neo4j.dump" is in "../import" directory (this file should be download before)
* The version of neo4j is 2025.05-community-bullseye
* The http port to use is 7474
* The bolt port to use is 7687

Command :
bash ./setup_neo4j.sh --container-name neo4j_pangenome_test


Usage exemple (non default configuration)
---------
* the neo4j data are stored in "~/work/project/neo4j". 
* The dump file "base.dump" is in "~/work/project/neo4j/impot" directory (this file should be download before)
* The version of neo4j is 2025.05-community-bullseye
* The http port to use is 7474
* The bolt port to use is 7687
* max memory to use : 32g
* max swap : 34g

Command :
bash ./setup_neo4j.sh --base-dir \~/work/project/neo4j --image neo4j:2025.05-community-bullseye --dump \~/work/project/neo4j/import/db.dump  --http-port 7474 --bolt-port 7687 --container-name neo4j_pangenome_test --max-mem 32g --max-swap 34g

Usage of neo4j container :
---------

Once the setup script terminated, the neo4j database could be managed with the following command (replace $HTTP_PORT by the value defined in the setup script, by default 7474) :
* Launching database : docker start $CONTAINER_NAME 
* Stopping database : docker stop $CONTAINER_NAME
* Restarting database : docker r$CONTAINER_NAME 
* Show logs : docker logs -f $CONTAINER_NAME
* Get the status : docker exec -it $CONTAINER_NAME neo4j status
* Delete the container (without deleting data) : docker rm -f $CONTAINER_NAME
* Dump the database (docker must be running) :
	- docker run --rm --user=neo4j -v $PATH_TO_DATA/data:/data -v $PATH_TO_DATA/import:/import neo4j:2025.05.0-community-bullseye  neo4j-admin database dump --to-path=/import neo4j
