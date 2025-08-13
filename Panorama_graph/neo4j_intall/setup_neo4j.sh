#!/bin/bash
set -e

set -euo pipefail

# --- DEFAULT VALUES ---
DOCKER_IMAGE="neo4j:2025.05-community-bullseye"
HTTP_PORT="7474"
BOLT_PORT="7687"
CONF_SOURCE_FILE="./conf/neo4j.conf"
CONF_FILE="../db_conf.json"
NEO4J_BASE_DIR="../data"
DUMP_SOURCE_FILE_DEFAULT="../import/neo4j.dump"
DUMP_SOURCE_FILE=""

echo "NEO4J_BASE_DIR $NEO4J_BASE_DIR"

# --- HELP FUNCTION ---
function usage() {
  cat <<EOF
Usage: $0 --base-dir <dir> [--dump <dump_file>] [--image <neo4j_image>] [--http-port <port>] [--bolt-port <port>]

Options:
  --base-dir    Neo4j base directory (data, logs, conf, import, plugins) - required
  --dump        Neo4j dump file to import (optional)
  --image       Neo4j Docker image (default: $DOCKER_IMAGE)
  --http-port   HTTP port to expose (default: $HTTP_PORT)
  --bolt-port   Bolt port to expose (default: $BOLT_PORT)
  --container-name   Name of the docker container
  -h, --help    Show this help message

EOF
}

# --- OPTION PARSING ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir)
      NEO4J_BASE_DIR=$(realpath "$2")
      shift 2
      ;;
    --dump)
      DUMP_SOURCE_FILE=$(realpath "$2")
      shift 2
      ;;
    --image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --http-port)
      HTTP_PORT="$2"
      shift 2
      ;;
    --bolt-port)
      BOLT_PORT="$2"
      shift 2
      ;;
	--container-name)
      CONTAINER_NAME="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "❌ Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# --- VALIDATION ---
if [ -z "$CONTAINER_NAME" ]; then
  echo "❌ The --container-name parameter is required."
  usage
  exit 1
fi
NEO4J_AUTH="neo4j/Administrateur"
NEO4J_LOGIN="neo4j"
NEO4J_PASSWORD="Administrateur"
DUMP_DEST_DIR="../data/import"
DUMP_DEST_FILE="$DUMP_DEST_DIR/neo4j.dump"
PLUGINS_DIR="$NEO4J_BASE_DIR/plugins"
CONF_DIR="$NEO4J_BASE_DIR/conf"

# --- EXTRACT VERSION ---
NEO4J_VERSION=$(echo "$DOCKER_IMAGE" | cut -d':' -f2)
if [ -z "$NEO4J_VERSION" ]; then
  echo "⚠️ Unable to extract Neo4j version from Docker image; skipping APOC installation."
  APOC_VERSION=""
else
  APOC_VERSION="${NEO4J_VERSION}.0"
fi

APOC_JAR_NAME="apoc-${APOC_VERSION}-core.jar"
APOC_JAR_PATH="$PLUGINS_DIR/$APOC_JAR_NAME"
APOC_URL="https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/${APOC_JAR_NAME}"

# --- FUNCTION TO CHECK AND PREPARE data/databases/neo4j ---
function check_and_prepare_data_dir() {
  local data_db_dir="$NEO4J_BASE_DIR/data/databases/neo4j"
  if [ -d "$data_db_dir" ] && [ "$(ls -A "$data_db_dir")" ]; then
    echo "⚠️ Neo4j data already exists in $data_db_dir."
    read -p "Do you want to delete existing data to import the dump? (yes/no) " answer
    if [[ "$answer" =~ ^(yes|y|oui|o)$ ]]; then
      echo "🧹 Deleting old data..."
      rm -rf "$NEO4J_BASE_DIR/data"
      rm -rf "$NEO4J_BASE_DIR/conf"
	  rm -rf "$NEO4J_BASE_DIR/logs"
	  rm -rf "$NEO4J_BASE_DIR/plugin"
    else
      echo "❌ Import cancelled. Existing data kept."
      exit 1
    fi
  fi
}

echo "📁 Creating Neo4j directories (data, logs, conf, import, plugins)..."
for dir in data logs conf import plugins; do
  mkdir -p "$NEO4J_BASE_DIR/$dir"
done

if [ -d "$DUMP_DEST_DIR" ]; then
    # If DUMP_SOURCE_FILE is set (not empty) and the file exists
    if [ -n "$DUMP_SOURCE_FILE" ] && [ -f "$DUMP_SOURCE_FILE" ]; then
        echo "Moving $DUMP_SOURCE_FILE to $DUMP_DEST_FILE"
        mv "$DUMP_SOURCE_FILE" "$DUMP_DEST_FILE"
    # Otherwise, if the default dump file exists
    elif [ -f "$DUMP_SOURCE_FILE_DEFAULT" ]; then
        echo "Moving $DUMP_SOURCE_FILE_DEFAULT to $DUMP_DEST_FILE"
        mv "$DUMP_SOURCE_FILE_DEFAULT" "$DUMP_DEST_FILE"
    elif [ ! -f "$DUMP_DEST_FILE" ]; then
        echo "No dump file found. The database will start with no data."
    fi
else
    echo "Destination directory '$DUMP_DEST_FILE' does not exist."
    exit 
fi



if docker ps -a --format '{{.Names}}' | grep -q "^$CONTAINER_NAME$"; then
  echo "🧹 Removing existing container ($CONTAINER_NAME)..."
  docker rm -f $CONTAINER_NAME > /dev/null
fi


if [ -f "$CONF_SOURCE_FILE" ]; then
  echo " Copying conf file to conf directory🔧"
  cp "$CONF_SOURCE_FILE" "$CONF_DIR/"
else
  echo "⚠️ File $CONF_SOURCE_FILE doesn't exist"
fi



if [ -f "$DUMP_DEST_FILE" ]; then
  echo "📥 Preparing to import dump..."
  check_and_prepare_data_dir
  echo "🛑 Stopping Neo4j before import..."
  docker exec "$CONTAINER_NAME" neo4j stop || echo "ℹ️ Neo4j already stop"

  echo "📂 Importing dump..."
  docker run --rm \
  -v "$DUMP_DEST_DIR":/import \
  -v "$NEO4J_BASE_DIR/data":/data \
  "$DOCKER_IMAGE" \
  neo4j-admin database load neo4j --from-path=/import --overwrite-destination=true || {
    echo "❌ Fail loading dump"
    exit 1
  }
  
  echo "▶️ Restarting Neo4j..."
  docker start $CONTAINER_NAME
fi


# --- IMPORT CSV FILES IF PRESENT ---
CSV_NODES_FILE="$DUMP_DEST_DIR/nodes.csv"
CSV_RELATIONSHIPS_FILE="$DUMP_DEST_DIR/relations.csv"

if [ -f "$CSV_NODES_FILE" ] && [ -f "$CSV_RELATIONSHIPS_FILE" ]; then
  echo "📥 Detected CSV files for import (nodes.csv and relations.csv)"
  check_and_prepare_data_dir
  echo "🛑 Stopping Neo4j before CSV import..."
  docker exec "$CONTAINER_NAME" neo4j stop || echo "ℹ️ Neo4j already stopped"

  echo "📂 Importing CSV into new database..."
  docker run --rm \
    -v "$NEO4J_BASE_DIR/data":/data \
    -v "$DUMP_DEST_DIR":/import \
	-e NEO4J_AUTH=$NEO4J_AUTH \
    "$DOCKER_IMAGE" \
    neo4j-admin database import full \
      --verbose \
	  --nodes=Node=/import/nodes.csv \
	  --nodes=Sequence=/import/sequences.csv \
      --relationships=/import/relations.csv || {
        echo "❌ Failed to import CSV files"
        exit 1
      }
#  rm -f "$NEO4J_BASE_DIR/data/databases/neo4j/neostore.transaction.db*"
#  rm -rf "$NEO4J_BASE_DIR/data/transactions/neo4j"
#  echo "▶️ Restarting Neo4j..."
#  docker start $CONTAINER_NAME
fi

echo "🚀 Starting Neo4j ($DOCKER_IMAGE) on HTTP port $HTTP_PORT and Bolt port $BOLT_PORT..."
docker run -d \
  --name $CONTAINER_NAME \
  -e NEO4J_AUTH=$NEO4J_AUTH \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT=yes \
  -e NEO4J_apoc_export_file_enabled=true \
  -e NEO4J_apoc_import_file_enabled=true \
  -e NEO4J_apoc_import_file_use__neo4j__config=true \
  -e NEO4J_PLUGINS=\[\"apoc\"\] \
  -p "$HTTP_PORT:7474" \
  -p "$BOLT_PORT:7687" \
  -v "$NEO4J_BASE_DIR/data":/data \
  -v "$NEO4J_BASE_DIR/logs":/logs \
  -v "$NEO4J_BASE_DIR/conf":/conf \
  -v "$NEO4J_BASE_DIR/import":/import \
  -v "$NEO4J_BASE_DIR/plugins":/plugins \
  "$DOCKER_IMAGE"

echo "⏳ Waiting 10 seconds for Neo4j to start..."
sleep 10

# Creating JSON conf file
cat <<EOF > $CONF_FILE
{
  "container_name": "$CONTAINER_NAME",
  "http_port": $HTTP_PORT,
  "bolt_port": $BOLT_PORT,
  "login":"$NEO4J_LOGIN",
  "password":"$NEO4J_PASSWORD"
}
EOF


echo "✅ Neo4j is ready!"
echo "🌍 HTTP: http://localhost:$HTTP_PORT"
echo "🔗 BOLT: bolt://localhost:$BOLT_PORT"
echo "👤 Login: neo4j / password"