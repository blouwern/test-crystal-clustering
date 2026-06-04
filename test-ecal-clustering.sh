#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export ECAL_CLUSTERING_DIR=$SCRIPT_DIR
export ECAL_CLUSTERING_DATA_DIR=$SCRIPT_DIR/data
# export ECAL_CLUSTERING_SRC_DIR_CPP=$SCRIPT_DIR/algorithm/cpp
# export ECAL_CLUSTERING_SRC_DIR_PY=$SCRIPT_DIR/algorithm/python
