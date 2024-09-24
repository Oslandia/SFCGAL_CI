#!/bin/bash

export LD_LIBRARY_PATH=/SFCGAL/build/src:${LD_LIBRARY_PATH}
export LIBRARY_PATH=/SFCGAL/build/src:${LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=/SFCGAL/build/include:${CPLUS_INCLUDE_PATH}
export C_INCLUDE_PATH=/SFCGAL/build/include:${C_INCLUDE_PATH}

exec "$@"
