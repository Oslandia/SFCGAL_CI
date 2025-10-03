#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXAMPLES_DIR="${PROJECT_DIR}/examples"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if module is built
if [ ! -f "${EXAMPLES_DIR}/sfcgal.wasm" ]; then
    echo -e "${RED}Error: SFCGAL module not built yet${NC}"
    echo -e "Run ${YELLOW}./scripts/build.sh${NC} first"
    exit 1
fi

# Find available port
PORT=8888
while lsof -Pi :$PORT -sTCP:LISTEN -t >/dev/null 2>&1; do
    PORT=$((PORT + 1))
done

echo -e "${GREEN}Starting web server on port ${PORT}...${NC}"
echo -e "Open ${YELLOW}http://localhost:${PORT}${NC} in your browser"
echo
echo "Press Ctrl+C to stop the server"

cd "${EXAMPLES_DIR}"
python3 -m http.server ${PORT} --bind 127.0.0.1