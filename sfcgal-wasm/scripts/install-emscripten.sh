#!/bin/bash

echo "=== Installing Emscripten SDK ===" echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EMSDK_DIR="${PROJECT_DIR}/build/emsdk"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

if [ -d "${EMSDK_DIR}" ]; then
    echo -e "${YELLOW}Emscripten SDK already installed at ${EMSDK_DIR}${NC}"
    echo -e "${YELLOW}To reinstall, remove the directory first: rm -rf ${EMSDK_DIR}${NC}"
    exit 0
fi

echo -e "${YELLOW}Downloading Emscripten SDK...${NC}"
mkdir -p "${PROJECT_DIR}/build"
cd "${PROJECT_DIR}/build"

git clone https://github.com/emscripten-core/emsdk.git
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to clone emsdk repository${NC}"
    exit 1
fi

cd emsdk

echo -e "${YELLOW}Installing latest Emscripten...${NC}"
./emsdk install latest
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to install Emscripten${NC}"
    exit 1
fi

echo -e "${YELLOW}Activating Emscripten...${NC}"
./emsdk activate latest
if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to activate Emscripten${NC}"
    exit 1
fi

echo
echo -e "${GREEN}✅ Emscripten SDK installed successfully!${NC}"
echo
echo "To use Emscripten in your current shell, run:"
echo -e "${YELLOW}source ${EMSDK_DIR}/emsdk_env.sh${NC}"
echo
echo "Or add this to your ~/.bashrc or ~/.zshrc:"
echo -e "${YELLOW}source ${EMSDK_DIR}/emsdk_env.sh > /dev/null 2>&1${NC}"
