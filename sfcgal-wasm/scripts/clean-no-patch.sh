#!/bin/bash

echo "=== Cleaning No-Patch Build ==="
echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${YELLOW}This will remove:${NC}"
echo "  - build/sfcgal-wasm-no-patch/"
echo "  - examples/sfcgal.js and sfcgal.wasm (if from no-patch build)"
echo
read -p "Continue? (y/N) " -n 1 -r
echo

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled"
    exit 0
fi

# Remove no-patch SFCGAL build
if [ -d "${PROJECT_DIR}/build/sfcgal-wasm-no-patch" ]; then
    echo -e "${YELLOW}Removing build/sfcgal-wasm-no-patch...${NC}"
    rm -rf "${PROJECT_DIR}/build/sfcgal-wasm-no-patch"
    echo -e "${GREEN}✅ Removed${NC}"
fi

# Remove binding outputs (they could be from either method)
if [ -f "${PROJECT_DIR}/examples/sfcgal.js" ]; then
    echo -e "${YELLOW}Removing examples/sfcgal.js and sfcgal.wasm...${NC}"
    rm -f "${PROJECT_DIR}/examples/sfcgal.js"
    rm -f "${PROJECT_DIR}/examples/sfcgal.wasm"
    echo -e "${GREEN}✅ Removed${NC}"
fi

echo
echo -e "${GREEN}✅ Clean complete${NC}"
echo
echo "To rebuild with no-patch method:"
echo "  ./scripts/build-all-no-patch.sh"
echo
echo "Or step by step:"
echo "  ./scripts/build-deps-no-patch.sh     # If deps need rebuilding"
echo "  ./scripts/build-sfcgal-wasm-no-patch.sh"
echo "  ./scripts/build.sh"
