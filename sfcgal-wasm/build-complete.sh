#!/bin/bash

# Parse command-line arguments
CLEAN=false
for arg in "$@"; do
    case $arg in
        --clean)
            CLEAN=true
            shift
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: $0 [--clean]"
            echo "  --clean: Remove all build artifacts before building"
            exit 1
            ;;
    esac
done

cat << 'EOF'
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                       ║
║           SFCGAL WebAssembly - Complete Build Script                 ║
║                    (No-Patch Method)                                  ║
║                                                                       ║
║  This script will:                                                    ║
║  1. Install Emscripten SDK                                            ║
║  2. Download and build dependencies (GMP, MPFR, Boost, CGAL)         ║
║  3. Build SFCGAL as WebAssembly static library                        ║
║  4. Build the final WebAssembly binding                               ║
║                                                                       ║
║  Options:                                                             ║
║    --clean    Clean all build artifacts before building              ║
║                                                                       ║
║  Expected time: 20-30 minutes                                         ║
║  Disk space needed: ~2GB                                              ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
EOF

echo
echo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Configuration
LOG_FILE="${SCRIPT_DIR}/build-complete.log"
START_TIME=$(date +%s)

# Functions
log_step() {
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}▶ $1${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_error() {
    echo -e "${RED}❌ $1${NC}"
}

log_info() {
    echo -e "${YELLOW}ℹ $1${NC}"
}

show_elapsed_time() {
    local END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))
    local MINUTES=$((ELAPSED / 60))
    local SECONDS=$((ELAPSED % 60))
    echo -e "${CYAN}⏱  Elapsed time: ${MINUTES}m ${SECONDS}s${NC}"
}

check_prerequisites() {
    log_step "Checking prerequisites"

    local MISSING=()

    command -v cmake >/dev/null 2>&1 || MISSING+=("cmake")
    command -v git >/dev/null 2>&1 || MISSING+=("git")
    command -v curl >/dev/null 2>&1 || MISSING+=("curl")
    command -v tar >/dev/null 2>&1 || MISSING+=("tar")
    command -v make >/dev/null 2>&1 || MISSING+=("make")
    command -v python3 >/dev/null 2>&1 || MISSING+=("python3")

    if [ ${#MISSING[@]} -gt 0 ]; then
        log_error "Missing required tools: ${MISSING[*]}"
        echo
        echo "Please install them first. For example:"
        echo "  Ubuntu/Debian: sudo apt-get install cmake git curl tar make python3"
        echo "  Fedora/RHEL:   sudo dnf install cmake git curl tar make python3"
        echo "  macOS:         brew install cmake git curl"
        exit 1
    fi

    log_success "All prerequisites found"
    echo
}

# Trap errors
set -e
set -o pipefail
trap 'log_error "Build failed at step: $BASH_COMMAND"; exit 1' ERR

# Start
echo | tee "$LOG_FILE"
log_info "Build log: $LOG_FILE"
log_info "Working directory: $SCRIPT_DIR"
echo

# Clean if requested
if [ "$CLEAN" = true ]; then
    log_step "Cleaning build artifacts"

    echo -e "${YELLOW}This will remove:${NC}"
    [ -d "build/sfcgal-wasm" ] && echo "  - build/sfcgal-wasm/"
    [ -f "examples/sfcgal.js" ] && echo "  - examples/sfcgal.js"
    [ -f "examples/sfcgal.wasm" ] && echo "  - examples/sfcgal.wasm"
    [ -d "deps/install" ] && echo "  - deps/install/"
    [ -d "deps/gmp-6.3.0" ] && echo "  - deps/gmp-6.3.0/"
    [ -d "deps/mpfr-4.2.2" ] && echo "  - deps/mpfr-4.2.2/"
    [ -d "deps/CGAL-6.0.2" ] && echo "  - deps/CGAL-6.0.2/"
    echo -e "${YELLOW}Note: Emscripten SDK (build/emsdk/) will be preserved${NC}"
    echo -e "${YELLOW}Note: Boost headers (deps/boost_*/) will be preserved${NC}"
    echo
    read -p "Continue with cleanup? (y/N) " -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled by user"
        exit 0
    fi

    echo
    log_info "Removing build artifacts..."

    rm -rf build/sfcgal-wasm 2>/dev/null && log_success "Removed build/sfcgal-wasm/"
    rm -f examples/sfcgal.js 2>/dev/null && log_success "Removed examples/sfcgal.js"
    rm -f examples/sfcgal.wasm 2>/dev/null && log_success "Removed examples/sfcgal.wasm"
    rm -rf deps/install 2>/dev/null && log_success "Removed deps/install/"
    rm -rf deps/gmp-6.3.0 2>/dev/null && log_success "Removed deps/gmp-6.3.0/"
    rm -rf deps/mpfr-4.2.2 2>/dev/null && log_success "Removed deps/mpfr-4.2.2/"
    rm -rf deps/CGAL-6.0.2 2>/dev/null && log_success "Removed deps/CGAL-6.0.2/"

    echo
    log_success "Cleanup complete"
    echo
fi

# Ask for confirmation
echo -e "${YELLOW}This script will download and compile several large libraries.${NC}"
echo -e "${YELLOW}Continue? (y/N)${NC}"
read -r REPLY
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled by user"
    exit 0
fi
echo

# Step 1: Prerequisites
check_prerequisites

# Step 2: Install Emscripten
log_step "Step 1/4: Installing Emscripten SDK"

if [ -f "build/emsdk/emsdk" ]; then
    log_info "Emscripten already installed, skipping..."
else
    ./scripts/install-emscripten.sh 2>&1 | tee -a "$LOG_FILE"
fi

log_success "Emscripten SDK ready"
show_elapsed_time
echo

# Step 3: Source Emscripten environment
log_step "Activating Emscripten environment"

# Source may return non-zero even on success, so check differently
if ! source build/emsdk/emsdk_env.sh 2>&1 | tee -a "$LOG_FILE"; then
    if [ -z "$EMSDK" ]; then
        log_error "Failed to activate Emscripten environment"
        exit 1
    fi
fi

log_success "Emscripten environment activated"
echo

# Step 4: Build dependencies
log_step "Step 2/4: Building dependencies (GMP, MPFR, Boost, CGAL)"
log_info "This step takes ~15 minutes"
echo

if ! ./scripts/build-deps.sh 2>&1 | tee -a "$LOG_FILE"; then
    log_error "Dependencies build failed"
    exit 1
fi

log_success "Dependencies built successfully"
show_elapsed_time
echo

# Step 5: Build SFCGAL WebAssembly library
log_step "Step 3/4: Building SFCGAL WebAssembly library"
log_info "This step takes ~8 minutes"
echo

if ! ./scripts/build-sfcgal-wasm.sh 2>&1 | tee -a "$LOG_FILE"; then
    log_error "SFCGAL build failed"
    exit 1
fi

log_success "SFCGAL library built successfully"
show_elapsed_time
echo

# Step 6: Build WebAssembly binding
log_step "Step 4/4: Building WebAssembly binding"
log_info "This step takes ~20 seconds"
echo

if ! ./scripts/build.sh 2>&1 | tee -a "$LOG_FILE"; then
    log_error "WebAssembly binding build failed"
    exit 1
fi

log_success "WebAssembly binding built successfully"
show_elapsed_time
echo

# Final summary
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))
TOTAL_MINUTES=$((TOTAL_ELAPSED / 60))
TOTAL_SECONDS=$((TOTAL_ELAPSED % 60))

echo
log_step "Build Complete! 🎉"
echo

log_success "All components built successfully using method"
echo
log_info "Total build time: ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s"
echo

echo -e "${CYAN}Output files:${NC}"
ls -lh examples/sfcgal.js examples/sfcgal.wasm 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo
echo -e "${CYAN}Build method used:${NC}"
echo -e "  ${GREEN}✅ No-Patch${NC} - Clean build without CGAL modifications"
echo

echo -e "${CYAN}To test the examples:${NC}"
echo -e "  ${YELLOW}./scripts/serve.sh${NC}"
echo -e "  Then open ${BLUE}http://localhost:8888${NC}"
echo

echo -e "${CYAN}Build log saved to:${NC}"
echo -e "  ${YELLOW}$LOG_FILE${NC}"
echo

cat << 'EOF'
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                       ║
║                    ✅  BUILD SUCCESSFUL  ✅                           ║
║                                                                       ║
║  Your SFCGAL WebAssembly module is ready to use!                     ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
EOF
