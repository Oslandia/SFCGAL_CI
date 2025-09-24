#!/bin/sh
#
# Doxygen documentation checker for specific files or directories
# POSIX-compatible script for FreeBSD, macOS, and Linux
# Usage:
#   ./doxygen-check-files.sh src/MyFile.cpp
#   ./doxygen-check-files.sh src/MyDirectory/
#   ./doxygen-check-files.sh src/File1.h src/File2.cpp
#   ./doxygen-check-files.sh src/  (check all files in src directory)
#

set -e

root_dir=$(dirname "$0")
# POSIX-compatible way to get absolute path
root_dir=$(cd "$root_dir/.." && pwd)
DOXYGEN_CONFIG_FILE="$root_dir/Doxyfile_ci.files"
DOXYGEN_WARNINGS_LOG="$root_dir/doxygen_warnings_files.log"

# Check if doxygen is available (POSIX way)
if ! command -v doxygen >/dev/null 2>&1; then
    echo "Error: doxygen is not available!" >&2
    echo "Please install doxygen to use this script." >&2
    exit 1
fi

# Function to display usage
show_usage() {
    echo "Usage: $0 <file_or_directory> [file_or_directory...]"
    echo ""
    echo "Examples:"
    echo "  $0 src/NURBSCurve.cpp"
    echo "  $0 src/NURBSCurve.h src/Curve.h"
    echo "  $0 src/algorithm/"
    echo "  $0 src/"
    echo ""
    echo "Supported file extensions: .cpp .hpp .c .h"
    exit 1
}

# Check if arguments provided
if [ $# -eq 0 ]; then
    echo "Error: No files or directories specified!" >&2
    show_usage
fi

# Clean up previous run
rm -f ${DOXYGEN_CONFIG_FILE} ${DOXYGEN_WARNINGS_LOG}

# Create Doxygen configuration
cat <<EOF > ${DOXYGEN_CONFIG_FILE}
QUIET = YES
WARNINGS = YES
WARN_IF_UNDOCUMENTED = YES
WARN_IF_DOC_ERROR = YES
WARN_NO_PARAMDOC = YES
RECURSIVE = YES
EXTRACT_PRIVATE = NO
GENERATE_HTML = NO
GENERATE_LATEX = NO
GENERATE_XML = YES
WARN_LOGFILE = ${DOXYGEN_WARNINGS_LOG}
DOT_GRAPH_MAX_NODES = 100
PREDEFINED = DOXYGEN_SHOULD_SKIP_THIS
ENABLE_PREPROCESSING = YES
EOF

# Function to check if file has valid extension (POSIX way)
is_valid_file() {
    file="$1"
    case "$file" in
        *.cpp|*.hpp|*.c|*.h) return 0 ;;
        *) return 1 ;;
    esac
}

# Function to get relative path (POSIX compatible)
get_relative_path() {
    target="$1"
    base="$2"

    # Simple relative path calculation
    case "$target" in
        "$base"/*)
            echo "${target#$base/}"
            ;;
        *)
            echo "$target"
            ;;
    esac
}

# Function to collect files from arguments
collect_files() {
    files_to_check=""
    temp_file="/tmp/doxygen_check_files.$$"

    for arg in "$@"; do
        # Convert to absolute path for consistency
        case "$arg" in
            /*)
                target="$arg"
                ;;
            *)
                target="$root_dir/$arg"
                ;;
        esac

        # Check if argument exists
        if [ ! -e "$target" ]; then
            echo "Warning: '$arg' does not exist, skipping..." >&2
            continue
        fi

        if [ -f "$target" ]; then
            # It's a file
            if is_valid_file "$target"; then
                # Convert back to relative path for doxygen
                rel_path=$(get_relative_path "$target" "$root_dir")
                files_to_check="$files_to_check $rel_path"
            else
                echo "Warning: '$arg' is not a supported file type, skipping..." >&2
            fi
        elif [ -d "$target" ]; then
            # It's a directory - find all valid files recursively
            find "$target" -type f \( -name "*.cpp" -o -name "*.hpp" -o -name "*.c" -o -name "*.h" \) > "$temp_file"
            while IFS= read -r file; do
                if is_valid_file "$file"; then
                    rel_path=$(get_relative_path "$file" "$root_dir")
                    files_to_check="$files_to_check $rel_path"
                fi
            done < "$temp_file"
            rm -f "$temp_file"
        else
            echo "Warning: '$arg' is neither a file nor directory, skipping..." >&2
        fi
    done

    echo "$files_to_check"
}

# Collect files to check
FILES_TO_CHECK=$(collect_files "$@")

# Remove leading/trailing whitespace
FILES_TO_CHECK=$(echo $FILES_TO_CHECK | xargs)

# Exit if no valid files found
if [ -z "$FILES_TO_CHECK" ]; then
    echo "Error: No valid files found to check!" >&2
    echo "Supported extensions: .cpp .hpp .c .h" >&2
    exit 1
fi

echo "Files to check:"
for file in $FILES_TO_CHECK; do
    echo "  $file"
done
echo ""

# Add files to doxygen config
echo "INPUT = $FILES_TO_CHECK" >> ${DOXYGEN_CONFIG_FILE}

# Run doxygen
echo "Running doxygen documentation check..."
if doxygen ${DOXYGEN_CONFIG_FILE} > /dev/null 2>&1; then
    # Check results
    if [ -s ${DOXYGEN_WARNINGS_LOG} ]; then
        echo "❌ Doxygen documentation errors detected:"
        echo "================================================"
        cat ${DOXYGEN_WARNINGS_LOG}
        echo "================================================"

        # Count warnings/errors
        warning_count=$(grep -c "warning:" ${DOXYGEN_WARNINGS_LOG} 2>/dev/null || echo 0)
        error_count=$(grep -c "error:" ${DOXYGEN_WARNINGS_LOG} 2>/dev/null || echo 0)

        echo ""
        echo "Summary: $warning_count warning(s), $error_count error(s) found."
        exit 1
    else
        echo "✅ No Doxygen documentation errors detected."
        echo "All specified files have proper documentation!"
        exit 0
    fi
else
    echo "❌ Doxygen execution failed!" >&2
    exit 1
fi