#!/usr/bin/env bash

root_dir=$(dirname $0)
root_dir=$(realpath $root_dir/..)
DOXYGEN_CONFIG_FILE="$root_dir/Doxyfile_ci.check"
DOXYGEN_WARNINGS_LOG="$root_dir/doxygen_warnings.log"

if ! type -p doxygen >/dev/null; then
    echo "doxygen is not available!" >> /dev/stderr
    exit
fi

rm -f ${DOXYGEN_CONFIG_FILE} ${DOXYGEN_WARNINGS_LOG}

cat <<EOF >> ${DOXYGEN_CONFIG_FILE}
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

if [ -n "$1" ]; then
    CHANGED_FILES="$*"
elif [ -n "$CI_MERGE_REQUEST_DIFF_BASE_SHA" ]; then
    CHANGED_FILES=$(git diff-tree --diff-filter=d --name-only -r $CI_MERGE_REQUEST_DIFF_BASE_SHA $CI_COMMIT_SHA)

else
    CHANGED_FILES=$(git diff --diff-filter=AM --name-only $(git merge-base HEAD master))
fi

if [ "$?" != 0 ]; then
    echo "Git operation failed!" >/dev/stderr
    exit 1
fi

# remove empty lines
CHANGED_FILES=$(echo $CHANGED_FILES | sed '/^[\s\n\r]*$/d')

# filter by extensions
NEW_FILES=""
for f in $CHANGED_FILES;
do
    if [[ "$f" =~ \.cpp$ ]] || [[ "$f" =~ \.hpp$ ]] || [[ "$f" =~ \.c$ ]] || [[ "$f" =~ \.h$ ]]; then
        NEW_FILES="$NEW_FILES $f"
    fi
done
CHANGED_FILES=$NEW_FILES

# exit if nothing to do
if [ -z "$CHANGED_FILES" ]; then
    echo "No valid changed files to document."
    exit 0
fi

echo "Files to check:"
echo $CHANGED_FILES

echo "INPUT = "${CHANGED_FILES} >> ${DOXYGEN_CONFIG_FILE}

doxygen ${DOXYGEN_CONFIG_FILE}

if [ -s ${DOXYGEN_WARNINGS_LOG} ]; then
    echo "Doxygen errors have been detected."
    cat ${DOXYGEN_WARNINGS_LOG}
    exit 1
else
    echo "No Doxygen errors detected."
    exit 0
fi
