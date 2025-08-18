#!/bin/sh

DOXYGEN_CONFIG_FILE="Doxyfile_ci.check"
DOXYGEN_WARNINGS_LOG="doxygen_warnings.log"

rm -f ${DOXYGEN_CONFIG_FILE} ${DOXYGEN_WARNINGS_LOG}

cat <<EOF >> ${DOXYGEN_CONFIG_FILE}
QUIET = YES
WARNINGS = YES
WARN_IF_UNDOCUMENTED = YES
WARN_IF_DOC_ERROR = YES
WARN_NO_PARAMDOC = YES
RECURSIVE = YES
EXTRACT_ALL = YES
GENERATE_HTML = NO
GENERATE_LATEX = NO
GENERATE_XML = YES
WARN_LOGFILE = ${DOXYGEN_WARNINGS_LOG}
DOT_GRAPH_MAX_NODES = 100
PREDEFINED = DOXYGEN_SHOULD_SKIP_THIS
EOF

echo "changed files"
CHANGED_FILES=$(git diff-tree --diff-filter=d --name-only -r $CI_MERGE_REQUEST_DIFF_BASE_SHA $CI_COMMIT_SHA 'src/***.cpp' 'src/***.hpp' 'src/***.c' 'src/***.h')
echo $CHANGED_FILES;

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
