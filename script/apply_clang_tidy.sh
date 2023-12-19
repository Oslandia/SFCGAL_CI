#!/usr/bin/env bash

# This script generates a commit for each clang-tidy check using various clang utilities.
#
# By default, it uses Linux paths.
# You can execute it by defining the variables CLANG_TIDY_BIN, RUN_CLANG_TIDY_BIN, 
# and APPLY_CLANG_TIDY_BIN when running the script:
#
# > env CLANG_TIDY_BIN=/usr/local/bin/clang-tidy15 \
#       RUN_CLANG_TIDY_BIN=/usr/local/bin/run-clang-tidy15 \
#       APPLY_CLANG_TIDY_BIN=/usr/local/bin/clang-apply-replacements15 \
#       ./apply_clang_tidy.sh
#
# There isn't a clang-tidy CI in this project as of now.
# This script is useful for potential changes in .clang-tidy.
# We encourage developers to apply the fixes suggested by clang-tidy.


CLANG_TIDY_BIN=${CLANG_TIDY_BIN:-/usr/bin/clang-tidy}
RUN_CLANG_TIDY_BIN=${RUN_CLANG_TIDY_BIN:-/usr/bin/run-clang-tidy}
APPLY_CLANG_TIDY_BIN=${APPLY_CLANG_TIDY_BIN:-/usr/bin/clang-apply-replacements}

run_clang_tidy() {
    checks=$(${CLANG_TIDY_BIN} --list-checks | awk '/  /{print $1}')

    # Iterate over each check and run clang-tidy with fixes
    while IFS= read -r check; do
        ${RUN_CLANG_TIDY_BIN} -p build -clang-tidy-binary=${CLANG_TIDY_BIN} -clang-apply-replacements-binary=${APPLY_CLANG_TIDY_BIN} -checks="-=*,$check" -fix

        # Check if changes were made
        if [ -n "$(git status -s)" ]; then
            # Commit the changes with a relevant message
            git add .
            git commit -m "Fix clang-tidy check: $check"
        else
            echo "No changes for clang-tidy check: $check"
        fi
    done <<< "$checks"
}

run_clang_tidy
