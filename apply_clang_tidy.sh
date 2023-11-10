#!/usr/bin/env bash

CLANG_TIDY_BIN=${CLANG_TIDY_BIN:-/usr/bin/clang-tidy}
RUN_CLANG_TIDY_BIN=${RUN_CLANG_TIDY_BIN:-/usr/bin/run-clang-tidy}
APPLY_CLANG_TIDY_BIN=${APPLY_CLANG_TIDY_BIN:-/usr/bin/clang-apply-replacements}

run_clang_tidy() {
    checks=$(clang-tidy15 --list-checks | awk '/  /{print $1}')

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
