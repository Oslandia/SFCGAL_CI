#!/bin/sh

# Configuration
PROJECT_ID="19674165"
NEWS_FILE="NEWS"
TMPDIR=$(mktemp -d /tmp/sfcgal_news.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT INT HUP TERM

# Output files by category
FEATURES="$TMPDIR/features"
FIXES="$TMPDIR/fixes"
IMPROVES="$TMPDIR/improves"
CI="$TMPDIR/ci"
DOCS="$TMPDIR/docs"
TESTS="$TMPDIR/tests"
CHORE="$TMPDIR/chore"
MISC="$TMPDIR/misc"

# 1. Extract version from NEWS (first line like "2.1.0 (2025-05-14):")
VERSION=$(head -n 1 "$NEWS_FILE" | sed -n 's/^\([0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*\).*/\1/p')
TAG="v$VERSION"

if [ -z "$VERSION" ]; then
  echo "❌ Could not extract version number from $NEWS_FILE" >&2
  exit 1
fi

# 2. Get tag date
TAG_DATE=$(git for-each-ref --format="%(creatordate:iso8601-strict)" "refs/tags/$TAG")
if [ -z "$TAG_DATE" ]; then
  echo "❌ Git tag $TAG not found" >&2
  exit 1
fi

echo "📌 Tag: $TAG"
echo "📅 Fetching MRs merged since: $TAG_DATE"

# 3. Fetch merged MRs via GitLab API
MR_JSON=$(curl --silent \
  "https://gitlab.com/api/v4/projects/$PROJECT_ID/merge_requests?state=merged&updated_after=$TAG_DATE&per_page=100")

# 4. Parse and categorize MRs
echo "$MR_JSON" | jq -r '.[] | "- !\(.iid): \(.title) (\(.author.name))"' \
  | while IFS= read -r line; do
    case "$line" in
      *feat*|*Feat*|*feature*|*Feature*) echo "$line" >> "$FEATURES" ;;
      *fix*|*Fix*) echo "$line" >> "$FIXES" ;;
      *refactor*|*update*|*Update*|*improve*|*Improve*) echo "$line" >> "$IMPROVES" ;;
      *ci*|*CI*) echo "$line" >> "$CI" ;;
      *doc*|*Doc*|*capi:*) echo "$line" >> "$DOCS" ;;
      *test*|*Test*) echo "$line" >> "$TESTS" ;;
      *chore*|*build*|*release*|*sonar*|*clang-tidy*|*docker*) echo "$line" >> "$CHORE" ;;
      *) echo "$line" >> "$MISC" ;;
    esac
done

# 5. Print grouped output
print_section() {
  title="$1"
  file="$2"
  if [ -s "$file" ]; then
    echo "## $title"
    cat "$file"
    echo ""
  fi
}

print_section "New Features (Feat)" "$FEATURES"
print_section "Bug Fixes (Fix)" "$FIXES"
print_section "Improvements (Improve/Update)" "$IMPROVES"
print_section "Continuous Integration (CI)" "$CI"
print_section "Documentation (Docs)" "$DOCS"
print_section "Tests (Tests)" "$TESTS"
print_section "Chore / Maintenance (Chore)" "$CHORE"
print_section "Misc" "$MISC"
