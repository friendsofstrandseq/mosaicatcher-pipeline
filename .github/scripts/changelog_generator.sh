#!/bin/bash
# changelog_generator.sh - Generate changelog with PR references and commit hashes
# Adapted for MosaiCatcher Pipeline with assembly-specific container tags

# Get range (can be tag, commit, or branch names)
PREVIOUS_TAG=${1:-""}
CURRENT_REF=${2:-HEAD}

# If no previous tag provided, find the last tag before current
if [ -z "$PREVIOUS_TAG" ]; then
  # Get all tags sorted by version, exclude current tag, get the first one (most recent)
  if [ "$CURRENT_REF" != "HEAD" ]; then
    # For a specific tag, find the previous tag
    PREVIOUS_TAG=$(git tag -l --sort=-version:refname | grep -v "^${CURRENT_REF}$" | grep -v "beta" | head -1)
  else
    # For HEAD, get the most recent tag
    PREVIOUS_TAG=$(git describe --tags --abbrev=0 HEAD^ 2>/dev/null || git tag -l --sort=-version:refname | head -1)
  fi
fi

CURRENT_VERSION=${CURRENT_REF#v}
PREVIOUS_VERSION=${PREVIOUS_TAG#v}
REPO=${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}
MAX_ITEMS=5

echo "Generating changelog from ${PREVIOUS_TAG} to ${CURRENT_REF}..." >&2

# Start with version comparison
echo "## What's Changed"
echo ""
echo "**Full Changelog**: https://github.com/$REPO/compare/${PREVIOUS_TAG}...${CURRENT_REF}"
echo ""

# Container Images section
echo "### Container Images 🐳"
echo ""
echo "All images available on **GitHub Container Registry**:"
echo ""
echo "\`\`\`plaintext"
echo "# Base image (all conda environments, no BSgenome)"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-$CURRENT_VERSION"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-latest"
echo ""
echo "# Assembly-specific images (base + BSgenome)"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-$CURRENT_VERSION"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg19-$CURRENT_VERSION"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-$CURRENT_VERSION"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-$CURRENT_VERSION"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm39-$CURRENT_VERSION"
echo ""
echo "# Latest tags (updated on stable releases)"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-latest"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg19-latest"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-latest"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-latest"
echo "ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm39-latest"
echo "\`\`\`"
echo ""

# Get commits (short hash for summary)
if [ -z "$PREVIOUS_TAG" ]; then
  COMMITS=$(git log --pretty=format:"* %s [%h]" "$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
  ALL_COMMITS=$(git log --pretty=format:"%H|%s" "$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
else
  COMMITS=$(git log --pretty=format:"* %s [%h]" "$PREVIOUS_TAG".."$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
  ALL_COMMITS=$(git log --pretty=format:"%H|%s" "$PREVIOUS_TAG".."$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
fi

# Extract PR numbers and add links
CHANGELOG=$(echo "$COMMITS" | sed -E 's|\(#([0-9]+)\)|([#\1](https://github.com/'"$REPO"'/pull/\1))|g')

# Helper function to show limited items from a category
show_limited_items() {
  local items="$1"
  local count=$(echo "$items" | wc -l)

  if [ "$count" -le "$MAX_ITEMS" ]; then
    echo "$items"
  else
    echo "$items" | head -n "$MAX_ITEMS"
    echo "* ... and $((count - MAX_ITEMS)) more"
  fi
}

echo "### Changes 📜"
echo ""

# Features
FEATURES=$(echo "$CHANGELOG" | grep -E "^\\* (feat|feature|add):" || echo "")
if [ -n "$FEATURES" ]; then
  echo "#### New Features ✨"
  echo ""
  show_limited_items "$FEATURES"
  echo ""
fi

# Bug fixes
FIXES=$(echo "$CHANGELOG" | grep -E "^\\* (fix|bug|issue):" || echo "")
if [ -n "$FIXES" ]; then
  echo "#### Bug Fixes 🐛"
  echo ""
  show_limited_items "$FIXES"
  echo ""
fi

# Improvements
IMPROVEMENTS=$(echo "$CHANGELOG" | grep -E "^\\* (refactor|perf|style|improve|update|enhance):" || echo "")
if [ -n "$IMPROVEMENTS" ]; then
  echo "#### Improvements 🚀"
  echo ""
  show_limited_items "$IMPROVEMENTS"
  echo ""
fi

# Breaking changes - look in commit bodies
if [ -n "$PREVIOUS_TAG" ]; then
  BREAKING=$(git log --pretty=format:"%b" "$PREVIOUS_TAG".."$CURRENT_REF" | grep -i "BREAKING CHANGE:" || echo "")
  if [ -n "$BREAKING" ]; then
    echo "#### Breaking Changes ⚠️"
    echo ""
    echo "$BREAKING"
    echo ""
  fi
fi

# Chores
CHORES=$(echo "$CHANGELOG" | grep -E "^\\* (chore|build|ci):" || echo "")
if [ -n "$CHORES" ]; then
  echo "#### Chores 🧹"
  echo ""
  show_limited_items "$CHORES"
  echo ""
fi

# Documentation
DOCS_COMMITS=$(echo "$CHANGELOG" | grep -E "^\\* (docs):" || echo "")
if [ -n "$DOCS_COMMITS" ]; then
  echo "#### Documentation Updates 📚"
  echo ""
  show_limited_items "$DOCS_COMMITS"
  echo ""
fi

# Other changes
OTHER=$(echo "$CHANGELOG" | grep -v -E "^\\* (feat|feature|add|fix|bug|issue|refactor|perf|style|improve|update|enhance|chore|build|ci|docs):" || echo "")
if [ -n "$OTHER" ]; then
  echo "#### Other Changes 📝"
  echo ""
  show_limited_items "$OTHER"
  echo ""
fi

# Full commit history with GitHub links
echo ""
echo "<details>"
echo "<summary>📝 Full Commit History (click to expand)</summary>"
echo ""
echo ""

while IFS='|' read -r hash message; do
  # Extract PR numbers and add links to the message
  message_with_pr=$(echo "$message" | sed -E 's|\(#([0-9]+)\)|([#\1](https://github.com/'"$REPO"'/pull/\1))|g')
  echo "* [\`${hash:0:7}\`](https://github.com/$REPO/commit/$hash) - $message_with_pr"
done <<< "$ALL_COMMITS"

echo ""
echo "</details>"
echo ""

# Documentation section
echo "### Documentation 📖"
echo ""
echo "For more details, please refer to:"
echo "- [MosaiCatcher Documentation](https://friendsofstrandseq.github.io/mosaicatcher-docs/)"
echo "- [Version Management Guide](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/docs/version-management.md)"
echo "- [Release Workflow](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/docs/release-workflow.md)"
echo "- [Container Usage](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/workflow/containers/README.md)"
