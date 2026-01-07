#!/bin/bash
# changelog_generator.sh - Generate changelog with PR references and commit hashes
# Adapted for MosaiCatcher Pipeline with assembly-specific container tags

# Get range (can be tag, commit, or branch names)
PREVIOUS_TAG=${1:-""}
CURRENT_REF=${2:-HEAD}

# If no previous tag provided, try to find the last release tag (excluding current)
if [ -z "$PREVIOUS_TAG" ] && [ "$CURRENT_REF" != "HEAD" ]; then
  PREVIOUS_TAG=$(git tag -l --sort=-version:refname | grep -v "^${CURRENT_REF}$" | head -1)
fi

# Final fallback to git describe if still no tag found
if [ -z "$PREVIOUS_TAG" ]; then
  PREVIOUS_TAG=$(git describe --tags --abbrev=0 HEAD^ 2>/dev/null || echo "")
fi

CURRENT_VERSION=${CURRENT_REF#v}

# Start with Container Images section
echo "### Container Images 🐳"
echo ""
echo "All images available on **GitHub Container Registry** and **Docker Hub**:"
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
echo ""
echo "# Also available on Docker Hub"
echo "weber8thomas/mosaicatcher-pipeline:<tag>"
echo "\`\`\`"
echo ""

# Get commits
if [ -z "$PREVIOUS_TAG" ]; then
  COMMITS=$(git log --pretty=format:"* %s [%h]" "$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
else
  COMMITS=$(git log --pretty=format:"* %s [%h]" "$PREVIOUS_TAG".."$CURRENT_REF" | grep -v "Merge" | sed '/^$/d')
fi

# Extract PR numbers and add links
CHANGELOG=$(echo "$COMMITS" | sed -E 's|\(#([0-9]+)\)|([#\1](https://github.com/'"${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}"'/pull/\1))|g')

echo -e "\n<details>\n<summary>Click to expand the changelog for $CURRENT_VERSION</summary>\n"

echo "### Changes 📜"
echo ""

# Features
FEATURES=$(echo "$CHANGELOG" | grep -E "^\\* (feat|feature|add):" || echo "")
if [ -n "$FEATURES" ]; then
  echo "#### New Features ✨"
  echo ""
  echo "$FEATURES"
  echo ""
fi

# Bug fixes
FIXES=$(echo "$CHANGELOG" | grep -E "^\\* (fix|bug|issue):" || echo "")
if [ -n "$FIXES" ]; then
  echo "#### Bug Fixes 🐛"
  echo ""
  echo "$FIXES"
  echo ""
fi

# Improvements
IMPROVEMENTS=$(echo "$CHANGELOG" | grep -E "^\\* (refactor|perf|style|improve|update|enhance):" || echo "")
if [ -n "$IMPROVEMENTS" ]; then
  echo "#### Improvements 🚀"
  echo ""
  echo "$IMPROVEMENTS"
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
  echo "$CHORES"
  echo ""
fi

# Documentation
DOCS_COMMITS=$(echo "$CHANGELOG" | grep -E "^\\* (docs):" || echo "")
if [ -n "$DOCS_COMMITS" ]; then
  echo "#### Documentation Updates 📚"
  echo ""
  echo "$DOCS_COMMITS"
  echo ""
fi

# Other changes
OTHER=$(echo "$CHANGELOG" | grep -v -E "^\\* (feat|feature|add|fix|bug|issue|refactor|perf|style|improve|update|enhance|chore|build|ci|docs):" || echo "")
if [ -n "$OTHER" ]; then
  echo "#### Other Changes 📝"
  echo ""
  echo "$OTHER"
  echo ""
fi

echo -e "\n</details>\n"

# Documentation section
echo "### Documentation 📖"
echo ""
echo "For more details, please refer to:"
echo "- [MosaiCatcher Documentation](https://friendsofstrandseq.github.io/mosaicatcher-docs/)"
echo "- [Version Management Guide](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/docs/version-management.md)"
echo "- [Release Workflow](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/docs/release-workflow.md)"
echo "- [Container Usage](https://github.com/${GITHUB_REPOSITORY:-friendsofstrandseq/mosaicatcher-pipeline}/blob/master/workflow/containers/README.md)"
