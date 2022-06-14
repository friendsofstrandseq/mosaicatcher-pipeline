#!/bin/bash

OLD_TAG="$1" &&
    NEW_TAG=$(echo "$OLD_TAG" | sed 's/v//g') &&
    echo "$OLD_TAG" "$NEW_TAG" &&
    git tag "$NEW_TAG" "$OLD_TAG" &&
    git tag -d "$OLD_TAG" &&
    git push $2 "$NEW_TAG" :"$OLD_TAG"
