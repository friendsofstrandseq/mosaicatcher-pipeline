#!/usr/bin/env python3
"""
Simple script to check if container hashes match conda environment files.
Similar to a linting check - compares existing vs computed hashes.
"""

import hashlib
import sys
from pathlib import Path

import yaml


def calculate_env_hash(env_file_path):
    """Calculate hash of environment file content."""
    with open(env_file_path, "rb") as f:
        content = f.read()
    return hashlib.sha256(content).hexdigest()[:12]


def main():
    # Read container configuration
    containers_yaml = Path("workflow/containers/containers.yaml")
    with open(containers_yaml, "r") as f:
        config = yaml.safe_load(f)

    all_ok = True

    print("Checking container hashes...")
    print("-" * 60)

    for env_name, container_ref in config["containers"].items():
        # Calculate actual hash from conda environment file
        env_file = Path(f"workflow/envs/{env_name}.yaml")
        computed_hash = calculate_env_hash(env_file)

        # Extract hash from container reference
        if ":env-" in container_ref:
            existing_hash = container_ref.split(":env-")[-1]
        else:
            print(f"❌ {env_name}: No hash tag found in container reference")
            all_ok = False
            continue

        # Compare hashes
        if existing_hash == computed_hash:
            print(f"✅ {env_name}: Hash matches ({existing_hash})")
        else:
            print(f"❌ {env_name}: Hash mismatch!")
            print(f"   Existing: env-{existing_hash}")
            print(f"   Computed: env-{computed_hash}")
            all_ok = False

    print("-" * 60)

    if not all_ok:
        print("\n❌ Container hash check failed!")
        print(
            "To fix: python workflow/scripts/utils/update_container_references.py --use-hash"
        )
        sys.exit(1)
    else:
        print("\n✅ All container hashes are correct!")
        sys.exit(0)


if __name__ == "__main__":
    main()
