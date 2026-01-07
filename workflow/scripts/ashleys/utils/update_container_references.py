#!/usr/bin/env python3
"""
Script to update container references in containers.yaml based on versioning scheme.
This script can be used to update containers to use different tags (version, hash, commit).
"""

import argparse
import hashlib
import os
import sys
from pathlib import Path

import yaml


def calculate_env_hash(env_file_path):
    """Calculate hash of environment file content."""
    with open(env_file_path, "rb") as f:
        content = f.read()
    return hashlib.sha256(content).hexdigest()[:12]


def update_container_references(config_file, version=None, commit=None, use_hash=False):
    """Update container references in containers.yaml."""

    # Read current configuration
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    # Get version from config.yaml if not provided
    if not version:
        config_yaml_path = (
            Path(__file__).parent.parent.parent.parent / "config" / "config.yaml"
        )
        with open(config_yaml_path, "r") as f:
            main_config = yaml.safe_load(f)
        version = f"v{main_config['version']}"

    # Update container references
    for env_name in config["containers"]:
        base_image = f"ghcr.io/friendsofstrandseq/ashleys-qc-pipeline/envs/{env_name}"

        if use_hash:
            # Calculate hash for the environment file
            env_file_path = (
                Path(__file__).parent.parent.parent / "envs" / f"{env_name}.yaml"
            )
            env_hash = calculate_env_hash(env_file_path)
            tag = f"env-{env_hash}"
        elif commit:
            tag = commit
        else:
            tag = version

        config["containers"][env_name] = f"{base_image}:{tag}"

    # Write updated configuration
    with open(config_file, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Updated container references in {config_file}")
    for env_name, image in config["containers"].items():
        print(f"  {env_name}: {image}")


def main():
    parser = argparse.ArgumentParser(description="Update container references")
    parser.add_argument(
        "--config",
        default="workflow/containers/containers.yaml",
        help="Path to containers.yaml file",
    )
    parser.add_argument("--version", help="Version tag to use (e.g., v2.3.5)")
    parser.add_argument("--commit", help="Commit hash to use")
    parser.add_argument(
        "--use-hash",
        action="store_true",
        help="Use environment file hash for container tags",
    )

    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: Config file {args.config} not found")
        sys.exit(1)

    update_container_references(args.config, args.version, args.commit, args.use_hash)


if __name__ == "__main__":
    main()
