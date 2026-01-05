#!/usr/bin/env python3
"""
Update container references in Snakemake rule files based on current hashes.
This script reads the container hashes from containers.yaml and updates all
rule files to use the current container images.
"""

import hashlib
import re
from pathlib import Path

import yaml


def calculate_env_hash(env_file_path):
    """Calculate hash of environment file content."""
    with open(env_file_path, "rb") as f:
        content = f.read()
    return hashlib.sha256(content).hexdigest()[:12]


def get_current_containers():
    """Get current container references with updated hashes."""
    containers = {}

    # Calculate current hashes for each environment
    for env_name in ["ashleys_base", "ashleys_rtools"]:
        env_file = Path(f"workflow/envs/{env_name}.yaml")
        if env_file.exists():
            env_hash = calculate_env_hash(env_file)
            containers[
                env_name
            ] = f"ghcr.io/friendsofstrandseq/ashleys-qc-pipeline/envs/{env_name}:env-{env_hash}"

    return containers


def update_rule_file(file_path, containers):
    """Update container references in a single rule file."""
    with open(file_path, "r") as f:
        content = f.read()

    updated = False

    for env_name, container_ref in containers.items():
        # Pattern to match container directives
        pattern = rf'container:\s*"ghcr\.io/friendsofstrandseq/ashleys-qc-pipeline/envs/{env_name}:env-[a-f0-9]+"'
        replacement = f'container:\n        "{container_ref}"'

        new_content = re.sub(pattern, replacement, content, flags=re.MULTILINE)
        if new_content != content:
            content = new_content
            updated = True

    if updated:
        with open(file_path, "w") as f:
            f.write(content)
        print(f"Updated {file_path}")

    return updated


def main():
    # Get current container references
    containers = get_current_containers()

    print("Current container references:")
    for env_name, container_ref in containers.items():
        print(f"  {env_name}: {container_ref}")

    # Update all rule files
    rule_files = list(Path("workflow/rules").glob("*.smk"))
    total_updated = 0

    for rule_file in rule_files:
        if update_rule_file(rule_file, containers):
            total_updated += 1

    print(f"\nUpdated {total_updated} rule files")

    # Also update containers.yaml
    containers_yaml = Path("workflow/containers/containers.yaml")
    config = {"containers": containers}

    with open(containers_yaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Updated {containers_yaml}")


if __name__ == "__main__":
    main()
