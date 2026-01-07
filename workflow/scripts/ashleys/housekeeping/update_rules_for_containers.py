#!/usr/bin/env python3
"""
Script to update Snakemake rules to use containers or conda environments
"""

import re
from pathlib import Path

# Mapping of conda environment files to environment names
ENV_MAPPING = {
    "../envs/ashleys_base.yaml": "ashleys_base",
    "../envs/ashleys_rtools.yaml": "rtools",
}


def update_rule_file(file_path):
    """
    Update a Snakemake rule file to use get_container_or_conda function
    """
    print(f"Updating {file_path}...")

    with open(file_path, "r") as f:
        content = f.read()

    # Track if any changes were made
    changed = False

    # Replace conda environment declarations
    for conda_env, env_name in ENV_MAPPING.items():
        # Pattern to match conda environment declarations
        pattern = rf'(\s+)conda:\s*\n\s*["\']?{re.escape(conda_env)}["\']?'
        replacement = rf'\\1{env_name}_env = get_container_or_conda("{env_name}")\n\\1exec({env_name}_env)'

        if re.search(pattern, content):
            content = re.sub(pattern, replacement, content)
            changed = True

    # Write back if changed
    if changed:
        with open(file_path, "w") as f:
            f.write(content)
        print(f"  ✓ Updated {file_path}")
    else:
        print(f"  - No changes needed for {file_path}")

    return changed


def main():
    """
    Main function to update all rule files
    """
    # Find all Snakemake rule files
    rules_dir = Path("workflow/rules")
    rule_files = list(rules_dir.glob("*.smk"))

    print(f"Found {len(rule_files)} rule files to check:")
    for rule_file in rule_files:
        print(f"  - {rule_file}")

    print("\\nUpdating rule files...")

    total_changed = 0
    for rule_file in rule_files:
        if update_rule_file(rule_file):
            total_changed += 1

    print(f"\\nSummary: Updated {total_changed} out of {len(rule_files)} files")

    # Also need to update the configuration
    print("\\nUpdating configuration...")
    update_config()


def update_config():
    """
    Update the main configuration to add use_containers option
    """
    config_file = Path("config/config.yaml")

    if not config_file.exists():
        print("  ✗ config.yaml not found")
        return

    with open(config_file, "r") as f:
        content = f.read()

    # Add use_containers option if not present
    if "use_containers:" not in content:
        # Add after the version line
        content = content.replace(
            "version: 2.3.5",
            "version: 2.3.5\\n\\n# Container configuration\\nuse_containers: false",
        )

        with open(config_file, "w") as f:
            f.write(content)

        print("  ✓ Added use_containers option to config.yaml")
    else:
        print("  - use_containers option already present")


if __name__ == "__main__":
    main()
