#!/usr/bin/env python3
"""
Script to add container directives to all Snakemake rules
"""

import re
from pathlib import Path

# Environment mapping
ENV_MAPPING = {
    "../envs/ashleys_base.yaml": "ashleys_base",
    "../envs/ashleys_rtools.yaml": "rtools",
}


def update_rule_file(file_path):
    """Update a single rule file to include container directives"""
    print(f"Processing {file_path}...")

    with open(file_path, "r") as f:
        content = f.read()

    original_content = content

    # Pattern to match conda directives and add container directives
    for conda_env, env_name in ENV_MAPPING.items():
        # Look for conda directive with the environment path
        pattern = rf'(\s+)conda:\s*\n\s*["\']?{re.escape(conda_env)}["\']?'

        # Replacement includes both conda and container directives
        replacement = (
            rf"\1conda:\n"
            rf'\1    "{conda_env}"\n'
            rf"\1container:\n"
            rf'\1    CONTAINERS["{env_name}"]'
        )

        # Apply replacement
        content = re.sub(pattern, replacement, content, flags=re.MULTILINE)

    # Write back if content changed
    if content != original_content:
        with open(file_path, "w") as f:
            f.write(content)
        print(f"  ✓ Updated {file_path}")
        return True
    else:
        print(f"  - No changes needed for {file_path}")
        return False


def main():
    """Update all rule files"""
    rules_dir = Path("workflow/rules")

    # Find all .smk files except common.smk
    rule_files = [f for f in rules_dir.glob("*.smk") if f.name != "common.smk"]

    print(f"Found {len(rule_files)} rule files to update:")
    for rule_file in rule_files:
        print(f"  - {rule_file}")

    print("\\nUpdating rules...")

    updated_count = 0
    for rule_file in rule_files:
        if update_rule_file(rule_file):
            updated_count += 1

    print(f"\\nSummary: Updated {updated_count} out of {len(rule_files)} rule files")


if __name__ == "__main__":
    main()
