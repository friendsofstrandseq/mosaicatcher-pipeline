#!/usr/bin/env python3
"""
Convert hardcoded container references to function calls in Snakemake rules.
"""

import re
from pathlib import Path


def update_rule_file(file_path):
    """Update container references in a single rule file to use get_container function."""
    with open(file_path, "r") as f:
        content = f.read()

    # Pattern to match hardcoded container references
    pattern_base = r'container:\s*"ghcr\.io/friendsofstrandseq/ashleys-qc-pipeline/envs/ashleys_base:env-[a-f0-9]+"'
    pattern_rtools = r'container:\s*"ghcr\.io/friendsofstrandseq/ashleys-qc-pipeline/envs/ashleys_rtools:env-[a-f0-9]+"'

    # Replace with function calls
    content = re.sub(
        pattern_base,
        'container:\n            get_container("ashleys_base")',
        content,
        flags=re.MULTILINE,
    )
    content = re.sub(
        pattern_rtools,
        'container:\n            get_container("ashleys_rtools")',
        content,
        flags=re.MULTILINE,
    )

    # Write updated content
    with open(file_path, "w") as f:
        f.write(content)

    print(f"Updated {file_path}")


def main():
    # Update all rule files
    rule_files = list(Path("workflow/rules").glob("*.smk"))

    for rule_file in rule_files:
        if rule_file.name != "common.smk":  # Skip common.smk
            update_rule_file(rule_file)

    print(f"Updated {len(rule_files)-1} rule files")


if __name__ == "__main__":
    main()
