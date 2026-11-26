#!/usr/bin/env bash
# Regression test for rule prepare_haplotag_input_bam (workflow/rules/haplotagging.smk).
#
# The rule must:
#   1. never change the mtime of the upstream bam/{cell}.sort.mdup.bam
#      (otherwise every rule that consumed it looks outdated and a resume re-runs
#      the whole downstream part);
#   2. never produce an output older than its inputs (snakemake's
#      check_output_mtime raises "Output ... has older modification time");
#   3. leave the upstream BAM intact when its temp() outputs are deleted.
#
# It builds a minimal workflow that mirrors the real chain
#   mark_duplicates -> index_input_bam -> check_sm_tag -> symlink_selected_bam
#   -> prepare_haplotag_input_bam (extracted verbatim from RULE_FILE) -> haplotag_bams
# and checks the points above, plus a resume and ITER forced re-runs of the rule.
#
# Usage: test_haplotag_input_mtime.sh RULE_FILE WORKDIR SEED_BAM [MODE] [ITER]
#   RULE_FILE  .smk file that contains rule prepare_haplotag_input_bam
#   WORKDIR    scratch directory (deleted and recreated)
#   SEED_BAM   small coordinate-sorted BAM carrying an @RG SM tag
#   MODE       symlink (default): selected/ is a symlink, as with use_light_data: False
#              copy:     selected/ is a copy, as with use_light_data: True
#              reheader: the BAM carries a pool SM tag, rewritten with fix_sm_tag: True
#   ITER       number of forced re-runs of the rule (default 300)
# Needs snakemake, samtools and python3 on PATH. Exit status = number of failed checks.
set -euo pipefail

RULE_FILE=$(readlink -f "$1")
WD=$2
SEED=$(readlink -f "$3")
MODE=${4:-symlink}
ITER=${5:-300}

FAILS=0
pass() { echo "PASS: $*"; }
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
mtime() { python3 -c 'import os, sys; print(os.stat(sys.argv[1]).st_mtime_ns)' "$1"; }
sm_of() { samtools view -H "$1" | grep -m1 '^@RG' | tr '\t' '\n' | sed -n 's/^SM://p'; }

SAMPLE=$(sm_of "$SEED")
CELL="${SAMPLE}TESTCELL"

rm -rf "$WD"
mkdir -p "$WD"
cd "$WD"
D="$WD/data"

FIX_SM=False
if [ "$MODE" = "reheader" ]; then
    # Pooled library: the BAM carries the pool name, not the sample.
    samtools reheader -c "sed 's/\tSM:[^\t]*/\tSM:TESTPOOL/'" "$SEED" > seed_pool.bam
    SEED="$WD/seed_pool.bam"
    FIX_SM=True
fi

# The rule declares conda: "../envs/mc_bioinfo_tools.yaml" relative to its own file;
# snakemake reads that file to record job metadata, so mirror the layout (conda
# itself is not used here: samtools comes from PATH).
mkdir -p rules envs
printf 'channels: [conda-forge, bioconda]\ndependencies: [samtools]\n' > envs/mc_bioinfo_tools.yaml
awk '/^rule prepare_haplotag_input_bam:/ {p = 1; print; next} p && /^rule / {exit} p' \
    "$RULE_FILE" > rules/prepare_rule.smk
grep -q 'rule prepare_haplotag_input_bam' rules/prepare_rule.smk || {
    echo "rule prepare_haplotag_input_bam not found in $RULE_FILE"
    exit 99
}

cat > Snakefile <<EOF
D = config["data_location"]
SEL = "cp" if config["mode"] == "copy" else "ln -s"


rule all:
    input:
        f"{D}/${SAMPLE}/haplotag/bam/${CELL}.bam.htg",
        f"{D}/${SAMPLE}/checks/${CELL}.sm_check.ok",


rule mark_duplicates:
    input:
        config["seed"],
    output:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    shell:
        "cp {input} {output}"


rule index_input_bam:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    shell:
        "samtools index {input}"


rule check_sm_tag:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/checks/{cell}.sm_check.ok",
    shell:
        "samtools view -H {input} | grep -q '^@RG' && echo ok > {output}"


rule symlink_selected_bam:
    input:
        bam="{folder}/{sample}/bam/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    output:
        bam="{folder}/{sample}/selected/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/selected/{cell}.sort.mdup.bam.bai",
    params:
        cmd=SEL,
    shell:
        "{params.cmd} {input.bam} {output.bam} && {params.cmd} {input.bai} {output.bai}"


include: "rules/prepare_rule.smk"


rule haplotag_bams:
    input:
        bam="{folder}/{sample}/haplotag/input/{cell}.bam",
        bai="{folder}/{sample}/haplotag/input/{cell}.bam.bai",
    output:
        "{folder}/{sample}/haplotag/bam/{cell}.bam.htg",
    shell:
        "samtools quickcheck {input.bam} && samtools view -c {input.bam} > {output}"
EOF

# Targets and flags go FIRST: --config (and -q in snakemake 9) take a variable
# number of values and would swallow any target placed after them.
# --max-checksum-file-size 0: snakemake ignores a newer mtime when the input is
# under the checksum threshold (1 MB by default) and its checksum is unchanged.
# Real BAMs are far above it, so disable checksums to behave like production.
smk() {
    snakemake "$@" --cores 1 --nolock --max-checksum-file-size 0 \
        --config "data_location=$D" "seed=$SEED" "mode=$MODE" "fix_sm_tag=$FIX_SM"
}
BAM="$D/$SAMPLE/bam/$CELL.sort.mdup.bam"
SEL_BAM="$D/$SAMPLE/selected/$CELL.sort.mdup.bam"
OK="$D/$SAMPLE/checks/$CELL.sm_check.ok"
HIN="$D/$SAMPLE/haplotag/input/$CELL.bam"
HTG="$D/$SAMPLE/haplotag/bam/$CELL.bam.htg"

plan() {
    # Rules a dry-run would execute, one per line (empty if nothing to do).
    # A failing dry-run prints "DRYRUN_FAILED" so it can never pass as "nothing to do".
    { smk -n 2>&1 || echo "DRYRUN_FAILED 0"; } |
        awk '/^DRYRUN_FAILED/ {print "DRYRUN_FAILED"; next} /^Job stats/ {p = 1; next} p && /^total/ {p = 0; next} p && NF == 2 && $2 ~ /^[0-9]+$/ {print $1}' |
        sort -u # snakemake 9 prints the job stats table twice
}

echo "== rule file: $RULE_FILE | mode: $MODE | workdir: $WD"

# Step 1: upstream part only (stands for the alignment/QC part of a run).
smk "$OK" "$SEL_BAM" > step1.log 2>&1 || { cat step1.log; exit 98; }
M0=$(mtime "$BAM")
sleep 1.1

# Step 2: full run, including prepare_haplotag_input_bam and the temp() cleanup.
if smk > step2.log 2>&1; then pass "full run succeeded"; else fail "full run failed (see $WD/step2.log)"; fi
M1=$(mtime "$BAM")
if [ "$M0" = "$M1" ]; then
    pass "upstream bam mtime unchanged by the full run"
else
    fail "upstream bam mtime changed by the full run ($M0 -> $M1)"
fi

# Step 3: temp() outputs removed, upstream BAM still there and intact.
if [ -e "$HIN" ] || [ -L "$HIN" ]; then
    fail "temp output $HIN still present after the run"
else
    pass "temp outputs removed"
fi
if [ -f "$BAM" ] && samtools quickcheck "$BAM"; then
    pass "upstream bam exists and passes quickcheck after temp() cleanup"
else
    fail "upstream bam missing or corrupt after temp() cleanup"
fi

# Step 4: nothing left to do.
P=$(plan)
if [ -z "$P" ]; then pass "dry-run after full run: nothing to be done"; else fail "dry-run after full run plans: $(echo $P)"; fi

# Step 5: resume after losing a downstream output.
rm -f "$HTG"
P=$(plan | tr '\n' ' ')
if [ "$P" = "all haplotag_bams prepare_haplotag_input_bam " ]; then
    pass "resume plans only: $P"
else
    fail "resume plans: $P (expected: all haplotag_bams prepare_haplotag_input_bam)"
fi
if smk > step5.log 2>&1; then pass "resume run succeeded"; else fail "resume run failed (see $WD/step5.log)"; fi
M2=$(mtime "$BAM")
if [ "$M0" = "$M2" ]; then pass "upstream bam mtime unchanged by the resume"; else fail "upstream bam mtime changed by the resume"; fi
P=$(plan)
if [ -z "$P" ]; then pass "dry-run after resume: nothing to be done"; else fail "dry-run after resume plans: $(echo $P)"; fi

# Step 6: ITER forced re-runs of the rule alone; snakemake's own check_output_mtime
# fails the run if an output ends up older than an input.
older=0
other=0
moved=0
for i in $(seq 1 "$ITER"); do
    if ! smk --force "$HIN" "$HIN.bai" > iter.log 2>&1; then
        if grep -q "older modification time" iter.log; then older=$((older + 1)); else other=$((other + 1)); fi
        cp iter.log "iter_fail_$i.log"
    fi
    [ "$(mtime "$BAM")" = "$M0" ] || moved=$((moved + 1))
done
if [ "$older" -eq 0 ] && [ "$other" -eq 0 ]; then
    pass "$ITER forced re-runs: 0 'older modification time', 0 other errors"
else
    fail "$ITER forced re-runs: $older 'older modification time', $other other errors"
fi
if [ "$moved" -eq 0 ]; then
    pass "$ITER forced re-runs: upstream bam mtime never changed"
else
    fail "$ITER forced re-runs: upstream bam mtime changed in $moved runs"
fi

# Step 7: the BAM handed to whatshap (explicit targets are kept despite temp()).
if samtools quickcheck "$HIN" && [ "$(samtools view -c "$HIN")" = "$(samtools view -c "$BAM")" ]; then
    pass "haplotag input bam is valid and has the same read count as upstream"
else
    fail "haplotag input bam invalid or read count differs"
fi
if [ "$(sm_of "$HIN")" = "$SAMPLE" ]; then
    pass "haplotag input bam has SM:$SAMPLE"
else
    fail "haplotag input bam has SM:$(sm_of "$HIN"), expected $SAMPLE"
fi
echo "   haplotag input: $(stat -c '%F' "$HIN"); rule log: $(tail -1 "$D"/log/prepare_haplotag_input_bam/*/*.log)"

echo "== $FAILS failed check(s)"
exit "$FAILS"
