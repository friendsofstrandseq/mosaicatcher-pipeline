# Python Club — Snakemake: What's new, how it compares to Nextflow, and HPC at EMBL

**Format:** 40-min seminar + Q&A
**Audience:** EMBL Python Club — mixed (some Snakemake v7 users, some Nextflow users, some never used a workflow manager)
**Tone:** High-level intro → low-level technical detail; strictly objective on the Snakemake vs Nextflow comparison
**Deliverable type:** Outline + speaker notes (this document)

---

## 0. Title options + organizer email

### Title candidates

1. **"Snakemake in 2026: what changed since v7, and how it stacks up against Nextflow on EMBL HPC"**
   *(neutral, descriptive — recommended for a Python Club)*

2. **"From Make to Plugin Soup: Snakemake v8/v9, a Nextflow side-by-side, and EMBL HPC field notes"**
   *(more playful; signals depth without overselling)*

3. **"Reproducible pipelines for Python people: Snakemake's last two major releases, an objective Nextflow comparison, and what we learned running on EMBL HPC"**
   *(longest, but maps cleanly to the three sections)*

### Draft email to the organizer (≈80 words)

> Subject: Python Club proposal — Snakemake v8/v9, Nextflow comparison, EMBL HPC tips
>
> Hi <name>,
>
> I'd like to propose a 40-min Python Club session on Snakemake. The talk has three parts: (1) a brief refresher on workflow managers and FAIR principles, (2) what actually changed in Snakemake between v7 and v9 (executor & storage plugins, software-deployment unification, new Python API), and (3) a strictly objective Snakemake vs Nextflow comparison on low-level technical points (dataflow model, resume/cache semantics, cloud, ecosystem). I'll close with a checklist of things to verify when running either tool on EMBL HPC, based on our recent work on MosaiCatcher. Audience: Python users with mixed pipeline experience. Brief (~5 min) live demo on SLURM included.
>
> Best,
> <you>

### Compact bullets if they want a 2-line abstract

- **What:** A neutral, technical walk-through of Snakemake's v8/v9 re-architecture, a like-for-like comparison with Nextflow DSL2/nf-core, and an EMBL-HPC checklist drawn from our MosaiCatcher experience.
- **Why now:** Snakemake 8 (Jan 2024) and 9 (Feb 2025) restructured the tool into plugins; many v7 workflows and CI configs at EMBL silently still target the old surface.

---

## 1. Timing budget (40 min, with 5 min Q&A)

| # | Section | Time | Cumulative |
|---|---|---|---|
| 1 | Hook + intro: WM, FAIR, Köster framing | 5 min | 5 |
| 2 | Snakemake core concepts (quick refresher) | 3 min | 8 |
| 3 | v7 → v8 → v9: what actually changed | 10 min | 18 |
| 4 | Snakemake vs Nextflow — low-level comparison | 9 min | 27 |
| 5 | Live demo on EMBL HPC (`dev/snakemake-v9-features/`) | 5 min | 32 |
| 6 | EMBL HPC checklist + gotchas | 3 min | 35 |
| 7 | Q&A | 5 min | 40 |

**Tradeoff if running long:** drop the live demo to a 30-second screen-grab and reclaim 4 min for Q&A — the demo material lives in `dev/snakemake-v9-features/` and can be circulated.

---

## 2. Section 1 — Hook + intro (5 min)

### 2.1 Hook (30 s)

Open with a familiar pain point, not a definition:

> "Raise your hand if your last pipeline was a chain of `bash` + `nohup` + a README with the commands in the right order. Keep it up if you've ever lost a week because someone re-ran step 3 with new inputs and forgot to re-run step 4."

### 2.2 What is a workflow manager? (1 min)

Three properties to anchor:

- **Declarative**: describe *what* should exist (the outputs), not *how* to get there step by step.
- **DAG-based**: the tool figures out the dependency graph and what is stale.
- **Portable**: same workflow definition runs on a laptop, SLURM, Kubernetes, AWS Batch — the executor is configuration, not code.

### 2.3 Why bother — the FAIR framing (1.5 min)

Map the FAIR principles (originally for *data*, but adopted for *workflows* by the WorkflowHub / EOSC / nf-core / Snakemake communities) onto what a workflow manager gives you:

| FAIR | Workflow-manager equivalent |
|---|---|
| **Findable** | Workflow Catalog entries; standard repo layout; CFF/Zenodo citation |
| **Accessible** | Git-hosted, container/conda specs bundled, can be re-run by anyone |
| **Interoperable** | Same workflow on laptop / SLURM / cloud; standard report formats |
| **Reproducible** | Pinned envs, content-addressed outputs, machine-readable provenance |

Köster's own framing in the rolling F1000 paper (Mölder et al. 2021, v3 Sep 2025) extends "reproducibility" into a **triad**: *reproducibility + adaptability + transparency*. Worth quoting verbatim, since the audience can immediately see why each leg matters:

> *"We postulate that it is equally important to ensure adaptability and transparency. The former describes the ability to modify the analysis to answer extended or slightly different research questions. The latter describes the ability to understand the analysis in order to judge whether it is not only technically, but methodologically valid."* — Mölder et al., F1000Research 10:33 (2021/2025)

### 2.4 A few numbers (30 s)

- Snakemake: founding paper Köster & Rahmann, *Bioinformatics* 2012; >1M downloads; >2000 citations; ~11 new citations/week (Köster, CERN 2024).
- Nextflow: founding paper Di Tommaso et al., *Nat. Biotechnol.* 2017; nf-core (*Genome Biology* July 2025): **124 curated pipelines, 1400+ modules, 70+ subworkflows, ~2600 contributors, ~10,800 Slack members**.

Frame: both tools are mature, both are mainstream in bioinformatics, this is not a "which is better" talk.

### 2.5 Further reading slide (cuts here in the deck)

Links to keep on a closing slide, not to read aloud:

- Mölder et al. 2021 (v3 2025): <https://doi.org/10.12688/f1000research.29032.3>
- Köster's slides hub: <https://slides.com/johanneskoester>
- ISMB/ECCB 2021 keynote: <https://www.youtube.com/watch?v=LchfsoH1ITg>
- FAIR Data Infrastructure Conf 2020 (most FAIR-focused): <https://www.youtube.com/watch?v=_dG9b3a9zkk>
- Snakemake Workflow Catalog: <https://snakemake.github.io/snakemake-workflow-catalog/>

---

## 3. Section 2 — Snakemake core concepts in 3 minutes (3 min)

The point of this section is to put everyone on the same page so the v8/v9 differences land. **Don't teach Snakemake here**; show *just enough* to read the rest of the slides.

### 3.1 The Snakefile object model

A `Snakefile` is just Python with extra top-level constructs. Three things on a slide:

```python
rule align:
    input:  "fastq/{sample}.fq.gz"     # wildcards
    output: "bam/{sample}.bam"
    threads: 8
    resources: mem_mb=16_000, runtime=120
    conda:  "envs/bwa.yaml"             # per-rule environment
    shell:  "bwa mem -t {threads} ref {input} | samtools sort -o {output}"
```

Three concepts visible: **rule**, **wildcard** (`{sample}`), **per-rule software spec** (`conda:` / `container:`).

### 3.2 Pull-based DAG (Make-style)

Snakemake builds the DAG **statically and up-front**. You ask for a target file; Snakemake recursively asks "to produce X, what input do I need, and which rule produces it?" until it hits existing files. `--dry-run` prints the whole plan before anything runs.

This is the single most important concept to keep visible — it shows up again in Section 4 as the deepest semantic difference from Nextflow.

### 3.3 Reproducibility primitives in one slide

- `conda:` / `container:` per rule
- `--report report.html` for a single-file shareable provenance report
- `.snakemake/` metadata tracks code, params, input set, env hashes for rerun triggers

---

## 4. Section 3 — Snakemake v7 → v8 → v9: what actually changed (10 min)

The headline to put on the title slide of this section:

> **v8 (Jan 2024) was a massive re-architecture; v9 (Feb 2025) is a near-invisible cleanup release. The disruption was the v7→v8 jump.**

Then walk through these axes — one slide each.

### 4.1 The executor plugin split (2 min)

**v7:** cluster integration was baked in. `--cluster "sbatch ..."`, `--drmaa`, `--kubernetes`, `--google-lifesciences`, `--tibanna`, `--tes`, `--cluster-generic-submit-cmd` — all shipped in core.

**v8:** every execution backend is now a separately pip-installable plugin selected by `--executor <name>`. Plugin authors subclass `RemoteExecutor` from `snakemake-interface-executor-plugins`; settings are declared as a `@dataclass` and become CLI flags automatically.

```bash
pip install snakemake-executor-plugin-slurm
snakemake --executor slurm --jobs 100 ...
```

Plugins as of late 2025: `slurm`, `lsf`, `htcondor`, `kubernetes`, `googlebatch` (replaces deprecated Google Life Sciences), `azure-batch`, `flux`, `cluster-generic` (back-compat wrapper). Catalog: <https://snakemake.github.io/snakemake-plugin-catalog/>.

**Worth pointing out:** the SLURM plugin submits *Snakemake itself* as the job script on the compute node — necessary so `run:` blocks and wrappers work remotely. Installing the executor plugin also installs a companion `jobstep` plugin used on the worker.

### 4.2 The storage plugin split (1.5 min)

**v7:** import per-protocol remote providers inside the Snakefile:

```python
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()
rule x:
    input: S3.remote("mybucket/in.txt")
```

**v8:** all remote IO moved out of core into **storage plugins** (`-s3`, `-gcs`, `-azure`, `-http`, `-ftp`, `-zenodo`, `-fs`, …). NCBI and EGA became Snakemake wrappers instead. The new directive:

```python
storage:
    provider="s3",
    max_requests_per_second=10,

rule example:
    input: storage("s3://mybucket/example.txt")
    output: "example.txt"
```

For project-wide remote IO: `--default-storage-provider s3 --default-storage-prefix s3://mybucket/`.

**The HPC-relevant one is `snakemake-storage-plugin-fs`** — used to stage data between shared FS and local scratch on compute nodes (see Section 6).

### 4.3 Software deployment unified (1 min)

**v7:** `--use-conda`, `--use-singularity`, `--use-envmodules`, `--conda-frontend mamba`.

**v8/v9:** one flag, `--software-deployment-method` (alias `--sdm`):

```bash
snakemake --sdm conda
snakemake --sdm apptainer        # was --use-singularity (singularity → apptainer rename)
snakemake --sdm conda apptainer  # both: container is the OS, conda installs tools inside
snakemake --sdm env-modules
```

`--conda-frontend mamba` is **gone** because `libmamba` is the conda default solver. Future direction (per the docs roadmap): software deployment itself will become a plugin, with a likely Pixi backend.

### 4.4 CLI rename cheat sheet (1 min)

Make this a single table slide and leave it on screen for a beat — people will photograph it:

| v7 | v8/v9 |
|---|---|
| `--cluster "sbatch ..."` | `--executor cluster-generic --cluster-generic-submit-cmd "sbatch ..."` |
| `--drmaa` | removed (use a specific executor plugin) |
| `--kubernetes` | `--executor kubernetes` |
| `--google-lifesciences` | `--executor googlebatch` |
| `--tibanna` | `--executor tibanna` |
| `--use-conda` | `--sdm conda` |
| `--use-singularity` | `--sdm apptainer` |
| `--singularity-args` | `--apptainer-args` |
| `--conda-frontend mamba` | removed (libmamba default) |
| `--restart-times N` | `--retries N` |

**Profiles are now versioned.** A profile directory may contain `config.yaml` (legacy) plus `config.v8+.yaml`; Snakemake picks the highest-version file the running engine understands. Existing v7 profiles **do not** work as-is on v8.

### 4.5 Reports got nicer (30 s)

- `category` / `subcategory` / `labels` can be functions of wildcards and (optionally) `input` / `output` / `params` — per-sample dynamic labelling.
- Toggle switches appear automatically when items differ on one labelled dimension.
- First-class Datavzrd embedding via `report(directory(...), htmlindex="index.html")`.
- `--report-stylesheet custom.css` for branding.
- `--report-after-run` (added in 9.x) auto-generates the report on successful completion.

### 4.6 Python API rewritten (1 min — only if audience cares about programmatic use)

**v7:** `from snakemake import snakemake; snakemake("Snakefile", ...)` with ~80 kwargs.

**v8/v9:** three-layer dataclass API: `SnakemakeApi → WorkflowApi → DAGApi`, mirroring load → resolve → execute. Verbose for simple cases (see issue [#2792](https://github.com/snakemake/snakemake/issues/2792)) but type-safe and IDE-discoverable. Worth one slide if anyone in the room embeds Snakemake in a larger tool.

### 4.6b Snakemake 2026 — Munich hackathon + v9.17→v9.21 (1.5 min)

The Snakemake community ran a hackathon at **MDSI / TU Munich, 9–13 March 2026** (40+ participants, fully booked). The BioHackrXiv post-event report ([10.37044/osf.io/h6zqj_v1](https://index.biohackrxiv.org/2026/04/27/h6zqj.html)) calls out three roadmap themes: **core performance on heterogeneous HPC, plugin extensions for domain-specific needs, and lowering the entry barrier for novices.** A burst of releases v9.17 → v9.21 landed tightly around and after the hackathon week.

Concrete additions (verified from the releases.atom feed):

| Version | Date | What it adds |
|---|---|---|
| 9.17.0 | 2026-03-13 | **Pluggable metadata persistence** (files or DB); **lambda rule `priority`**; topologically ordered job table; `on...` directive on modules |
| 9.18.0 | 2026-03-25 | **Storage plugins expose checksums** → cross-workflow caching correctness |
| 9.19.0 | 2026-03-28 | `profile.yaml` is the **new default filename**; you can pass a YAML file directly instead of a profile directory |
| 9.20.0 | 2026-05-02 | **Named multi-output caching**; `subpath(..., with_suffix=...)`; ⚠ **`runtime` CLI/profile values now interpreted as minutes (was seconds)** — migration footgun; plugins without settings auto-deploy |
| 9.21.0 | 2026-05-14 | Filename-prefix helper; tuned SQLite PRAGMAs; **auto-detected network filesystem type** |

The plugin architecture now spans **five orthogonal categories**: executor, storage, logger, **report**, and **scheduler**. Notable plugins this cycle: `snakemake-executor-plugin-cannon` (Harvard auto-partition selection), `snakemake-report-plugin-rocrate` (Workflow-Run RO-Crate provenance, Leo et al. 2024), `snakemake-report-plugin-nanopub`, community MILP scheduler plugin proposed via issue [#3646](https://github.com/snakemake/snakemake/issues/3646).

> **Honesty caveat for the talk:** The slides.com deck `snakemake-intro-updates-2026` is unlisted/403; the items above are reconstructed from the official Snakemake releases feed and the BioHackrXiv hackathon report, not transcribed from Köster's slides. Be explicit about that if asked.

### 4.7 So what's actually new in v9 specifically? (30 s)

Per the official migration page:

> *"Between Snakemake 8 and Snakemake 9, there is only a single breaking change in how custom loggers are provided, such that hardly any user should be affected."*

That change: **logger plugins** (`snakemake-logger-plugin-<name>`), interface in `snakemake-interface-logger-plugins`. Motivation: the v7/v8 synchronous log handler could crash entire workflows on monitoring-backend hiccups; the plugin model decouples I/O. Existing plugin: `snakemake-logger-plugin-snkmt` (SQLite/observability backend).

Cumulative point releases across 9.x also added: conda pinnings/post-deploy scripts/env-modules included in rerun-trigger detection, version-agnostic pickling of pandas/polars/numpy `params`, robust topological ordering when touching group outputs.

### 4.8 Things that bite you in v7 → v8 migration (1.5 min)

A "watch out" slide — these all surfaced in MosaiCatcher CI:

1. `dynamic()` outputs removed → use **checkpoints**.
2. `version:` directive removed → use `conda:` / `container:`.
3. `subworkflow:` directive removed → use `module:` + `use rule * from foo as foo_*`.
4. `from snakemake.remote.* import ...` — gone; rewrite to `storage`.
5. Profiles must be updated; old `config.yaml` references dead flags.
6. CI must install the executor plugin explicitly — `pip install snakemake-executor-plugin-slurm` is no longer implicit.
7. NCBI / EGA remote providers became **wrappers**, not storage plugins (different API entirely).
8. **9.20 (May 2026)**: `runtime` values in profiles/CLI now interpreted as **minutes**, not seconds. Audit before upgrading — `runtime: 60` silently changes from 1 min to 1 h.

---

## 5. Section 4 — Snakemake vs Nextflow, low-level comparison (9 min)

Frame this section on the title slide:

> **No "which is better." Both are mature, both are in production at EMBL. We are comparing them on technical primitives so you can pick the right one for your problem.**

Lead with this honest framing slide — it sets the tone for the rest.

### 5.1 The deepest semantic difference: dataflow model (2 min)

This is the slide that explains every other downstream difference. Spend time here.

- **Snakemake — pull-based, target-driven, Make-style.** You declare rules whose `output:` are file paths (often with wildcards). Execution starts from a *target* (first rule, `--until X`, explicit file path) and Snakemake recursively asks "to produce this, what do I need?" until the chain bottoms out at existing inputs.
- **Nextflow — push-based dataflow, Kahn/CSP-style** (implemented on the GPars library). You declare `process` definitions and connect them with `channel`s. Values flow asynchronously from emitters into consumers; a process fires as soon as one input record is available on every input channel. **No global target** — execution proceeds by token availability.

Pragmatic phrasing for the slide:
- *Snakemake:* "what file do I want?"
- *Nextflow:* "what stream of records flows into what process?"

DAG construction is a direct consequence:
- Snakemake builds the DAG **statically up front** (`--dry-run` shows the full plan; `--dag`/`--rulegraph` emit Graphviz). Conditional flows are expressed via Python returning input lists from config; **checkpoints** handle the case where an output's *contents* determine subsequent inputs.
- Nextflow builds the DAG **dynamically at runtime**: scatter/gather emerges from operators (`splitFastq`, `groupTuple`); `-preview` and `-with-dag` give an approximation, the true graph only fully exists once channels emit.

### 5.2 Language and tooling (45 s)

| | Snakemake | Nextflow |
|---|---|---|
| Host language | Python-embedded DSL | Groovy DSL2 on the JVM |
| Runtime | CPython | **Java 17+** (Nextflow 25.04, May 2025, dropped Java 11; supports up to Java 26) |
| Debugging | pdb + `--debug` / `--printshellcmds` | Log-driven via per-task `work/.command.sh`, `.command.log` |
| IDE | Python ecosystem + `snakefmt` | Seqera-maintained VS Code Nextflow Language Server |

> **Myth-buster:** "Nextflow needs Java" — *true*, and as of May 2025 it specifically needs Java 17+.

### 5.3 Resume / cache semantics (1.5 min)

This is the one most prone to common myths — handle carefully.

**Snakemake:** file **timestamps + tracked metadata** in `.snakemake/`. A job reruns if (a) the output is missing, (b) any input is newer than the output, or (c) tracked metadata flags changes: `--rerun-triggers {mtime,params,input,software-env,code}`. The classic Make pitfall (touching an input invalidates downstream) applies.

**Nextflow:** content-based **hash** stored in `.nextflow/cache/<session-id>/` (LevelDB). On `-resume`, each task's hash is recomputed from: the process script body, input values, **input file attributes (path, size, last-modified) by default — not file contents**, unless you opt in with `cache 'deep'`. The resolved container/conda directive participates. Debug with `-dump-hashes`.

> **Myth-buster:** "Nextflow `-resume` is content-hashed and therefore immune to timestamp issues." — *partially false*. By default it hashes file *attributes* (including mtime), not content. You must opt into `cache 'deep'` for true content hashing, and shared-FS timestamp jitter is a documented cause of spurious cache misses. Mitigation: `cache 'lenient'`.

Net: **both tools are sensitive to filesystem metadata by default**, both have escape hatches. Neither is magic.

### 5.4 Executor / scheduler integration (1 min)

| Backend | Snakemake (v9, plugin) | Nextflow (built-in module) |
|---|---|---|
| SLURM | ✅ `snakemake-executor-plugin-slurm` | ✅ `slurm` |
| LSF | ✅ `-lsf` | ✅ `lsf` |
| SGE / PBS / OAR / Moab / NQSII | community plugin status | ✅ all built-in |
| HTCondor | ✅ `-htcondor` | ✅ `condor` |
| Kubernetes | ✅ `-kubernetes` | ✅ `k8s` |
| Google Batch | ✅ `-googlebatch` | ✅ `google-batch` |
| Azure Batch | ✅ `-azure-batch` | ✅ `azurebatch` |
| AWS Batch | ❌ no official plugin (use Kubernetes-on-EKS, TES, or community `pcluster-slurm`) | ✅ `awsbatch` (first-class) |
| Flux | ✅ `-flux` | ✅ `flux` |
| TES | community plugin | ✅ `tes` |

> **Myth-buster:** "Snakemake can't do AWS Batch." — there's no official AWS-Batch executor plugin, but Kubernetes-on-EKS, TES, and community-maintained options exist. It is genuinely a thinner offering than Nextflow's `awsbatch` executor.

### 5.5 Software deployment & containers (45 s)

Largely parity in capability, different conventions:

- **Snakemake:** per-rule `conda:` (YAML) and `container:` (Docker/Apptainer URI). `--sdm conda apptainer` composes them (container = OS, conda = tools inside). No registry-mirror.
- **Nextflow:** per-process `conda` / `container` directives + engines (`docker.enabled`, `singularity.enabled`, …). nf-core conventions: Bioconda `environment.yml` + matching BioContainer (single-tool) or **mulled** BioContainer (multi-tool) from `quay.io/biocontainers`. nf-core has been migrating to **Seqera Containers**, built on the open-source **Wave** service, which builds containers on demand from conda specs and supports both `amd64` and `arm64`. **Fusion** mounts S3/GS/Azure Blob as POSIX inside containers, skipping stage-in/stage-out.

Both Wave and Fusion are Seqera-developed; Wave is open-source, Fusion is part of the commercial platform.

### 5.6 Reports & observability (30 s)

- **Snakemake:** `--report report.html` (single-file, self-contained, with provenance, rule graph, `report()`-annotated outputs). `benchmark:` directive writes per-job TSV.
- **Nextflow:** `-with-report` (resources), `-with-trace` (per-task TSV), `-with-timeline` (HTML Gantt), `-with-dag` (Graphviz/Mermaid), `-with-weblog` (HTTP POST of run events). **Seqera Platform** (formerly Nextflow Tower) consumes the same telemetry for persistent multi-user views, live monitoring, fleet management. Tower OSS edition exists under MPL-2.0; the full platform is closed-source.

### 5.7 Ecosystem & community (45 s)

- **nf-core** (*Genome Biology* July 2025): 124 curated pipelines, 1400+ modules, 70+ subworkflows, ~2600 contributors, ~10,800 Slack members. All modules carry a Bioconda env, container, `nf-test` test, and metadata. Standardized template + lint/sync CLI.
- **Snakemake ecosystem** is more **federated**: Workflow Catalog is a statically generated index of repos meeting structural criteria; `snakemake-workflows` GitHub org hosts ~28 reference pipelines (RNA-seq, WGS, ChIP-seq, …); **snakemake-wrappers** provides hundreds of versioned per-tool wrappers. No equivalent strict template-conformance regime, no central dashboard with module counts.

Honest summary: **nf-core sets a higher bar for module/pipeline standardization; Snakemake's ecosystem is broader at the wrapper level but less curated at the pipeline level.**

### 5.8 Testing (30 s)

- **Snakemake:** `--generate-unit-tests` scaffolds `pytest` test cases from a real run with intermediates preserved. Byte-exact comparison by default; overridable via an `OutputChecker` per file extension.
- **Nextflow:** **nf-test** (Forer & Schönherr, *GigaScience* 2025) — Groovy `when`/`then` blocks, channel assertions, **snapshot testing**, tag-based test selection. The official testing standard across nf-core modules/pipelines since 2024.

### 5.9 One-page summary slide — when does each shine?

Resist "which is better." Frame as fit:

- **Snakemake fits when:** team is Python-native; pipeline is target-/file-driven; debugging needs interactive Python; you want the DAG up-front; you want one Snakefile to read end-to-end.
- **Nextflow fits when:** team has Java/Groovy/Scala people; pipeline is stream-/event-driven (e.g., real-time ingest, scatter-on-discovery); cloud-native (AWS Batch, Fusion) is on the roadmap; you want curated nf-core pipelines off the shelf.

---

## 6. Section 5 — Live demo on EMBL HPC (5 min)

### 6.1 What to demo (the repo already has the material)

Path: `dev/snakemake-v9-features/` in this repo. It is a minimal 5-rule workflow exercising:

1. The `fs` storage provider (auto stage between shared FS and `/scratch`)
2. Rule grouping (`group: "processing"` keyword)
3. SLURM efficiency reporting
4. Scalability from 3 → 100 samples via config

### 6.2 The 5-minute script

| t | action | what to say |
|---|---|---|
| 0:00 | Open `Snakefile`, scroll the `process_step1` / `process_step2` rules | "Note `group: "processing"` — these two rules will run in the **same** SLURM job." |
| 1:00 | Open `profile/config.yaml` | "Five lines are the entire HPC story: `executor: slurm`, `default-storage-provider: fs`, `local-storage-prefix: /scratch/$USER/...`, `slurm-efficiency-report: true`, `max-jobs-per-second: 10`." |
| 2:00 | `pixi run snakemake --profile profile --dry-run` | "Static DAG — Snakemake shows the full plan before submitting anything." |
| 2:30 | `pixi run snakemake --profile profile -j 10` (kick off real run, leave running) | While it runs, talk over it. |
| 3:00 | Show one job's `/scratch/$USER/snakemake_test/$JOBID/` while it's live | "Files stage in, computation happens here, outputs are copied back." |
| 4:00 | When done: open the efficiency report | "Anything under 80% utilization is flagged — feedback loop for resource right-sizing." |
| 4:30 | `pixi run snakemake --profile profile --report report.html` | "And this is what gets shared with collaborators." |

### 6.3 Demo safety net

The live SLURM demo can fail (queue full, account expired, transient FS issue). Have a screen recording of a successful run ready. If the queue is slow, drop straight to the recording at t=2:30.

### 6.4 What NOT to demo

Don't try to demo the v8 → v9 migration of a real workflow live. Reference our `workflow/snakemake_profiles/.../v7/` directory in this repo as "see the diff in CI: `main.yaml` vs `main_v8.yaml`" — that's enough.

---

## 7. Section 6 — EMBL HPC checklist + gotchas (3 min)

This section condenses `docs/hpc-pipeline-optimization_ongoing.md`. The full draft is 1,241 lines; pick the **5 highest-impact items**, leave the rest as a circulated handout / Notion page link.

### 7.1 The 5-point checklist (one slide)

1. **Throttle scheduler submissions** — `max-jobs-per-second: 10`. Prevents the EMBL SLURM controller from rejecting submissions under load.
2. **Stage compute through scratch** — `default-storage-provider: fs` + `local-storage-prefix: /scratch/$USER/<workflow>`. Forces all intermediate I/O off the parallel FS, where it's a noisy neighbor.
3. **Remove `input-output` from `shared-fs-usage`** — otherwise Snakemake will skip the staging and read/write directly on the shared FS, defeating step 2.
4. **Enable the Snakemake output cache** — `export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/<group>/shared/snakemake_cache` + `cache: True` per rule. Computed indexes (BWA, STAR, reference dictionaries) are reused across users.
5. **Share conda & apptainer caches** — point `conda-prefix` and `apptainer-prefix` to a group-writable parallel-FS directory with the **setgid bit set** (`chmod 2775`). Stops every user re-solving the same env.

### 7.2 Storage tier cheat sheet (slide, leave on screen)

| Data type | Tier |
|---|---|
| Reference genomes | Cached parallel FS (read 100× during alignment) |
| BAM / CRAM intermediates | Parallel FS (high-throughput sequential writes) |
| Computed indexes (BWA, STAR) | Cached parallel FS (repeated reads) |
| Software environments (conda, containers) | Parallel FS (shared across users) |
| Temporary / sorting buffers | Local scratch (job-local, ultra-fast) |

### 7.3 Common pitfalls (verbal — 60 s)

Walk through three quickly:

1. **Stale `.snakemake/` on shared FS** — when two users run the same workflow from the same checkout, they fight over lock files. Fix: each user has their own checkout, or set `--workflow-profile` per user.
2. **Conda env builds racing each other** — if two users run the workflow simultaneously and the env hasn't been solved yet, they'll race. Fix: pre-build envs once with `--conda-create-envs-only`, then use the shared `conda-prefix`.
3. **`$JOBID` not expanded in scratch paths** — bug fixed in Snakemake **9.3.4+**. If you see literal `$JOBID` in your scratch path, upgrade.
4. **`--conda-frontend mamba` in your profile** — harmless warning on v8, removed on v9. Delete it.
5. **Old v7 `--cluster "sbatch ..."` in CI / profile** — silently broken on v8. Replace with `--executor slurm` (or `--executor cluster-generic --cluster-generic-submit-cmd "sbatch ..."` if you genuinely need the old behavior).

### 7.4 Verify-before-you-run handout

Promise to circulate the long-form HPC doc (`docs/hpc-pipeline-optimization_ongoing.md`) after the talk. Keep that on the closing slide.

---

## 8. Section 7 — Q&A (5 min)

### 8.1 Likely questions and prepared answers

**Q: "We just migrated to v8, should we move to v9 now?"**
*A:* Yes, the v8→v9 jump is essentially free — one breaking change in custom logger plugins, which most users don't have. v9 adds rerun-trigger improvements you want anyway.

**Q: "Should I switch from Snakemake to Nextflow / vice versa?"**
*A:* No, not because of this talk. Pick on the fit criteria from slide 5.9. The cost of a migration is in the months, the benefit is rarely structural.

**Q: "Does the live demo workflow work on any SLURM, or only EMBL?"**
*A:* Generic SLURM. EMBL-specific bits are the `/scratch/$USER/...` paths and the storage tier locations — everything else is portable.

**Q: "What about Pixi vs conda?"**
*A:* Pixi is the package manager *invoking* Snakemake. Inside Snakemake, software deployment is still `--sdm conda`. Future direction: a Pixi SDM backend is on the roadmap, not shipped.

**Q: "What about WDL / Cromwell / CWL?"**
*A:* Out of scope; happy to compare offline.

**Q: "Where is the EMBL HPC checklist published?"**
*A:* `docs/hpc-pipeline-optimization_ongoing.md` in this repo; will be circulated.

---

## 9. Appendix A — Sanity-check log on the Snakemake/Nextflow comparison

Per the brief, I double-checked every potentially controversial claim against current docs. Things I corrected vs. common folklore:

- ✅ **Nextflow Java requirement** — *true*, and Java 17+ specifically since Nextflow 25.04 (May 2025).
- ✅ **Nextflow `-resume` is content-hashed** — *incomplete*. Default is **file attributes** (path/size/mtime), not content. `cache 'deep'` opts into true content hashing.
- ✅ **Snakemake can't do cloud** — *false*. Official plugins for Kubernetes, Google Batch, Azure Batch, plus storage plugins for S3/GCS/Azure. Genuinely thinner than Nextflow for AWS Batch (no official plugin).
- ✅ **Snakemake has no resume** — *false*. Rerun triggers are configurable: `--rerun-triggers {mtime,params,input,software-env,code}`.
- ✅ **`--conda-frontend mamba` still needed for speed** — *false*. `libmamba` is the conda default solver; the flag is removed in v9.
- ✅ **Subworkflows are the right way to compose Snakemake workflows** — *false*, deprecated since v6.0. Use `module:` + `use rule * from foo as foo_*`.
- ✅ **Singularity / Apptainer** — Snakemake renamed the flag, the underlying tool is Apptainer (the OSS fork of Singularity).
- ✅ **Wave and Fusion are both commercial** — *partially false*. **Wave is open source** (Apache-2.0, source on GitHub); Fusion is the commercial one.
- ✅ **nf-core has X pipelines** — pinned to the *Genome Biology* July 2025 paper numbers (124 pipelines, 1400+ modules) rather than guessing.

## 10. Appendix B — Material checklist before the talk

- [ ] Confirm slide format with organizer (Marp / Quarto / Google Slides / reveal.js)
- [ ] Record fallback video of the live demo (5 min) and stash on a public link
- [ ] Pre-create the conda envs for the demo: `pixi run snakemake --profile dev/snakemake-v9-features/profile --conda-create-envs-only`
- [ ] Pre-build the apptainer images if any
- [ ] Smoke-test the demo on the actual EMBL SLURM partition you'll demo from, **within the 24 h before the talk** (queue state matters)
- [ ] Push `docs/hpc-pipeline-optimization_ongoing.md` to a circulated location (Notion / mailing list) ahead of the talk so people can read along
- [ ] Have the rolling F1000 paper PDF cached locally in case wifi dies
