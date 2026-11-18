# MEANtools — Software Maintenance & Modernization

> **Project purpose:** Improve the reliability, compatibility, reproducibility, and maintainability of MEANtools through systematic bug fixes, dependency updates, Python runtime modernization, automated testing, and clearer documentation—while preserving scientifically meaningful results.

**Project type:** Software maintenance and engineering improvement  
**Primary codebase:** [MEANtools development repository](https://github.com/kumarsaurabh20/Meantools_dev)  
**Reference documentation:** [MEANtools documentation](https://meantools.readthedocs.io/en/latest/)  
**Scope:** Existing software functionality and infrastructure, not a redesign of the scientific methodology  
**Status:** Ongoing; individual tasks and milestones are tracked through GitHub Issues and Projects

## 1. Background and motivation

MEANtools is a Python-based multi-omics analysis framework for connecting transcriptomic, metabolomic, and genomic evidence to the investigation of plant-specialized metabolism and biosynthetic pathways. The software includes workflows for gene–metabolite correlation, functional-cluster analysis, mass-transition processing, database preparation, reaction and pathway inference, and downstream visualization.

As research software evolves, small bugs accumulate, dependencies change, interfaces become outdated, and operating environments diverge. These changes can affect installation, runtime stability, reproducibility, and the ability of other researchers to reuse the software. In computational biology, an apparently minor change can also affect scientific interpretation if it alters filtering, ordering, missing-value handling, numerical precision, or output formats.

This project establishes a coordinated maintenance effort. Rather than introducing major new biological inference methods, it strengthens the existing MEANtools foundation so that current analyses are dependable and future methodological development is easier to implement and validate.

## 2. Project objectives

1. **Fix reproducible bugs:** Identify, prioritize, reproduce, and resolve defects in existing workflows.
2. **Modernize dependencies:** Review Python packages, remove unnecessary dependencies, replace deprecated APIs, and adopt supported compatible versions.
3. **Improve runtime compatibility:** Document and test supported Python versions and the environments used locally and on HPC systems.
4. **Establish automatic quality checks:** Use GitHub Actions for syntax validation, linting, unit tests, and reproducible integration tests.
5. **Protect scientific correctness:** Introduce regression tests so maintenance work does not silently change valid biological or numerical results.
6. **Improve maintainability:** Reduce duplicated logic, improve error messages, rationalize configuration, and strengthen modularity without unnecessary large rewrites.
7. **Improve reproducibility and onboarding:** Make installation, test data, database acquisition, and common commands easier to understand.
8. **Prepare stable maintenance releases:** Track changes, compatibility decisions, validation evidence, and known limitations.

## 3. Scope of work

### 3.1 Bug fixes and robustness

Examples of tasks include:

- Fix crashes caused by unexpected empty inputs, malformed files, missing columns, or inconsistent identifiers.
- Improve argument validation and actionable error messages.
- Handle missing data, duplicate identifiers, and unexpected column ordering explicitly.
- Audit path handling, temporary directories, output naming, and file overwriting.
- Fix incorrect assumptions about input types or file extensions.
- Address platform-specific issues that arise on macOS, Linux, or supported HPC environments.
- Confirm that error handling does not hide scientific or data-quality problems.

Each issue should include reproducible steps, observed behavior, expected behavior, relevant input characteristics, and a small test case where possible.

### 3.2 Dependency and packaging modernization

- Inventory dependencies declared in `pyproject.toml`, `setup.py`, `environment.yml`, requirements files, and scripts where present.
- Identify duplicated, obsolete, unpinned, incompatible, or unmaintained requirements.
- Replace deprecated APIs where behavior can be preserved and verified.
- Define a **supported version range**, rather than upgrading every package blindly to the newest release.
- Test package resolution in a clean environment, documenting any platform-specific limitations.
- Clarify which packages are runtime, optional, development, or HPC-only dependencies.
- Evaluate consolidation of packaging metadata to avoid contradictory dependency definitions.
- Add automated dependency review or update proposals only after the baseline test suite is reliable.

**Policy:** A dependency upgrade is complete only when the affected features have been tested and any behavior-changing differences are documented.

### 3.3 Python and operating-environment compatibility

- Select a primary supported Python version and one or more additional supported versions based on tested dependencies.
- Record operating-system and architecture assumptions, including supported macOS/Linux environments where applicable.
- Avoid implicit reliance on developer-specific local paths or manually configured environment variables.
- Make command-line entry points and configuration options predictable.
- Document environment creation and installation commands after they have been tested on a clean machine.
- Distinguish errors attributable to MEANtools from errors attributable to missing external tools or databases.

The supported-version matrix is a **decision to be established by testing**, not a promise that a particular Python version already works.

### 3.4 Routine correlation workflow

Areas to review include:

- Input transcriptomics/metabolomics matrix orientation and sample matching.
- Correlation calculation, thresholding, ranking, and filtering behavior.
- Numerical handling of ties, constants, missing values, and very small sample sizes.
- Deterministic output ordering and consistent output schemas.
- Memory efficiency when working with large paired-omics datasets.
- Clear documentation of expected inputs, options, outputs, and assumptions.

**Initial target:** `corrMultiomics.py` (verify its path in the development repository).

### 3.5 HPC correlation workflow

Areas to review include:

- Compatibility between routine and HPC correlation implementations.
- Correct job configuration, parameter propagation, and output directory handling.
- Consistent results when analyses are partitioned, parallelized, or recombined.
- Reproducible behavior under different worker counts and task scheduling arrangements.
- Informative logging, error propagation, and failure recovery.
- Separation of local syntax/unit checks from tests that require an actual scheduler or computing cluster.

**Initial targets:** `HPC_test/corrMultiomics_HPC.py` and `HPC_test/corrMultiomics.sh`.

**Important:** Passing `bash -n` validates shell syntax, **not** SLURM directive correctness, cluster compatibility, or successful HPC execution.

### 3.6 Other existing MEANtools modules

Maintenance should be extended incrementally to the current modules for database formatting, functional-cluster merging, mass-transition processing, structure mapping, reaction scoring, pathway inference, and graphics. Prioritize components that are frequently used, generate known defects, or block reproducible installations.

The project does **not** assume that every module must be rewritten. Targeted, verified improvements are preferred.

### 3.7 Documentation and developer experience

- Maintain clear setup instructions, dependency notes, usage examples, and troubleshooting guidance.
- Document required reference databases separately from the small test datasets included in Git.
- Describe expected file formats and common validation errors.
- Provide guidance for running routine versus HPC analyses.
- Keep developer documentation, CI instructions, and release notes up to date.
- Use short, reproducible examples with expected outputs wherever feasible.

## 4. Out of scope

Unless approved as a separate issue or milestone, this maintenance project does not cover:

- Major redesign of reaction-scoring or global pathway-inference algorithms.
- New MS/MS, single-cell, or spatial-omics methodological research modules.
- New large-scale experimental data generation.
- Uploading full reference databases or confidential datasets to GitHub.
- Large code rewrites without a concrete maintenance or correctness benefit.

Feature development may continue in parallel, but should be tracked in dedicated issues or projects and integrated only when tests and compatibility requirements are met.

## 5. Workstreams and deliverables

| Workstream | Primary focus | Evidence of completion |
| --- | --- | --- |
| Bugs and reliability | Reproducible defects, safe failure modes | Bug-reproducing test passes and issue is resolved |
| Dependencies | Compatible supported packages and packaging | Clean installation and documented dependency changes |
| Runtime compatibility | Python versions and operating environments | Passing compatibility matrix and published support policy |
| Correlation workflows | Routine/HPC consistency | Synthetic and representative regression tests pass |
| HPC integration | Submission scripts and scheduler assumptions | Local shell validation plus cluster-level validation when available |
| Testing and CI | Automated quality gates | Workflow results, test reports, and reviewable failures |
| Documentation | Clear setup and troubleshooting | Reproducible instructions verified by another user/environment |
| Releases | Stable maintenance update | Changelog, version/tag, smoke tests, and documented limitations |

## 6. GitHub workflow and project management

### Recommended project fields

- **Status:** Backlog → Ready → In Progress → In Review → Done
- **Type:** Bug / Dependency / Compatibility / Refactor / Test / Documentation / Performance
- **Priority:** P0 Critical / P1 High / P2 Medium / P3 Low
- **Component:** Correlation / HPC / Databases / Pathways / Packaging / CI / Docs
- **Assignee:** Responsible contributor
- **Milestone:** Baseline, Compatibility, Regression Testing, Maintenance Release

### Issue expectations

Every nontrivial task should have a GitHub Issue describing:

1. **Problem:** What is broken, outdated, or difficult to maintain?
2. **Impact:** Which commands, users, platforms, or scientific results are affected?
3. **Reproduction:** How can the problem be demonstrated with minimal data?
4. **Proposed change:** What is the smallest reasonable fix?
5. **Acceptance criteria:** How will we know the issue is resolved?
6. **Validation:** Which automated tests and manual checks are required?

Use labels consistently so issues can be filtered by component and priority. Link each Pull Request to the corresponding Issue where possible.

### Branch and review conventions

Use descriptive branch names such as `fix/correlation-missing-values`, `chore/update-dependencies`, `test/hpc-correlation-parity`, or `docs/install-guide`.

A typical contribution is: create an Issue → create a branch → implement the change → add/update tests → open a Pull Request → review CI results and scientific impact → merge to `main` → close the Issue.

Avoid committing substantial unreviewed changes directly to `main` once automated checks are established.

## 7. Testing and continuous integration

Testing will be introduced in layers, starting with fast checks and advancing toward scientific validation.

### Level 1 — Syntax and repository hygiene

- Compile Python sources to catch syntax errors.
- Parse shell scripts with `bash -n`.
- Check formatting/linting using an agreed tool such as Ruff.
- Exclude `__pycache__/`, `*.pyc`, temporary outputs, and large reference data from Git.

### Level 2 — Unit tests

Create lightweight `pytest` tests for functions that handle data parsing, correlation calculations, thresholds, normalization assumptions, error conditions, and transformations. Use tiny synthetic fixtures committed to the repository.

### Level 3 — Integration tests

Exercise small end-to-end workflows, checking expected output files, schemas, and numerical properties. Integration tests should not require hundreds of megabytes of reference data unless those assets are downloaded explicitly in a controlled environment.

### Level 4 — Scientific regression tests

Maintain reference results for approved examples and compare them after changes. Tests should identify:

- Numerical differences beyond defined tolerances.
- Changes to feature pairing, filtering, or ranking.
- Missing or duplicated outputs.
- Unexpected nondeterminism.
- Differences between routine and HPC implementations for equivalent inputs.

Scientific regression tests must accommodate documented floating-point tolerances and explicitly identify results where randomization or parallel scheduling affects ordering.

### Level 5 — HPC validation

On an authorized HPC system, validate scheduler submissions, environment activation, memory/runtime behavior, parallel consistency, log files, and output recovery. These are **separate from ordinary GitHub-hosted CI** unless a secure, approved HPC integration is configured.

### Proposed GitHub Actions workflows

Start with `.github/workflows/meantools-ci.yml` and add checks as the test suite grows. One YAML workflow may contain multiple steps or jobs. Separate workflows should be introduced when triggers, compute requirements, security, or release processes differ.

Recommended triggers:

- On pull requests: syntax, lint, unit tests, and lightweight integration checks.
- On pushes to `main`: the same essential checks.
- On manual dispatch or scheduled runs: optional compatibility or dependency tests.
- Before a release: full available regression and packaging checks.

**Passing CI is necessary but not sufficient for scientific correctness.** Acceptance also depends on documented numerical and biological validation where relevant.

## 8. Data management and repository hygiene

MEANtools relies on sizeable reference resources and may generate large intermediate files. The repository should remain source-code oriented.

**Commit to Git:** source code, scripts, configuration, documentation, small synthetic examples, stable metadata schemas, tests, and reasonably sized expected outputs.

**Keep outside Git:** full reaction/chemical reference databases, large omics matrices, HPC job output, scratch files, generated caches, and any sensitive or restricted data.

For externally stored reference datasets:

- Maintain a documented retrieval source and version where available.
- Provide checksums for reproducible downloads when possible.
- Record data licenses and appropriate citations.
- Separate test fixtures from production reference data.
- Avoid hard-coded developer-machine paths.
- Never add credentials, tokens, SSH keys, or sensitive research data to the repository.

Ignoring a file with `.gitignore` does not remove it from existing Git history; history cleanup must be managed separately if oversized files were previously committed.

## 9. Dependency upgrade policy

Changes to dependencies should follow an incremental, evidence-based process:

1. Record the currently working environment and baseline outputs.
2. Identify the dependency to update, the reason, and any upstream breaking changes.
3. Update a small related group of packages rather than everything simultaneously.
4. Run installation, syntax, unit, integration, and relevant regression tests.
5. Document changed behavior, minimum supported versions, and known platform constraints.
6. Merge only after the tests and scientific checks pass or exceptions are explicitly approved.

Prefer actively maintained dependencies, but avoid upgrades that destabilize validated calculations without a clear migration plan. Version constraints and any lock or environment files should reflect the project's actual supported configurations.

## 10. Proposed roadmap

### Phase A — Stabilize the baseline

- [ ] Inventory existing scripts, environments, and known issues.
- [ ] Record current installation and correlation-workflow behavior.
- [ ] Confirm `.gitignore` and remove generated files from active tracking.
- [ ] Add Python/Bash syntax checks for routine and HPC correlation scripts.
- [ ] Create small test fixtures and capture representative baseline outputs.

**Exit criteria:** A reproducible baseline with a documented list of known defects and CI checks running on new changes.

### Phase B — Modernize dependencies and runtime support

- [ ] Audit all dependency declarations and remove inconsistencies.
- [ ] Identify supported Python versions based on tested compatibility.
- [ ] Replace deprecated APIs that block supported environments.
- [ ] Test clean installation on the primary development environment and a Linux environment.
- [ ] Document version constraints and troubleshooting steps.

**Exit criteria:** A verified, documented environment that installs without undocumented manual workarounds.

### Phase C — Add regression testing and fix priority defects

- [ ] Introduce `pytest` for critical functions.
- [ ] Add edge-case handling and tests for correlation inputs.
- [ ] Compare routine and HPC outputs on a common small dataset.
- [ ] Fix high-priority defects, each with a regression test.
- [ ] Validate representative outputs from downstream modules after changes.

**Exit criteria:** Priority defects resolved and core numerical behavior protected by tests.

### Phase D — Prepare and publish a maintenance release

- [ ] Finalize user and developer documentation.
- [ ] Run clean-environment installation and end-to-end smoke tests.
- [ ] Complete available HPC validation and document remaining limitations.
- [ ] Review breaking changes, version number, and migration notes.
- [ ] Publish a changelog and tagged release when approved.

**Exit criteria:** A maintenance release with traceable changes, reproducibility evidence, and clear limitations.

*The phases describe a proposed sequence, not completed work or fixed delivery dates.*

## 11. Definition of done

A maintenance issue is **Done** when all applicable conditions are satisfied:

- [ ] The original defect or limitation is clearly described.
- [ ] The implementation fixes the issue without unrelated changes.
- [ ] Relevant tests are added or updated and pass locally and in CI.
- [ ] Scientific outputs are unchanged or intentional changes are explained and validated.
- [ ] Dependencies and compatibility requirements are updated where relevant.
- [ ] Documentation and examples reflect the new behavior.
- [ ] No large datasets, generated caches, or secrets are added to Git.
- [ ] The Pull Request has been reviewed and linked to the issue.
- [ ] Any remaining limitations or follow-up tasks have been recorded.

For changes that cannot be validated automatically—for example, cluster-specific runs—the missing verification must be stated explicitly rather than represented as complete.

## 12. Prioritization and risk management

| Priority | Typical examples | Response approach |
| --- | --- | --- |
| P0 — Critical | Scientifically incorrect outputs, data loss, security exposure, widely broken core workflow | Stop affected releases; investigate immediately |
| P1 — High | Reproducible crash, supported installation broken, core HPC workflow failure | Address in the next maintenance cycle |
| P2 — Medium | Recoverable error, confusing validation, partial compatibility defect | Schedule after core blockers |
| P3 — Low | Cosmetic improvement, minor refactor, low-impact documentation issue | Tackle opportunistically |

Key risks include scientific output drift after package upgrades; hidden assumptions in large database formats; discrepancies between local and HPC results; dependency conflicts; and over-expanding maintenance into research feature development. Mitigate these by maintaining baseline outputs, using small targeted changes, recording environment versions, and reviewing scientific differences explicitly.

## 13. Measuring progress

The project should track a small set of interpretable indicators:

- Open and resolved issues by priority and component.
- Percentage of priority bugs with a reproducible regression test.
- CI pass rate on the default branch and pull requests.
- Number of tested Python/environment combinations.
- Installation success in a clean environment.
- Routine/HPC parity for defined benchmark inputs.
- Number of undocumented manual setup steps removed.
- Release readiness: critical issues, tests, docs, and validation evidence.

Avoid equating line coverage alone with scientific confidence; a smaller set of meaningful regression tests can be more valuable than broad but superficial coverage.

## 14. Contributing

Contributions are welcome in the form of well-scoped bug reports, reproducible test cases, documentation improvements, dependency fixes, and code changes.

Before substantial work, open or claim an Issue to avoid duplicated effort. Prefer small Pull Requests that are easy to test and review. For scientific calculations, explain whether results should remain identical and attach a small comparison example. Do not submit large datasets or proprietary inputs; construct synthetic reproductions when possible.

## 15. References and related resources

- [MEANtools development repository](https://github.com/kumarsaurabh20/Meantools_dev)
- [MEANtools public repository](https://github.com/kumarsaurabh20/meantools)
- [MEANtools documentation](https://meantools.readthedocs.io/en/latest/)
- [Published MEANtools research article](https://doi.org/10.1371/journal.pbio.3003307)
- [GitHub Actions documentation](https://docs.github.com/en/actions)
- [GitHub Projects documentation](https://docs.github.com/en/issues/planning-and-tracking-with-projects)

---

**Guiding principle:** Modernize the implementation without compromising the interpretability, reproducibility, or scientific validity of MEANtools results.
