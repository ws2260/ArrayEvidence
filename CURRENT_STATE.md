# ArrayEvidence — Current State and AI Handoff

Last updated: 2026-10-23
Status: Historical reconstruction; verify live state before execution.
Repository: ~/projects/ArrayEvidence-publication

## 1. Purpose

This document is the operational handoff for the ArrayEvidence
publication repository.

It records verified historical milestones, pending work, decisions,
and the protocol for resuming development.

It is not a substitute for checking Git, files, logs, and the
current execution environment.

## 2. Architecture

The agreed modular architecture is:

- Module 1: Technical genotype evidence, assay-aware allele
  reconciliation, genotype projection, QC, and provenance.
- Module 2A: Purpose-specific genotype utilities and exports,
  including demultiplexing, imputation preparation, and sc-eQTL.
- Module 2B: Optional variant annotation and external evidence.
- Module 3: Optional biomedical interpretation and adjudication.

Module 2B is not post-imputation QC.

### Explicit export-boundary decision (2026-10-23)

Module 1 establishes authoritative technical evidence through
versioned evidence tables, QC, and provenance.

Genotype-format exports, including generic VCF/BCF and
application-specific imputation, demultiplexing, and sc-eQTL
inputs, belong to Module 2A rather than Module 1.

Module 2A is under consideration, not implemented or committed
to a release. A limited Qdenga imputation-preparation workflow
is a feasibility objective, not an established capability.

Preparing genotype inputs for an external imputation workflow
does not mean ArrayEvidence performs genotype imputation.

Module 2B and Module 3 are also under consideration.

Documentation must distinguish available functionality,
active development, and capabilities under consideration.
An architectural proposal is not a delivery commitment.

Authoritative architecture:
docs/ARCHITECTURE.md
docs/decisions/0002-modular-purpose-scoped-pipeline.md

## 3. Release history and plans

### v0.1.0 — Released frozen Module 1 baseline

Historical evidence:
- Release tag v0.1.0 points to commit dda3e5e.
- 402 automated tests passed, 2 skipped (verified 2026-10-23).
- AIDA: 48 sample evidence bundles built, verified, and QC audited.
- Qdenga UAT-001: 21/21 independent DRAGEN Array samples passed
  build, verify, and qc.
- Clean-install reproducibility confirmed on one Qdenga sample:
  byte-identical evidence tables and identical aggregate QC.

These results establish technical reproducibility, not
independent genotype accuracy.

### v0.1.1 — Proposed maintenance release

Candidate improvements:
- Structured QC JSON output.
- Installation and quick-start documentation.
- Minor usability and error-handling fixes.
- Performance and resource benchmarking.

Preserve backward compatibility and scientific behavior.

### v0.2.0 — Proposed feature release

Candidate improvements:
- Batch orchestration.
- Aggregate cohort QC reporting.
- Other substantial functionality extensions.

No agreed v0.1.2 scope has been recovered.

Do not treat proposals as implemented features.

## 4. Grant progress report priorities

Immediate priority: reproducible technical and scientific evidence
for the dengue sc-eQTL grant progress report.

Completed:
- AIDA Module 1 technical validation.
- Frozen v0.1.0 release.
- Independent Qdenga Module 1 UAT-001.

Pending or requiring fresh verification:
- Finalize and commit public validation documentation.
- Qdenga imputation preparation and GeTH submission.
- Limited Module 2A proof of concept, if feasible.
- Qdenga WGS cohort QC and variant calling.
- WGS-array genotype concordance, if comparable calls become available.
- Grant report figures and accurate validation statements.

## 5. Qdenga WGS status

Historical checkpoint:
- 21 indexed BAMs on MU-AI HPC.
- Sequential Parabricks HaplotypeCaller job: 369139.
- Planned gVCF directory:
  ~/NL/Qdenga/gVCF_ALT/

Job completion and resulting files have NOT been verified
in this handoff.

Do not claim WGS-array concordance until appropriate
WGS genotype calls and matching QC have been assessed.

## 6. Documentation checkpoint

Documentation review completed through Step 289 on 2026-10-23.

Reviewed the Module 1/2A export boundary, release claims,
validation documentation, and AI continuity files.
The standard pytest run returned 402 passed, 2 skipped.

Eight documentation files are modified or untracked.
They remain unstaged and uncommitted pending final review.
Do not discard or commit them without explicit approval.

## 7. Scientific guardrails

- Preserve participant-level genetic data outside public Git.
- Do not publish proprietary or restricted vendor resources.
- Do not conflate assay reconciliation with genotype accuracy.
- Do not interpret technical MISMATCH or PARTIAL_MATCH
  classifications as independently confirmed genotype errors.
- Do not change reconciliation algorithms merely to increase
  MATCH percentages.
- Require independent evidence, regression tests, and
  documented scientific decisions for semantic changes.
- Keep Module 1 and Module 2A responsibilities separate.
- Preserve the frozen v0.1.0 release.

## 8. AI continuity protocol

At the beginning of a new development session:

1. Read CURRENT_STATE.md.
2. Read docs/ARCHITECTURE.md and docs/ROADMAP.md.
3. Read applicable architecture decisions and validation documents.
4. Inspect Git status, relevant commits, and the environment.
5. Distinguish verified facts from historical reports and plans.
6. Identify the last completed step before proposing another.
7. Proceed one command at a time and wait for user output.
8. Never repeat completed operations without a reason.

If this document is missing, inaccessible, outdated, or
contradicted by other evidence, explicitly ask the user
for the latest CURRENT_STATE.md and relevant logs.

Never pretend to remember missing project details.

After significant milestones, update this document with:
- Date and verified repository state.
- Commands/tests executed and their actual outcomes.
- Files changed.
- Decisions and unresolved issues.
- Exact next action.

## 9. Immediate next action

Documentation consistency, release provenance, and the
tracked-file diff have been reviewed. The next action is
to obtain explicit approval before staging or committing.
After documentation is finalized, return to the dengue
sc-eQTL grant progress report priorities.

Do not automatically commit any files.
