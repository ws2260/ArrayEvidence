# ArrayEvidence

**ArrayEvidence** is a deterministic, privacy-first, assay-aware
pipeline for evaluating genetic evidence from genotyping arrays.

The initial focus is Illumina Global Screening Array (GSA)
data and DRAGEN Array VCF outputs.

## Goals

ArrayEvidence aims to:

1. Reconcile array assay definitions with genotype calls.
2. Validate allele identity against a reference genome.
3. Preserve genotype uncertainty, assayability, and provenance.
4. Produce technically interpretable genotype datasets.
5. Support downstream genetic analyses through independent modules.
6. Enable optional, purpose-scoped biomedical interpretation.

ArrayEvidence starts downstream of genotype calling. It does not
replace DRAGEN Array or other primary genotype callers.

## Modular architecture

### Module 1 — Technical genotype evidence

The standalone core provides:

- Manifest-to-VCF reconciliation
- Reference and strand validation
- SNV and indel allele matching
- Multiallelic genotype projection
- Technical QC and provenance
- Policy-filtered genotype exports

Planned outputs include evidence tables, VCF/BCF, and QC reports.

Module 1 does not require clinical or disease annotation.

### Module 2A — Genotype utility

Optional application-specific exports for:

- Single-cell donor demultiplexing
- Genotype imputation
- sc-eQTL and genetic association analyses

### Module 2B — Variant annotation

Optional, authorization-dependent biomedical annotations.

### Module 3 — Biomedical interpretation

Optional phenotype-specific prioritization, evidence adjudication,
and research-oriented reporting.

## Purpose-limited analysis

Not every study permits broad genetic interpretation.

ArrayEvidence is designed to support technical-only workflows
and explicitly scoped biomedical analyses.

Selecting a module does not itself establish ethics, consent,
or data-use authorization.

## Current development status

The publication-track implementation currently includes:

- Manifest and VCF parsing
- Reference-aware allele validation
- SNV and multiallelic reconciliation
- Multiallelic genotype projection
- Indel event reconstruction and normalization
- Reference-normalized indel candidate discovery
- Allele-aware indel reconciliation and genotype projection

The next milestone is a versioned technical evidence schema
and structured evidence export, developed using test-driven development.

A reconciled genotype is not automatically more accurate than
the original call. Accuracy claims require independent validation.

## Development and validation

- Synthetic fixtures: software correctness and regression tests
- AIDA: authorized development and technical characterization
- Held-out paired GSA-WGS data: planned independent validation
  after method freeze

Participant-derived genetic data must not be committed to
the publication source repository.

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Technical evidence contract](docs/EVIDENCE_CONTRACT.md)
- [Development roadmap](docs/ROADMAP.md)
- [Data policy](docs/DATA_POLICY.md)
- [ADR-0001: Indel candidate discovery](docs/decisions/0001-indel-candidate-discovery.md)
- [ADR-0002: Modular purpose-scoped pipeline](docs/decisions/0002-modular-purpose-scoped-pipeline.md)

## Project status

Research software under active development.

Technical reconciliation is the first independently releasable
component. Biomedical annotation and interpretation are
optional future extensions.
