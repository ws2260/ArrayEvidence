# ArrayEvidence Architecture

Status: Proposed architecture, accepted design direction
Updated: 2026-10-21

## Mission

ArrayEvidence is a deterministic, privacy-first, assay-aware
genetic evidence pipeline.

Its purpose is to transform array-derived genotype calls into
technically interpretable, reproducible evidence suitable for
authorized downstream genetic analyses.

The pipeline must preserve uncertainty, provenance, assayability,
and the distinction between technical and biomedical evidence.

## Architectural principles

1. Every major module must provide independent scientific utility.
2. Technical genotype processing must not require clinical annotation.
3. Biomedical interpretation must be explicitly optional.
4. Analysis must respect the authorized purpose and data scope.
5. Modules must communicate through documented, versioned interfaces.
6. The original genotype calls must remain unchanged.
7. Missing evidence must never be converted into positive evidence.
8. Decisions must be deterministic and independently testable.
9. Participant-derived data must remain within authorized environments.
10. Each module must preserve provenance and uncertainty.

## Module 1: Array evidence and genotype reconciliation

Responsibilities:

- Manifest parsing and assay identification
- Reference assembly and strand validation
- SNV and indel allele reconciliation
- Multiallelic genotype projection
- Assayability and missingness assessment
- Technical genotype and assay QC
- Evidence provenance and ambiguity reporting

Standalone outputs:

- Structured technical evidence table
- Technical QC report
- Machine-readable provenance

Module 1 establishes technical evidence and does not generate
general-purpose genotype-format exports in v0.1.0.

Any downstream VCF/BCF export must preserve source-call provenance,
represent uncertainty explicitly, and undergo documented,
purpose-appropriate validation.

Module 1 must never require disease annotation.

## Module 2A: Genotype utility and export profiles

Potential genotype-format export capabilities under consideration:

- General-purpose VCF/BCF under explicit export policies
- Single-cell donor demultiplexing
- Preparation of inputs for external genotype imputation
- sc-eQTL and genetic association analysis
- Other authorized population-genetic analyses

Different applications may require different export policies.

Export profiles must not silently change the underlying
technical evidence.

## Module 2B: Variant annotation

Optional biomedical evidence enrichment:

- Variant and gene annotations
- Population allele frequencies
- ClinVar and ClinGen evidence
- Gene-disease validity
- Disease mechanisms and inheritance information

Annotation must respect the approved analysis scope.

Annotations must not overwrite technical genotype evidence.

## Module 3: Biomedical interpretation

Optional, authorization-dependent functions:

- Phenotype-specific variant prioritization
- Inheritance-aware assessment
- Evidence adjudication
- Research-oriented interpretation reports

Technical uncertainty must remain visible.

ArrayEvidence is not a substitute for clinical confirmation
or professional diagnostic interpretation.

## Purpose-scoped execution

The pipeline must support technical-only operation.

Examples:

A. Demultiplexing-only study:
   Module 1 -> Module 2A

B. Phenotype-specific association study:
   Module 1 -> authorized Module 2A analyses

C. Authorized rare-disease study:
   Module 1 -> Module 2B -> Module 3

Module selection is not equivalent to legal or ethics authorization.

Where users are restricted to specific variants or regions,
data-access restrictions must be enforced before or during
processing as appropriate, not merely at the reporting stage.

Unrestricted intermediate genotype files must not be exposed
to users authorized to access only a restricted subset.

## Module interfaces

Each module should declare:

- Required inputs and reference versions
- Permitted data scope
- Output schema and version
- Evidence states and uncertainty semantics
- Provenance
- Validation requirements
- Failure behavior

Allele identity, genotype evidence, technical quality,
clinical annotation, and interpretation must remain
distinct evidence dimensions.

The current Module 1 reconciliation and versioned TSV export schema are documented in [Technical Evidence Contract](EVIDENCE_CONTRACT.md).

## Validation strategy

Synthetic tests:
- Deterministic software correctness
- Edge cases and invariants

AIDA:
- Authorized development and integration testing
- Technical characterization

Held-out GSA-WGS data:
- Independent validation after method freeze

Development-cohort concordance is not an independent
accuracy estimate.

## Publication strategy

Module 1 should be independently releasable and potentially
publishable as a computational genomics methods framework.

Biomedical annotation and interpretation should be optional
extensions, potentially supporting separate publications.

The development priority is a stable technical evidence
platform before broader interpretation capabilities.
